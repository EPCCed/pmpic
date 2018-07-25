!Module containing routines that act on semi-spectral data
! The data is spectral in the x and y directions, and is stored as real arrays, but with alternating
! elements corresponding to the real and imaginary parts of complex numbers
! So for an array x(2n), then x(1) is the real component of a complex number, and x(2) is the
! imaginary component et cetera.
!
! Current routines:
! - diffx (differentiate wrt x using wavenumber multiplication)
! - diffy (differentiate wrt x using wavenumber multiplication)
!
! Planned Routines:
! - tridiag - solve a tridiagonal series of equations
! - diffz - differentiate wrt z using tridiagonal
! - laplinv - invert laplacian using tridiagonal method
module fftops_mod
  use fftw_mod
  use state_mod
  use datadefn_mod, only: DEFAULT_PRECISION, PRECISION_TYPE
  use MPI


implicit none

real(kind=DEFAULT_PRECISION), parameter :: twopi=8.d0*atan(1.d0)

real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:) :: kx, ky, kz, k2, filter

real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:) :: left_sendbuff, right_sendbuff, &
                                                             down_sendbuff, up_sendbuff, &
                                                             left_recvbuff, right_recvbuff, &
                                                             down_recvbuff, up_recvbuff

integer :: x_start, y_start, x_stop, y_stop
integer :: left, right, up, down
integer :: nx, ny, nz
logical :: initialised = .false.
integer :: ierr

logical :: x_start_swap, x_end_swap, y_start_swap, y_end_swap !do we need to swap the first/last array elements with neighbouring processes

contains

!initialises data structures needed, communication buffers etc
  subroutine fftops_init(state,xs, ys, sizes)
    type(model_state_type), intent(inout) :: state
    integer, intent(in) :: xs, ys, sizes(:)
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: xval, yval

    if (initialised) then
      print *, "Error fft_ops is already initialised"
      call MPI_Finalize(ierr)
      stop
    endif

    ! set array size parameters and allocate k arrays
    x_start=xs
    y_start=ys
    nx=sizes(3)
    ny=sizes(2)
    nz=sizes(1)
    x_stop = x_start+nx-1
    y_stop = y_start+ny-1

    allocate(kx(nz,ny,nx), ky(nz,ny,nx), k2(nz,ny,nx))

    !determine structure of array (i.e. are all the complex numbers on the array, or are some truncated)

    if (mod(x_start,2) .eq. 0) then
      x_start_swap=.true.
    else
      x_start_swap=.false.
    endif

    if (mod(x_stop,2) .eq. 0) then
      x_end_swap=.false.
    else
      x_end_swap=.true.
    endif

    if (mod(y_start,2) .eq. 0) then
      y_start_swap=.true.
    else
      y_start_swap=.false.
    endif

    if (mod(y_stop,2) .eq. 0) then
      y_end_swap=.false.
    else
      y_end_swap=.true.
    endif

    ! allocate send/recv buffers if needed

    if (x_start_swap) then
      allocate(left_sendbuff(nz,ny), left_recvbuff(nz,ny))
    endif
    if (x_end_swap) then
      allocate(right_sendbuff(nz,ny),right_recvbuff(nz,ny))
    endif
    if (y_start_swap) then
      allocate(down_sendbuff(nz,nx),down_recvbuff(nz,nx))
    endif
    if (y_end_swap) then
      allocate(up_sendbuff(nz,nx),up_recvbuff(nz,nx))
    endif


    !set up wavenumber arrays

    do i=1,nx
      xval = ((x_start+i-2)/2) /(state%global_grid%resolution(3)*state%global_grid%size(3))
      do j=1,ny
        yval = ((y_start+j-2)/2) / (state%global_grid%resolution(2)*state%global_grid%size(2))
        do k=1,nz
          kx(k,j,i) = xval
          ky(k,j,i) = yval
          k2(k,j,i) = xval*xval + yval*yval
        enddo
      enddo
    enddo



    !determine neighbours

    down=state%local_grid%neighbours(2,1) ! in -y direction
    up=state%local_grid%neighbours(2,3) !in +y direction
    left=state%local_grid%neighbours(3,1) ! in -x direction
    right=state%local_grid%neighbours(3,3) ! in +x direction


    initialised = .true.


  end subroutine

  !gives spectral derivative in the x direction: out = i*kx*in
  !This requires us to swap the real and imaginary parts of the numbers around - essentially swap
  !even and odd array elements around. The catch is that these arrays are decomposed between
  !processes and some complex number pairs may be split between processes, so we sometimes need to
  !send messages between processes to swap these pairs round
  subroutine diffx(state,in,out)
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), intent(in) :: in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(out) :: out(:,:,:)
    integer :: left_sendrequest, left_recvrequest, right_sendrequest, right_recvrequest
    integer :: istart, iend
    integer :: i, j, k
    integer :: statuses(MPI_STATUS_SIZE,4)
    integer :: requests(4)=MPI_REQUEST_NULL

    istart=1
    iend=nx

    print *, "enter diffx"

    !send/recv start/end values of arrays if needed (non-blocking)
    if (x_start_swap) then
      left_sendbuff(:,:) = in(:,:,1)
      call MPI_isend(buf=left_sendbuff,&
                     count=nz*ny,&
                     datatype=PRECISION_TYPE,&
                     dest=left,&
                     tag=0,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(1),&
                     ierror=ierr)
      call MPI_irecv(buf=left_recvbuff,&
                     count=nz*ny,&
                     datatype=PRECISION_TYPE,&
                     source=left,&
                     tag=1,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(2),&
                     ierror=ierr)
      istart=2
    endif

    if (x_end_swap) then
      right_sendbuff(:,:) = in(:,:,nx)
      call MPI_isend(buf=right_sendbuff,&
                     count=nz*ny,&
                     datatype=PRECISION_TYPE,&
                     dest=right,&
                     tag=1,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(3),&
                     ierror=ierr)
      call MPI_irecv(buf=right_recvbuff,&
                     count=nz*ny,&
                     datatype=PRECISION_TYPE,&
                     source=right,&
                     tag=0,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(4),&
                     ierror=ierr)
      iend=nx-1
    endif

    !swap internal values (i.e. multiply by i=sqrt(-1))
    do i=istart,iend,2
      out(:,:,i)=-1.*in(:,:,i+1)
      out(:,:,i+1) = in(:,:,i)
    enddo

    !wait for comms to complete (if necessary)

      call MPI_Waitall(count=4,&
                       array_of_requests=requests,&
                       array_of_statuses=statuses,&
                       ierror=ierr)


    !set end values (if necessary)
    if (x_start_swap) then
      out(:,:,1) = left_recvbuff(:,:)
    endif
    if (x_end_swap) then
      out(:,:,nx) = -1.*right_recvbuff(:,:)
    endif

    !multiply by k

    out(:,:,:) = out(:,:,:) * kx(:,:,:)*twopi

  end subroutine

  !gives spectral derivative in the y direction: out = i*ky*in
  subroutine diffy(state,in,out)
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), intent(in) :: in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(out) :: out(:,:,:)
    integer :: jstart, jend
    integer :: i, j, k
    integer :: statuses(MPI_STATUS_SIZE,4)
    integer :: requests(4)=MPI_REQUEST_NULL

    jstart=1
    jend=nx

    print *, "enter diffy"

    !send/recv start/end values of arrays if needed (non-blocking)
    if (y_start_swap) then
      down_sendbuff(:,:) = in(:,1,:)
      call MPI_isend(buf=down_sendbuff,&
                     count=nz*nx,&
                     datatype=PRECISION_TYPE,&
                     dest=down,&
                     tag=0,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(1),&
                     ierror=ierr)
      call MPI_irecv(buf=down_recvbuff,&
                     count=nz*nx,&
                     datatype=PRECISION_TYPE,&
                     source=down,&
                     tag=1,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(2),&
                     ierror=ierr)
      jstart=2
    endif

    if (y_end_swap) then
      up_sendbuff(:,:) = in(:,ny,:)
      call MPI_isend(buf=up_sendbuff,&
                     count=nz*nx,&
                     datatype=PRECISION_TYPE,&
                     dest=up,&
                     tag=1,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(3),&
                     ierror=ierr)
      call MPI_irecv(buf=up_recvbuff,&
                     count=nz*nx,&
                     datatype=PRECISION_TYPE,&
                     source=up,&
                     tag=0,&
                     comm=state%parallel%monc_communicator,&
                     request=requests(4),&
                     ierror=ierr)
      jend=ny-1
    endif

    !swap internal values (i.e. multiply by i=sqrt(-1))
    do j=jstart,jend,2
      out(:,j,:)=-1.*in(:,j+1,:)
      out(:,j+1,:) = in(:,j,:)
    enddo

    !wait for comms to complete

      call MPI_Waitall(count=4,&
                       array_of_requests=requests,&
                       array_of_statuses=statuses,&
                       ierror=ierr)


    !set end values (if necessary)
    if (y_start_swap) then
      out(:,1,:) = down_recvbuff(:,:)
    endif
    if (y_end_swap) then
      out(:,ny,:) = -1.*up_recvbuff(:,:)
    endif

    !multiply by k

    out(:,:,:) = out(:,:,:) * ky(:,:,:) *twopi

  end subroutine
end module
