!Module containing routines that act on semi-spectral data
! The data is spectral in the x and y directions, and is stored as real arrays, but with alternating
! elements corresponding to the real and imaginary parts of complex numbers
! So for an array x(2n), then x(1) is the real component of a complex number, and x(2) is the
! imaginary component et cetera.
!
! Current routines:
! - diffx (differentiate wrt x using wavenumber multiplication)
! - diffy (differentiate wrt x using wavenumber multiplication)
! - tridiagonal - solve set of tridiagonal equations
! - diffz - differentiate wrt z using tridiagonal method
! - laplinv - invert laplacian using tridiagonal method
! - spectral_filter - filter spectral data according to the 2/3 rule
module fftops_mod
  use state_mod
  use datadefn_mod, only: DEFAULT_PRECISION, PRECISION_TYPE
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use MPI


implicit none

real(kind=DEFAULT_PRECISION), parameter :: twopi=8.d0*atan(1.d0)
real(kind=DEFAULT_PRECISION) :: dz

real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:) :: kx, ky, kz, k2, filter

real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:) :: left_sendbuff, right_sendbuff, &
                                                             down_sendbuff, up_sendbuff, &
                                                             left_recvbuff, right_recvbuff, &
                                                             down_recvbuff, up_recvbuff

integer :: x_start, y_start, x_stop, y_stop
integer :: left, right, up, down
integer :: nx, ny, nz
logical :: initialised = .false.
integer :: ierr, comm

integer :: diffxyhandle, diffzhandle, laplinvhandle, filterhandle

logical :: x_start_swap, x_end_swap, y_start_swap, y_end_swap !do we need to swap the first/last array elements with neighbouring processes



contains

!initialises data structures needed, communication buffers etc
  subroutine fftops_init(state,xs, ys, sizes)
    type(model_state_type), intent(inout) :: state
    integer, intent(in) :: xs, ys, sizes(:)
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: xval, yval, ksqmax, lx, ly

    if (initialised) then
      if (state%parallel%my_rank .eq. 0) print *, "Warning: fft_ops is already initialised"
      return
    endif

    ! set array size parameters and allocate k arrays
    x_start=xs
    y_start=ys
    nx=sizes(3)
    ny=sizes(2)
    nz=sizes(1)
    x_stop = x_start+nx-1
    y_stop = y_start+ny-1

    allocate(kx(nz,ny,nx), ky(nz,ny,nx), k2(nz,ny,nx), filter(nz,ny,nx))

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

   !set up wavenumber arrays K = n/l
    lx = state%global_grid%resolution(3)*state%global_grid%size(3)
    ly = state%global_grid%resolution(2)*state%global_grid%size(2)

    ksqmax = (state%global_grid%size(3)/2./lx)**2 + (state%global_grid%size(2)/2./ly)**2
    ksqmax = ksqmax * 2./9. !set limit to (2/3)^2 * kmax^2
    do i=1,nx
      xval = ((x_start+i-2)/2) /(lx)
      do j=1,ny
        yval = ((y_start+j-2)/2) / (ly)
        do k=1,nz
          kx(k,j,i) = xval
          ky(k,j,i) = yval
          k2(k,j,i) = xval*xval + yval*yval
        enddo
      enddo
    enddo

    !assign filter
    where( k2 .gt. ksqmax)
      filter = 0.
    elsewhere
      filter = 1.
    endwhere

    !needed normalisation
    kx(:,:,:) = kx(:,:,:)*twopi
    ky(:,:,:) = ky(:,:,:)*twopi
    k2(:,:,:) = k2(:,:,:)*twopi*twopi


    !determine neighbours

    down=state%local_grid%neighbours(2,1) ! in -y direction
    up=state%local_grid%neighbours(2,3) !in +y direction
    left=state%local_grid%neighbours(3,1) ! in -x direction
    right=state%local_grid%neighbours(3,3) ! in +x direction

    !cache monc communicator
    comm=state%parallel%monc_communicator



    !cache dz
    dz=state%global_grid%resolution(1)

    call register_routine_for_timing("diffxy", diffxyhandle, state)
    call register_routine_for_timing("diffz", diffzhandle, state)
    call register_routine_for_timing("laplinv", laplinvhandle, state)
    call register_routine_for_timing("spec_filter", filterhandle, state)




    initialised = .true.


  end subroutine

  !spectrally filters f according to the 2/3 rule:
  ! -wavenumbers greater than 2/3 of the maximum wavenumber are set to zero
  ! Arguments:
  ! - f               = array to be filtered
  ! - out (OPTIONAL)  = array to have filtered f returned to. If this input is not
  !                     present then the filter is done in-place
  subroutine spectral_filter(f,out)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: f
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out), optional :: out
    !$OMP SINGLE
    call timer_start(filterhandle)
    !$OMP END SINGLE

    if (present(out)) then
      !$OMP WORKSHARE
      out(:,:,:) = f(:,:,:)*filter(:,:,:)
      !$OMP END WORKSHARE
    else
      !$OMP WORKSHARE
      f(:,:,:) = f(:,:,:)*filter(:,:,:)
      !$OMP END WORKSHARE
    endif

    !$OMP SINGLE
    call timer_stop(filterhandle)
    !$OMP END SINGLE

  end subroutine

  !gives spectral derivative in the x direction: out = 2*pi*i*kx*in
  !This requires us to swap the real and imaginary parts of the numbers around - essentially swap
  !even and odd array elements around. The catch is that these arrays are decomposed between
  !processes and some complex number pairs may be split between processes, so we sometimes need to
  !send messages between processes to swap these pairs round
  subroutine diffx(in,out)
    real(kind=DEFAULT_PRECISION), intent(in) :: in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(out) :: out(:,:,:)
    integer :: left_sendrequest, left_recvrequest, right_sendrequest, right_recvrequest
    integer, save :: istart, iend
    integer :: i, j, k
    integer :: statuses(MPI_STATUS_SIZE,4)
    integer :: requests(4)=MPI_REQUEST_NULL




    !$OMP SINGLE
    call timer_start(diffxyhandle)
    istart=1
    iend=nx
    !send/recv start/end values of arrays if needed (non-blocking)
    if (x_start_swap) then
      left_sendbuff(:,:) = in(:,:,1)
      call MPI_isend(left_sendbuff,&
                     nz*ny,&
                     PRECISION_TYPE,&
                     left,&
                     0,&
                     comm,&
                     requests(1),&
                     ierr)
      call MPI_irecv(left_recvbuff,&
                     nz*ny,&
                     PRECISION_TYPE,&
                     left,&
                     1,&
                     comm,&
                     requests(2),&
                     ierr)
      istart=2
    endif

    if (x_end_swap) then
      right_sendbuff(:,:) = in(:,:,nx)
      call MPI_isend(right_sendbuff,&
                     nz*ny,&
                     PRECISION_TYPE,&
                     right,&
                     1,&
                     comm,&
                     requests(3),&
                     ierr)
      call MPI_irecv(right_recvbuff,&
                     nz*ny,&
                     PRECISION_TYPE,&
                     right,&
                     0,&
                     comm,&
                     requests(4),&
                     ierr)
      iend=nx-1
    endif
    !$OMP END SINGLE

    !swap internal values (i.e. multiply by i=sqrt(-1))
    !$OMP DO
    do i=istart,iend,2
      out(:,:,i)=-1.*in(:,:,i+1)
      out(:,:,i+1) = in(:,:,i)
    enddo
    !$OMP END DO

    !wait for comms to complete (if necessary)

    !$OMP SINGLE
      call MPI_Waitall(4,&
                       requests,&
                       statuses,&
                       ierr)
    !$OMP END SINGLE

    !set end values (if necessary)
    if (x_start_swap) then
      !$OMP WORKSHARE
      out(:,:,1) = left_recvbuff(:,:)
      !$OMP END WORKSHARE
    endif
    if (x_end_swap) then
      !$OMP WORKSHARE
      out(:,:,nx) = -1.*right_recvbuff(:,:)
      !$OMP END WORKSHARE
    endif

    !multiply by k
    !$OMP WORKSHARE
    out(:,:,:) = out(:,:,:) * kx(:,:,:)
    !$OMP END WORKSHARE

    !$OMP SINGLE
    call timer_stop(diffxyhandle)
    !$OMP END SINGLE

  end subroutine

  !gives spectral derivative in the y direction: out = 2*pi*i*ky*in
  subroutine diffy(in,out)
    real(kind=DEFAULT_PRECISION), intent(in) :: in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(out) :: out(:,:,:)
    integer, save :: jstart, jend
    integer :: i, j, k
    integer :: statuses(MPI_STATUS_SIZE,4)
    integer :: requests(4)=MPI_REQUEST_NULL



    !$OMP SINGLE
    call timer_start(diffxyhandle)

    jstart=1
    jend=ny
    !send/recv start/end values of arrays if needed (non-blocking)
    if (y_start_swap) then
      down_sendbuff(:,:) = in(:,1,:)
      call MPI_isend(down_sendbuff,&
                     nz*nx,&
                     PRECISION_TYPE,&
                     down,&
                     0,&
                     comm,&
                     requests(1),&
                     ierr)
      call MPI_irecv(down_recvbuff,&
                     nz*nx,&
                     PRECISION_TYPE,&
                     down,&
                     1,&
                     comm,&
                     requests(2),&
                     ierr)
      jstart=2
    endif

    if (y_end_swap) then
      up_sendbuff(:,:) = in(:,ny,:)
      call MPI_isend(up_sendbuff,&
                     nz*nx,&
                     PRECISION_TYPE,&
                     up,&
                     1,&
                     comm,&
                     requests(3),&
                     ierr)
      call MPI_irecv(up_recvbuff,&
                     nz*nx,&
                     PRECISION_TYPE,&
                     up,&
                     0,&
                     comm,&
                     requests(4),&
                     ierr)
      jend=ny-1
    endif
    !$OMP END SINGLE

    !swap internal values (i.e. multiply by i=sqrt(-1))
    !$OMP DO
    do j=jstart,jend,2
      out(:,j,:)=-1.*in(:,j+1,:)
      out(:,j+1,:) = in(:,j,:)
    enddo
    !$OMP END DO

    !wait for comms to complete
    !$OMP SINGLE
      call MPI_Waitall(4,&
                       requests,&
                       statuses,&
                       ierr)
    !$OMP END SINGLE

    !set end values (if necessary)
    if (y_start_swap) then
      !$OMP WORKSHARE
      out(:,1,:) = down_recvbuff(:,:)
      !$OMP END WORKSHARE
    endif
    if (y_end_swap) then
      !$OMP WORKSHARE
      out(:,ny,:) = -1.*up_recvbuff(:,:)
      !$OMP END WORKSHARE
    endif

    !multiply by k
    !$OMP WORKSHARE
    out(:,:,:) = out(:,:,:) * ky(:,:,:)
    !$OMP END WORKSHARE

    !$OMP SINGLE
    call timer_stop(diffxyhandle)
    !$OMP END SINGLE

  end subroutine



  ! 4th order accurate derivative in the z direction (returns df/dz)
  ! Solves the tridigonal problem:
  ! 1/6 dfdz(i-1) + 2/3 dfdz(i) + 1/6 dfdz(i+1) = (f(i+1)-f(i-1))/2 for i=2,n-1
  !
  ! If the optional arguments "bot" and "top" are left out then it assumes zero gradient
  ! the boundaries:
  ! dfdz(1) = 0
  ! dfdz(n) = 0
  !
  ! If "bot" and "top" are included, then it assumes that f must be zero at the boundaries,
  ! and the boundary values are:
  ! 2/3 dfdz(1) + 1/3 dfdz(2) = f(2)/dz - 1/6 bot(1)*dz
  ! 1/3 dfdz(n-1) + 2/3 dfdz(n) = -f(n-1)/dz + 1/6 top(n) * dz
  !
  subroutine diffz(f,dfdz,bot,top)
    real(kind=DEFAULT_PRECISION), intent(in) :: f(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(out) :: dfdz(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(in), optional :: bot(:,:), top(:,:)

    integer :: nz,nx,ny,i,j,k
    real(kind=DEFAULT_PRECISION), allocatable :: d(:), a(:), b(:), c(:)
    real(kind=DEFAULT_PRECISION), parameter :: a0=0., b0=1., c0=0.
    real(kind=DEFAULT_PRECISION), parameter :: a1=1./6., b1=2./3., c1=1./6.
    real(kind=DEFAULT_PRECISION), parameter :: an=0., bn=1., cn=0.
    real(kind=DEFAULT_PRECISION), parameter :: al = 1./3., cl=1./3.

    !$OMP SINGLE
    call timer_start(diffzhandle)
    !$OMP END SINGLE

    nz=size(f,1)
    ny=size(f,2)
    nx=size(f,3)

    allocate(a(nz), b(nz), c(nz), d(nz))

    !setup a, b and c arrays (assuming dfdz=0 on boundaries for now)
    a(1)=a0
    b(1)=b0
    c(1)=c0
    a(2:nz-1)=a1
    b(2:nz-1)=b1
    c(2:nz-1)=c1
    a(nz)=an
    b(nz)=bn
    c(nz)=cn

    if (present(bot) .and. present(top)) then
      c(1)=cl
      a(nz)=al
    endif

    !loop over all columns and calculate df/dz
    !$OMP DO
    do i=1,nx
      do j=1,ny
        !set up d array for column
        if (present(bot) .and. present(top)) then
          d(1) = f(2,j,i)/dz - a1* bot(j,i)*dz
          d(nz)=-f(nz-1,j,i)/dz + a1 * top(j,i) * dz
        else
          d(1)=0.d0
          d(nz)=0.d0
        endif
        do k=2,nz-1
          d(k)=(f(k+1,j,i)-f(k-1,j,i))/2./dz
        enddo
        !solve for df/dz
        call tridiagonal(a,b,c,d,dfdz(:,j,i))
      enddo
    enddo
    !$OMP END DO
    deallocate(a,b,c,d)

    !$OMP SINGLE
    call timer_stop(diffzhandle)
    !$OMP END SINGLE

  end subroutine

  !inverts B = laplacian(A) using the tridiagonal problem:
  ! (1/dz^2 +-K^2/12) A(i-1) + (-2/dz^2 -5K^2/6) A(i) + (1/dz^2 - K^2/12) A(i+1)
  ! = 1/12 B(i-1) + 5/6 B(i) + 1/12 B(i+1) for i=2,n-1
  ! where K^2 = k^2 + l^2 is the square of the horizontal wavenumbers k and l
  !
  ! Boundary conditions are:
  !
  ! - f_zero_on_boundary
  !   We assume A is zero on the boundaries. This results in
  !   A(1) = 0.
  !   A(n) = 0.
  !
  ! - df_zero_on_boundary
  !   We assume dA/dz is zero on the boundary.
  !   (-1/dz^2 -K^2/3) A(1)   + (1/dz^2 - K^2/6) A(2)  = 1/3 B(1) + 1/6 B(2)
  !   (1/dz^2 - K^2/6) A(n-1) + (-1/dz^2 - K^2/3) A(N) = 1/6 B(n-1) + 1/3 B(n)
  !
  ! This subroutine solves for the laplacian in place, so we pass in f as B, and
  ! A is returned in f

  subroutine laplinv(f,f_zero_on_boundary,df_zero_on_boundary)
    real(kind=DEFAULT_PRECISION), intent(inout) :: f(:,:,:)
    logical, intent(in), optional :: f_zero_on_boundary, df_zero_on_boundary

    logical :: df0, f0

    real(kind=DEFAULT_PRECISION), allocatable :: a(:), b(:), c(:), d(:)
    real(kind=DEFAULT_PRECISION), parameter ::f16=1./6., f23=2./3., f13=1./3.,&
                                              f112=1./12., f56=5./6.
    real(kind=DEFAULT_PRECISION) :: i2dz2, idz2
    integer :: nx, ny, nz
    integer :: i, j, k

    !$OMP SINGLE
    call timer_start(laplinvhandle)
    !$OMP END SINGLE


    df0=.false.
    f0=.false.

    !check that we have boundary conditions specified

    if (.not. present(f_zero_on_boundary) .and. .not. present(df_zero_on_boundary)) then
      write(*,*) "Error no boundary conditions specified for laplinv"
      call MPI_Finalize(ierr)
      stop
    endif

    ! get the boundary conditions
    if (present(f_zero_on_boundary)) f0 =f_zero_on_boundary
    if (present(df_zero_on_boundary)) df0 = df_zero_on_boundary

    ! check that we don't have multiple BCs
    if (f0 .AND. df0) then
     write(*,*) "Error multiple boundary conditions specified for laplinv"
     call MPI_Finalize(ierr)
     stop
    endif

    ! check that we have BCs    
    if ((.NOT.(f0 .OR. df0)) ) then
     write(*,*) "No boundary conditions specified for laplinv"
     call MPI_Finalize(ierr)
     stop
    endif
    
    !set any constants

    nx=size(f,3)
    ny=size(f,2)
    nz=size(f,1)
    i2dz2 = 2./dz/dz
    idz2 = 1./dz/dz

    !allocate tridiagonal coefficient arrays
    allocate(a(nz), b(nz), c(nz), d(nz))

    !loop over all columns
    !$OMP DO
    do i=1,nx
      do j=1,ny
        if (k2(1,j,i) .eq. 0) cycle !We don't want to calculate the case where k^2=0
        !specify a, b, c and d on column
        if (f0) then !f is zero on boundary
          a(1)=0.
          b(1)=1.
          c(1)=0.
          d(1)=0.
          a(nz)=0.
          b(nz)=1.
          c(nz)=0.
          d(nz)=0.
        else !df/dz is zero on boundary
          a(1) = 0.
          b(1) = -idz2 - f13*k2(1,j,i)
          c(1) =  idz2 - f16*k2(1,j,i)
          d(1) = f13*f(1,j,i) + f16*f(2,j,i)
          a(nz) = idz2 - f16*k2(nz,j,i)
          b(nz) = -idz2 - f13*k2(nz,j,i)
          c(nz) = 0.
          d(nz) = f13*f(nz,j,i) + f16*f(nz-1,j,i)
        endif
        a(2:nz-1) = idz2 - f112*k2(2:nz-1,j,i)
        b(2:nz-1) = -i2dz2 - f56*k2(2:nz-1,j,i)
        c(2:nz-1) = idz2 - f112*k2(2:nz-1,j,i)
        d(2:nz-1) = f112*f(1:nz-2,j,i) + f56*f(2:nz-1,j,i) + f112*f(3:nz,j,i)

        !invert using tridiagonal method
        call tridiagonal(a,b,c,d,f(:,j,i))
      enddo
    enddo
    !$OMP ENDDO
    !$OMP SINGLE
    call timer_stop(laplinvhandle)
    !$OMP END SINGLE
  end subroutine



  !solves a tridiagonal problem A*f = d, where A is a tridiagonal matrix with
  !coefficients a, b and c (a(1) and c(n) are zero as they do not exist)

  ! Algorithm:
  ! We have a(i)*f(i-1) + b(i)*f(i) + c(i)*f(i+1) = d(i), where we want to solve for f.
  !
  ! To do this we update c and d to:
  ! c'(i) = c(i) / (b(i)-a(i)*c'(i-1))
  ! d'(i) = (d(i) - a(i)*d'(i-1)) / (b(i)-a(i)*c'(i-1))
  !
  ! Then finally get f from:
  ! f(i) = d'(i) - c'(i)*f(i+1)
  !
  ! (see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm for more details)
  !
  subroutine tridiagonal(a,b,c,d,f)
    real(kind=DEFAULT_PRECISION), intent(in) :: a(:), b(:), c(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: d(:)
    real(kind=DEFAULT_PRECISION), intent(out):: f(:)

    integer :: i, n
    real(kind=DEFAULT_PRECISION), allocatable :: cp(:)
    real(kind=DEFAULT_PRECISION) :: denominator

    n=size(f)

    allocate(cp(n))

    ! we use cp(i) to represent c'(i-1)
    ! we use f(i) to represent d'(i)

    denominator=b(1)
    f(1)=d(1)/denominator

    do i=2,n
      cp(i) = c(i-1) / denominator !cp(i) is c'(i-1)
      denominator = b(i) - a(i) * cp(i)
      if (denominator .eq. 0) then
        print *, "Error in tridiagonal solver: divide by zero"
        error stop
      endif
      f(i) = (d(i) - a(i)*f(i-1))/(denominator)
    enddo


    do i=n-1,1,-1
      f(i) = f(i)-cp(i+1)*f(i+1)
    enddo

    deallocate(cp)

  end subroutine



end module
