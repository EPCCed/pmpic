!Halo swapping routines for MPIC functionality
!These halp swapping routines ar required in the par2grid and grid2par operations
!only certain halo locations need to be transferred, and occasionally summed rather than swapped
module MPIC_Haloswap_mod
  use state_mod, only: model_state_type
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use omp_lib
  use mpi
  use timer_mod


  implicit none

 !send/recv buffers
  real (kind=DEFAULT_PRECISION), allocatable, dimension (:,:) :: left_buf, right_buf, up_buf, down_buf
  real (kind=DEFAULT_PRECISION), allocatable, dimension(:) :: upper_corner_buf, lower_corner_buf

  integer :: requests(3)
  integer :: statuses(MPI_STATUS_SIZE,3)

  integer :: hx, hy
  integer :: nx,ny,nz
  integer :: xi,xf, yi, yf

  integer :: left, right, up, down, upper_corner, lower_corner

  integer, parameter :: xtag = 1, ytag=2, cortag=3

  integer :: comm, ierr

  integer :: g2p_handle, p2g_handle, mixing_handle

  logical :: mixing





contains

  subroutine MPIC_Haloswap_init(state)
    type(model_state_type), intent(inout) :: state

    real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:) :: data


    !get halo sizes and grid sizes
    hx=state%local_grid%halo_size(3)
    hy=state%local_grid%halo_size(2)

    nx=state%local_grid%size(3) !NOTE: only internal cell count (so not including halo cells)
    ny=state%local_grid%size(2)
    nz=state%local_grid%size(1)

    !start and end indices for the local grid data

    xi=hx+1
    xf=nx+hx

    yi=hy+1
    yf=ny+hy

    !allocate buffers

    allocate(left_buf(nz,ny), right_buf(nz,ny), up_buf(nz,nx), down_buf(nz,nx))
    allocate(upper_corner_buf(nz), lower_corner_buf(nz))

    ! define neighbouring processors
    down=state%local_grid%neighbours(2,1)
    up=state%local_grid%neighbours(2,3)
    left=state%local_grid%neighbours(3,1)
    right=state%local_grid%neighbours(3,3)

    lower_corner=state%local_grid%corner_neighbours(1,1)
    upper_corner=state%local_grid%corner_neighbours(4,1)

    comm=state%parallel%monc_communicator

    call register_routine_for_timing("par2grid_haloswap",p2g_handle,state)
    call register_routine_for_timing("grid2par_haloswap",g2p_handle,state)
    call register_routine_for_timing("mixing_haloswap",mixing_handle,state)

    mixing = .false.

    if (state%parallel%my_rank .eq. 0) print *, "MPIC_Haloswap initialised"
    ! if (state%parallel%my_rank .eq. 0) print *, xi, xf, yi, yf
    !
    ! return
    !
    !
    ! allocate(data(nz, ny + 2*hy, nx + 2*hx))
    !
    ! data(:,:,:) = 0.
    ! data(:,yi:yf,xi:xf) = state%parallel%my_rank+1
    !
    ! if (state%parallel%my_rank .eq. 0) then
    !   print *, "up=",up
    !   print *, "down=",down
    !   print *, "left=",left
    !   print *, "right=",right
    !   print *, "upcorn=", upper_corner
    !   print *, "lowcorn=", lower_corner
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf+1,xi-1), data(1, yf+1,xi),data(1, yf+1,xf),data(1, yf+1,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf,xi-1), data(1, yf,xi),data(1, yf,xf),data(1, yf,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi,xi-1), data(1, yi,xi),data(1, yi,xf),data(1, yi,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi-1,xi-1), data(1, yi-1,xi),data(1, yi-1,xf),data(1, yi-1,xf+1)
    ! endif
    !
    ! call grid2par_haloswap(state,data)
    !
    ! if (state%parallel%my_rank .eq. 0) then
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf+1,xi-1), data(1, yf+1,xi),data(1, yf+1,xf),data(1, yf+1,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf,xi-1), data(1, yf,xi),data(1, yf,xf),data(1, yf,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi,xi-1), data(1, yi,xi),data(1, yi,xf),data(1, yi,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi-1,xi-1), data(1, yi-1,xi),data(1, yi-1,xf),data(1, yi-1,xf+1)
    ! endif
    !
    ! data = state%parallel%my_rank+1
    !
    ! if (state%parallel%my_rank .eq. 0) then
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf+1,xi-1), data(1, yf+1,xi),data(1, yf+1,xf),data(1, yf+1,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yf,xi-1), data(1, yf,xi),data(1, yf,xf),data(1, yf,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi,xi-1), data(1, yi,xi),data(1, yi,xf),data(1, yi,xf+1)
    !   write(*,"(f2.0,1x,f2.0,1x,f2.0,1x,f2.0)") data(1,yi-1,xi-1), data(1, yi-1,xi),data(1, yi-1,xf),data(1, yi-1,xf+1)
    ! endif
    !
    ! call par2grid_haloswap(state,data)
    !
    !
    ! if (state%parallel%my_rank .eq. 0) then
    !   write(*,"(f3.0,1x,f3.0,1x,f3.0,1x,f3.0)") data(1,yf+1,xi-1), data(1, yf+1,xi),data(1, yf+1,xf),data(1, yf+1,xf+1)
    !   write(*,"(f3.0,1x,f3.0,1x,f3.0,1x,f3.0)") data(1,yf,xi-1), data(1, yf,xi),data(1, yf,xf),data(1, yf,xf+1)
    !   write(*,"(f3.0,1x,f3.0,1x,f3.0,1x,f3.0)") data(1,yi,xi-1), data(1, yi,xi),data(1, yi,xf),data(1, yi,xf+1)
    !   write(*,"(f3.0,1x,f3.0,1x,f3.0,1x,f3.0)") data(1,yi-1,xi-1), data(1, yi-1,xi),data(1, yi-1,xf),data(1, yi-1,xf+1)
    ! endif
    !
    !
    !
    ! call MPI_Finalize(ierr)
    ! stop




  end subroutine





! send grid data to halos of processes to the left and downwards
  subroutine grid2par_haloswap(state, array)
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: array
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, tag, src

    if (.not. mixing) call timer_start(g2p_handle)

    !copy data to buffers
    left_buf(:,:) = array(:,yi:yf,xi)
    down_buf(:,:) = array(:,yi,xi:xf)
    lower_corner_buf(:) = array(:,yi,xi)

    !send data (non-blocking)

    Call MPI_Isend(left_buf,nz*ny,PRECISION_TYPE,left,xtag,comm,requests(xtag),ierr)
    Call MPI_Isend(down_buf,nz*nx,PRECISION_TYPE,down,ytag,comm,requests(ytag),ierr)
    Call MPI_Isend(lower_corner_buf,nz,PRECISION_TYPE,lower_corner,cortag,comm,requests(cortag),ierr)



    !now receive the data from neighbouring processes

    do i=1,3

      !determine metadata of incoming message (who and where it's from)
      call MPI_probe(MPI_ANY_SOURCE,&
                     MPI_ANY_TAG,&
                     comm,&
                     status,&
                     ierr)

      src=status(MPI_SOURCE)
      tag=status(MPI_TAG)

      if ( tag .eq. xtag) then
        if (src .ne. right) then
          print *, "unexpected x source"
          error stop
        endif

        call MPI_Recv(right_buf,nz*ny,PRECISION_TYPE,right,xtag,comm,status,ierr)

        !unpack data

        array(:,yi:yf,xf+1) = right_buf(:,:)

      else if (tag .eq. ytag) then
        if (src .ne. up) then
          print *, "unexpected y source"
          error stop
        endif

        call MPI_Recv(up_buf,nz*nx,PRECISION_TYPE,up,ytag,comm,status,ierr)

        !unpack data

        array(:,yf+1,xi:xf) = up_buf(:,:)

      else if (tag .eq. cortag) then
        if (src .ne. upper_corner) then
          print *, "unexpected corner source"
          error stop
        endif

        call MPI_Recv(upper_corner_buf,nz,PRECISION_TYPE,upper_corner,cortag,comm,status,ierr)

        !unpack data

        array(:,yf+1,xf+1) = upper_corner_buf(:)

      else
        print *, "Invalid message request"
      endif

    enddo

    !make sure everything is swapped

    call MPI_Waitall(3,requests,statuses,ierr)

    if (.not. mixing) call timer_stop(g2p_handle)

    call MPI_Barrier(comm,ierr)

    

    !if (state%parallel%my_rank .eq. 0) print *, "grid2par haloswapping successful"

  end subroutine grid2par_haloswap


! send halo data to processes to the right and upwards, and appaned this to their edge cells
  subroutine par2grid_haloswap(state, array)
    implicit none
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: array
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, tag, src, xcount, ycount, ccount

    if (.not. mixing) call timer_start(p2g_handle)

    !copy data to buffers
    right_buf(:,:) = array(:,yi:yf,xf+1)
    up_buf(:,:) = array(:,yf+1,xi:xf)
    upper_corner_buf(:) = array(:,yf+1,xf+1)

    !send data (non-blocking)

    Call MPI_Isend(right_buf,nz*ny,PRECISION_TYPE,right,xtag,comm,requests(xtag),ierr)
    Call MPI_Isend(up_buf,nz*nx,PRECISION_TYPE,up,ytag,comm,requests(ytag),ierr)
    Call MPI_Isend(upper_corner_buf,nz,PRECISION_TYPE,upper_corner,cortag,comm,requests(cortag),ierr)



    !now receive the data from neighbouring processes

    xcount=0
    ycount=0
    ccount=0

    do i=1,3

      !determine metadata of incoming message (who and where it's from)
      call MPI_probe(MPI_ANY_SOURCE,&
                     MPI_ANY_TAG,&
                     comm,&
                     status,&
                     ierr)

      src=status(MPI_SOURCE)
      tag=status(MPI_TAG)

      if ( tag .eq. xtag) then
        if (src .ne. left) then
          print *, "unexpected x source"
          error stop
        endif
        if (xcount .gt. 1) error stop "xcount too high"

        call MPI_Recv(left_buf,nz*ny,PRECISION_TYPE,left,xtag,comm,status,ierr)

        !unpack data

        array(:,yi:yf,xi) = array(:,yi:yf,xi) + left_buf(:,:)
        xcount=xcount+1

      else if (tag .eq. ytag) then
        if (src .ne. down) then
          print *, "unexpected y source"
          error stop
        endif
        if (ycount .gt. 1) error stop "ycount too high"

        call MPI_Recv(down_buf,nz*nx,PRECISION_TYPE,down,ytag,comm,status,ierr)
        ycount=ycount+1
        !unpack data

        array(:,yi,xi:xf) = array(:,yi,xi:xf) + down_buf(:,:)

      else if (tag .eq. cortag) then
        if (src .ne. lower_corner) then
          print *, "unexpected corner source"
          error stop
        endif
        if (ccount .gt. 1) error stop "ccount too high"

        call MPI_Recv(lower_corner_buf,nz,PRECISION_TYPE,lower_corner,cortag,comm,status,ierr)

        !unpack data

        array(:,yi,xi) = array(:,yi,xi) + lower_corner_buf(:)

        ccount=ccount+1

      else
        print *, "Invalid message request"
        error stop
      endif

    enddo

    !make sure everything is swapped

    call MPI_Waitall(3,requests,statuses,ierr)

    if (.not. mixing) call timer_stop(p2g_handle)

    !We seem to need this barrier to prevent processes from running away from each other
    call MPI_Barrier(comm,ierr)

    

    !if (state%parallel%my_rank .eq. 0) print *, "par2grid haloswapping successful"

  end subroutine par2grid_haloswap


  ! a haloswap combining grid2par and par2grid haloswaps
  subroutine mixing_haloswap(state,grid)
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:) :: grid

    call timer_start(mixing_handle)
    mixing = .true.

    call par2grid_haloswap(state,grid)
    call grid2par_haloswap(state,grid)

    mixing = .false.
    call timer_stop(mixing_handle)

  end subroutine










  end module
