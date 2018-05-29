module parcel_haloswap_mod

  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use MPI

  implicit none

  logical :: initialised = .false.

  integer :: ierr

  !send buffers
  real (kind=DEFAULT_PRECISION), dimension(:,:), allocatable, target :: &
      N_buff, S_buff, E_buff, W_buff, NE_buff, NW_buff, SE_buff, SW_buff

  !index for which parcel to send where (value is either 0 (don't send) or a direction)
  integer, allocatable, dimension(:) :: index


  !directions to send in (N = +y, E = +x)
  integer, parameter :: N = 1
  integer, parameter :: NE = 2
  integer, parameter :: E = 3
  integer, parameter :: SE = 4
  integer, parameter :: S = 5
  integer, parameter :: SW = 6
  integer, parameter :: W = 7
  integer, parameter :: NW = 8

  !ranks to send to in each direction
  integer :: N_rank, NE_rank, E_rank, SE_rank, S_rank, SW_rank, W_rank, NW_rank



contains

  subroutine initialise_parcel_haloswapping(state)
    type(model_state_type), intent(inout) :: state

    if (initialised) then
      error stop "parcel haloswaping is initially initialised, cannot re-initialise"
    endif

    !determine neighbours

    S_rank=state%local_grid%neighbours(2,1) ! in -y direction
    N_rank=state%local_grid%neighbours(2,3) !in +y firection
    W_rank=state%local_grid%neighbours(3,1) ! in -x direction
    E_rank=state%local_grid%neighbours(3,3) ! in +x direction

    SW_rank=state%local_grid%corner_neighbours(1,1)
    SE_rank=state%local_grid%corner_neighbours(2,1)
    NW_rank=state%local_grid%corner_neighbours(3,1)
    NE_rank=state%local_grid%corner_neighbours(4,1)


    !allocate index array
    allocate(index(state%parcels%maxparcels_local))

    initialised = .true.

  end subroutine




  subroutine parcel_haloswap(state)
    implicit none
    type(model_state_type), intent(inout) :: state

    integer :: dir, nreceived, src
    integer :: nparcels_initial, nparcels_final, lastparcel
    integer :: i
    integer :: recv_status(MPI_STATUS_SIZE), send_statuses(MPI_STATUS_SIZE,8)
    !number of parcels to send/recv per direction and the request handle for sent messages
    integer :: nsend(8), nrecv(8), requests(8)
    !receive buffer
    real (kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: recv_buffer

    if (.not. initialised) error stop "parcel_haloswap not initialised"

    nparcels_initial = state%parcels%numparcels_local
    nparcels_final = nparcels_initial
    lastparcel=state%parcels%numparcels_local

    index(1:nparcels_initial) = 0
    index(nparcels_initial+1:state%parcels%maxparcels_local) = -1


    do dir=1,8
      call count_parcels_to_send(state,dir,nsend(dir))
    enddo

    print *, "to send=", nsend

    nparcels_final = nparcels_final - sum(nsend)

    !allocate the send buffers
    allocate(N_buff(nsend(N),state%parcels%n_properties))
    allocate(NE_buff(nsend(NE),state%parcels%n_properties))
    allocate(E_buff(nsend(E),state%parcels%n_properties))
    allocate(SE_buff(nsend(SE),state%parcels%n_properties))
    allocate(S_buff(nsend(S),state%parcels%n_properties))
    allocate(SW_buff(nsend(SW),state%parcels%n_properties))
    allocate(W_buff(nsend(W),state%parcels%n_properties))
    allocate(NW_buff(nsend(NW),state%parcels%n_properties))

    do dir=1,8
      call populate_send_buffer(state,dir)
    enddo

    do dir=1,8
      call send_buffers(state,dir,requests(dir),nsend(dir))
    enddo

    do i=1,8

      !see if a message has come in, who it came from and what its size is
      call check_for_message(state,dir,nreceived,src)

      allocate(recv_buffer(nreceived,state%parcels%n_properties))

      nparcels_final = nparcels_final + nreceived

      print *, "reading message from dir=", dir
      call MPI_Recv(buf=recv_buffer,&
                    count=nreceived*state%parcels%n_properties,&
                    datatype=PRECISION_TYPE,&
                    source=src,&
                    tag=dir,&
                    comm=state%parallel%monc_communicator,&
                    status=recv_status,&
                    ierror=ierr)

       !unpacks parcels into state%parcels, backfilling where possible
       call unpack_parcels(state,recv_buffer,nreceived, lastparcel,index)

       deallocate(recv_buffer)

     enddo

     ! if we end out with fewer parcels than we started with then we need to backfill
     if (nparcels_final .lt. nparcels_initial) then
       call backfill(state,index,lastparcel)
     endif

     !check that all our messages have been received
     call MPI_Waitall(count=8,&
                      array_of_requests=requests,&
                      array_of_statuses=send_statuses,&
                      ierror=ierr)

     !now we can safely deallocate the send buffers
     deallocate(N_buff)
     deallocate(NE_buff)
     deallocate(E_buff)
     deallocate(SE_buff)
     deallocate(S_buff)
     deallocate(SW_buff)
     deallocate(W_buff)
     deallocate(NW_buff)

     !update parcel number
     state%parcels%numparcels_local = nparcels_final

   end subroutine



   subroutine count_parcels_to_send(state,dir,nsend)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir
     integer, intent(out) :: nsend

     real (kind=DEFAULT_PRECISION) :: xmin, xmax, ymin, ymax
     real (kind=DEFAULT_PRECISION) :: x, y
     integer :: hy, hx, ny, nx
     integer :: i

     hy=state%local_grid%halo_size(2)
     hx=state%local_grid%halo_size(3)

     ny=state%local_grid%size(2) + 2*hy
     nx=state%local_grid%size(3) + 2*hx

     nsend=0

     ymax=y_coords(ny-hy+1)
     ymin=y_coords(hy)

     xmax=x_coords(nx-hx+1)
     xmin=x_coords(hx)

     if (dir .eq. N) then
       !looking for parcels y >= ymax within the x limits
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (y .ge. ymax) then
           if (x .ge. xmin .and. x .lt. xmax) then
             nsend=nsend+1
             if (index(i) .ne. 0) then
               print *, i,"index .ne. 0. Whoops!", dir,index(i)
             endif
             index(i) = dir
           endif
         endif
       enddo
     else if (dir .eq. S) then
       !looking for parcels y < ymin within the x limits
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (y .lt. ymin) then
           if (x .ge. xmin .and. x .lt. xmax) then
             nsend=nsend+1
             if (index(i) .ne. 0) then
               print *, i,"index .ne. 0. Whoops!", dir,index(i)
             endif
             index(i) = dir
           endif
         endif
       enddo
     else if (dir .eq. E) then
       !looking for parcels x >= xmax within the y limits
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .ge. xmax) then
           if (y .ge. ymin .and. y .lt. ymax) then
             nsend=nsend+1
             if (index(i) .ne. 0) then
               print *, i,"index .ne. 0. Whoops!", dir,index(i)
             endif
             index(i) = dir
           endif
         endif
       enddo
     else if (dir .eq. W) then
       !looking for parcels x < xmin within the y limits
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .lt. xmin) then
           if (y .ge. ymin .and. y .lt. ymax) then
             nsend=nsend+1
             if (index(i) .ne. 0) then
               print *, i,"index .ne. 0. Whoops!", dir,index(i)
             endif
             index(i) = dir
           endif
         endif
       enddo
     else if (dir .eq. NE) then
       !looking for parcels x >=xmax and y >= ymax
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .ge. xmax .and. y .ge. ymax) then
           nsend=nsend+1
           if (index(i) .ne. 0) then
             print *, i,"index .ne. 0. Whoops!", dir,index(i)
           endif
           index(i) = dir
         endif
       enddo
     else if (dir .eq. NW) then
       !looking for parcels x <xmin and y >= ymax
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .lt. xmin .and. y .ge. ymax) then
           nsend=nsend+1
           if (index(i) .ne. 0) then
             print *, i,"index .ne. 0. Whoops!", dir,index(i)
           endif
           index(i) = dir
         endif
       enddo
     else if (dir .eq. SE) then
       !looking for parcels x >=xmax and y < ymin
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .ge. xmax .and. y .lt. ymin) then
           nsend=nsend+1
           if (index(i) .ne. 0) then
             print *, i,"index .ne. 0. Whoops!", dir,index(i)
           endif
           index(i) = dir
         endif
       enddo
     else if (dir .eq. SW) then
       !looking for parcels x < xmin and y < ymin
       do i=1,state%parcels%numparcels_local
         x=state%parcels%x(i)
         y=state%parcels%y(i)
         if (x .lt. xmin .and. y .lt. ymin) then
           nsend=nsend+1
           if (index(i) .ne. 0) then
             print *, i,"index .ne. 0. Whoops!", dir,index(i)
           endif
           index(i) = dir
         endif
       enddo
     else
       print *, "wrong direction value", dir
       error stop
     endif

   end subroutine

   subroutine populate_send_buffer(state,dir)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir

     print *, "populate_send_buffer - nothing to see here yet"

   end subroutine

   subroutine send_buffers(state,dir,request,count)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir
     integer, intent(out) :: request
     integer, intent(in) :: count
     integer:: dest

     real(kind=DEFAULT_PRECISION), dimension(:,:), pointer :: buff

     if (dir .eq. N) then
       buff => N_buff
       dest = N_rank
     else if (dir .eq. NE) then
       buff => NE_buff
       dest=NE_rank
     else if (dir .eq. E) then
       buff => E_buff
       dest=E_rank
     else if (dir .eq. SE) then
       buff => SE_buff
       dest=SE_rank
     else if (dir .eq. S) then
       buff => S_buff
       dest=S_rank
     else if (dir .eq. SW) then
       buff => SW_buff
       dest=SW_rank
     else if (dir .eq. W) then
       buff => W_buff
       dest=W_rank
     else if (dir .eq. NW) then
       buff => NW_buff
       dest=NW_rank
     endif

     call MPI_ISend(buf=buff(1,1),&
                    count=count*state%parcels%n_properties,&
                    datatype=PRECISION_TYPE,&
                    dest=dest,&
                    tag=dir,&
                    comm=state%parallel%monc_communicator,&
                    request=request,&
                    ierror=ierr)

      print *, "sent in direction", dir

    end subroutine


    subroutine check_for_message(state,dir,nrecv,src)
      type(model_state_type), intent(in) :: state
      integer, intent(out) :: dir, nrecv, src
      integer :: status(MPI_STATUS_SIZE)

      call MPI_probe(source=MPI_ANY_SOURCE,&
                     tag=MPI_ANY_TAG,&
                     comm=state%parallel%monc_communicator,&
                     status=status,&
                     ierror=ierr)

      src=status(MPI_SOURCE)
      dir=status(MPI_TAG)

      call MPI_get_count(status=status,datatype=PRECISION_TYPE,count=nrecv,ierror=ierr)

      nrecv=nrecv/state%parcels%n_properties

      print *, "message: src, tag, size=", src, dir, nrecv



    end subroutine




    subroutine unpack_parcels(state,buff,nrecv,lastparcel,index)
      type(model_state_type), intent(inout) :: state
      real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: buff
      integer, intent(in) :: nrecv
      integer, intent(inout) :: lastparcel
      integer, dimension(:), intent(inout) :: index

      print *, "unpack parcels - nothing to see here yet"
    end subroutine

    subroutine backfill(state,index,lastparcel)
      type(model_state_type), intent(inout) :: state
      integer, dimension(:), intent(inout) :: index
      integer, intent(inout) :: lastparcel

      print *, "Backfill - nothing to be seen here yet"

    end subroutine

  !have array of maxparcels size (logical) containing t/f if parcel is to be removed or not

  !loop through parcels and identify parcels to be moved (another index array of integers to signify where?)

  !allocate buffers

  !populate buffers (can use OMP SECTIONS)

  !call MPI_Isend for each buffer

  !do i=1,4
     !call MPI_Probe (tag=any)
     !allocate recv_buffer
     !call mpi_recv

     !unpack parcels and backfill

  !do any necessary final backfilling





end module
