module parcel_haloswap_mod

  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, MPI_PARCEL_INT, PARCEL_INTEGER
  use MPI
  use timer_mod
  use parcel_interpolation_mod, only : zmin,zmax

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

  !parameters defining the index in the send/recv buffer that each parcel property belongs to
  integer, parameter :: X_INDEX=1
  integer, parameter :: Y_INDEX=2
  integer, parameter :: Z_INDEX=3
  integer, parameter :: P_INDEX=4
  integer, parameter :: Q_INDEX=5
  integer, parameter :: R_INDEX=6
  integer, parameter :: DXDT_INDEX=7
  integer, parameter :: DYDT_INDEX=8
  integer, parameter :: DZDT_INDEX=9
  integer, parameter :: DPDT_INDEX=10
  integer, parameter :: DQDT_INDEX=11
  integer, parameter :: DRDT_INDEX=12
  integer, parameter :: H_INDEX=13
  integer, parameter :: B_INDEX=14
  integer, parameter :: VOL_INDEX=15
  integer, parameter :: STRETCH_INDEX=16
  integer, parameter :: TAG_INDEX=17
  integer, parameter :: Q_START_INDEX=18

  integer :: destinations(8)

  integer :: handle
  integer :: handle_count
  integer :: handle_buffer
  integer :: handle_send
  integer :: handle_recv
  integer :: handle_unpack
  integer :: handle_backfill
  integer :: handle_sanity



contains

  subroutine initialise_parcel_haloswapping(state)
    type(model_state_type), intent(inout) :: state

    if (initialised) then
      error stop "parcel haloswaping is initially initialised, cannot re-initialise"
    endif

    call register_routine_for_timing("Par_HSWP_entire",handle,state)
    call register_routine_for_timing("Par_HSWP_count",handle_count,state)
    call register_routine_for_timing("Par_HSWP_buffer",handle_buffer,state)
    call register_routine_for_timing("Par_HSWP_send",handle_send,state)
    call register_routine_for_timing("Par_HSWP_recv",handle_recv,state)
    call register_routine_for_timing("Par_HSWP_unpack",handle_unpack,state)
    call register_routine_for_timing("Par_HSWP_backfill",handle_backfill,state)
    call register_routine_for_timing("Par_HSWP_sanity_chk",handle_sanity,state)


    !determine neighbours

    destinations(S)=state%local_grid%neighbours(2,1) ! in -y direction
    destinations(N)=state%local_grid%neighbours(2,3) !in +y firection
    destinations(W)=state%local_grid%neighbours(3,1) ! in -x direction
    destinations(E)=state%local_grid%neighbours(3,3) ! in +x direction

    destinations(SW)=state%local_grid%corner_neighbours(1,1)
    destinations(SE)=state%local_grid%corner_neighbours(2,1)
    destinations(NW)=state%local_grid%corner_neighbours(3,1)
    destinations(NE)=state%local_grid%corner_neighbours(4,1)


    !allocate index array
    allocate(index(state%parcels%maxparcels_local))

    initialised = .true.

  end subroutine




  subroutine parcel_haloswap(state)
    implicit none
    type(model_state_type), intent(inout) :: state

    integer :: dir, nreceived, src
    integer(kind=PARCEL_INTEGER) :: nparcels_initial, nparcels_final
    integer(kind=PARCEL_INTEGER) :: global_initial, global_final
    integer :: i
    integer :: recv_status(MPI_STATUS_SIZE), send_statuses(MPI_STATUS_SIZE,8)
    !number of parcels to send/recv per direction and the request handle for sent messages
    integer :: nsend(8), nrecv(8), requests(8)
    !receive buffer
    real (kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: recv_buffer

    integer(kind=PARCEL_INTEGER) :: lastparcel  !the last parcel we replaced (initially set to 1)
    integer :: q
    if (.not. initialised) error stop "parcel_haloswap not initialised"

    call timer_start(handle)

    nparcels_initial = state%parcels%numparcels_local
    nparcels_final = nparcels_initial
    lastparcel=state%parcels%numparcels_local

    global_initial = state%parcels%numparcels_global

    index(1:nparcels_initial) = 0
    !!index(nparcels_initial+1:state%parcels%maxparcels_local) = -1

    lastparcel=1
    
    call parcels_reflect_vertical_boundaries(state)
    
    !$OMP PARALLEL default(shared)

    !$OMP MASTER
    call timer_start(handle_count)
    !$OMP END MASTER

    !count the parcels to send, and also create an index if what parcel to send where
    !$OMP DO schedule(dynamic,1)
    do dir=1,8
      call count_parcels_to_send(state,dir,nsend(dir))
    enddo
    !$OMP END DO

    !$OMP MASTER
    call timer_stop(handle_count)


    nparcels_final = nparcels_final - sum(nsend)

    !allocate the send buffers
    allocate(N_buff(nsend(N),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(NE_buff(nsend(NE),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(E_buff(nsend(E),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(SE_buff(nsend(SE),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(S_buff(nsend(S),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(SW_buff(nsend(SW),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(W_buff(nsend(W),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))
    allocate(NW_buff(nsend(NW),state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))


    call timer_start(handle_buffer)
    !$OMP END MASTER

    !$OMP BARRIER


    do dir=1,8
      if (nsend(dir) .gt. 0) call populate_send_buffer(state,dir,nsend(dir))
    enddo


    !$OMP MASTER
    call timer_stop(handle_buffer)
    call timer_start(handle_send)
    do dir=1,8
      call send_buffers(state,dir,requests(dir),nsend(dir))
    enddo
    call timer_stop(handle_send)
    !$OMP END MASTER


    do i=1,8
      !see if a message has come in, who it came from and what its size is
      !$OMP MASTER
      call timer_start(handle_recv)
      call check_for_message(state,dir,nreceived,src)

      allocate(recv_buffer(nreceived,state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum))

      nparcels_final = nparcels_final + nreceived

      !print *, "reading message from dir=", dir
      call MPI_Recv(recv_buffer,&
                    nreceived*(state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum),&
                    PRECISION_TYPE,&
                    src,&
                    dir,&
                    state%parallel%monc_communicator,&
                    recv_status,&
                    ierr)

      call timer_stop(handle_recv)

      !$OMP END MASTER

      !$OMP BARRIER

       !unpacks parcels into state%parcels, backfilling where possible
       if (nreceived .gt. 0) call unpack_parcels(state,recv_buffer,nreceived, lastparcel,nparcels_initial)



       !$OMP SINGLE
       deallocate(recv_buffer)
       !$OMP END SINGLE

     enddo




     ! if we end out with fewer parcels than we started with then we need to backfill
     if (nparcels_final .lt. nparcels_initial) then
       !$OMP SINGLE
       call timer_start(handle_backfill)
       !$OMP END SINGLE
       call backfill(state,index,nparcels_initial-nparcels_final)
       !$OMP SINGLE
       call timer_stop(handle_backfill)
       !$OMP END SINGLE
     endif

     !$OMP MASTER
     !check that all our messages have been received
     call timer_start(handle_sanity)
     call MPI_Waitall(8,&
                      requests,&
                      send_statuses,&
                      ierr)

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



     !a sanity check. We're moving parcels around so the total number across all processes shouldn't change
     call MPI_Allreduce(state%parcels%numparcels_local,&
                        state%parcels%numparcels_global,&
                        1,&
                        MPI_PARCEL_INT,&
                        MPI_SUM,&
                        state%parallel%monc_communicator,&
                        ierr)

      call timer_stop(handle_sanity)

      if (global_initial .ne. state%parcels%numparcels_global) then
        print *, "parcel number not conserved"
        call MPI_Finalize(ierr)
        error stop
      endif

      !$OMP END MASTER
      !$OMP END PARALLEL

      call timer_stop(handle)

   end subroutine






   subroutine count_parcels_to_send(state,dir,nsend)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir
     integer, intent(out) :: nsend

     real (kind=DEFAULT_PRECISION) :: xmin, xmax, ymin, ymax
     real (kind=DEFAULT_PRECISION) :: x, y
     integer :: hy, hx, ny, nx
     integer(kind=PARCEL_INTEGER) :: i

     hy=state%local_grid%halo_size(2)
     hx=state%local_grid%halo_size(3)

     ny=state%local_grid%size(2) + 2*hy
     nx=state%local_grid%size(3) + 2*hx

     nsend=0

     ymax=y_coords(ny-hy+1)
     ymin=y_coords(hy+1)

     xmax=x_coords(nx-hx+1)
     xmin=x_coords(hx+1)

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


   subroutine parcels_reflect_vertical_boundaries(state)
     type(model_state_type), intent(inout) :: state


     real (kind=DEFAULT_PRECISION) :: z
     integer(kind=PARCEL_INTEGER) :: i

     do i=1,state%parcels%numparcels_local
       z=state%parcels%z(i)
       if (z<zmin) then
         state%parcels%z(i)=zmin+(zmin-z)
       else if (z>zmax) then
         state%parcels%z(i)=zmax-(z-zmax)
       endif         
     enddo

   end subroutine



   subroutine get_buffer_ptr(ptr,dir)
     real(kind=DEFAULT_PRECISION), dimension(:,:), pointer, intent(out) :: ptr
     integer, intent(in) :: dir

     if (dir .eq. N) then
       ptr => N_buff
     else if (dir .eq. NE) then
       ptr => NE_buff
     else if (dir .eq. E) then
       ptr => E_buff
     else if (dir .eq. SE) then
       ptr => SE_buff
     else if (dir .eq. S) then
       ptr => S_buff
     else if (dir .eq. SW) then
       ptr => SW_buff
     else if (dir .eq. W) then
       ptr => W_buff
     else if (dir .eq. NW) then
       ptr => NW_buff
     endif

   end subroutine






   subroutine populate_send_buffer(state,dir,npars)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir
     integer, intent(in) :: npars

     integer(kind=PARCEL_INTEGER), allocatable, dimension(:), save :: myindex

     real(kind=DEFAULT_PRECISION), dimension(:,:), pointer, save :: buff
     real(kind=DEFAULT_PRECISION), save :: xshift, yshift

     integer(kind=PARCEL_INTEGER) :: istart, i, num

     integer :: q

     !$OMP SINGLE

     xshift=0
     yshift=0

     !we want to create a new index (myindex) which says where from the main parcel list
     ! (state%parcels...) each parcel in the buffer comes from


     !print *, "npars=",npars

     allocate(myindex(npars))
     istart=1
     do num=1,npars
       do i=istart,state%parcels%numparcels_local
         if (index(i) .eq. dir) then
           myindex(num) = i
           istart=i+1
           exit
         endif
       enddo
     enddo

     !get the pointer to the buffer array
     call get_buffer_ptr(buff,dir)

     !if our process is on the edge of the domain and the send direction sends
     !off the end then we need to adjust the parcel positions so they wrap around
     if ((dir .eq. E ) .or. (dir .eq. NE) .or. (dir .eq. SE)) then
       if (state%parallel%my_coords(3) .eq. state%parallel%dim_sizes(3)-1) then
         xshift = state%global_grid%bottom(3)-state%global_grid%top(3)-&
         state%global_grid%resolution(3)
       endif
     else if ((dir .eq. W ) .or. (dir .eq. NW) .or. (dir .eq. SW)) then
       if (state%parallel%my_coords(3) .eq. 0) then
         xshift = state%global_grid%top(3)-state%global_grid%bottom(3)+&
         state%global_grid%resolution(3)
       endif
     endif
     if ((dir .eq. N ) .or. (dir .eq. NE) .or. (dir .eq. NW)) then
       if (state%parallel%my_coords(2) .eq. state%parallel%dim_sizes(2)-1) then
         yshift = state%global_grid%bottom(2)-state%global_grid%top(2)-&
         state%global_grid%resolution(2)
       endif
     else if ((dir .eq. S ) .or. (dir .eq. SW) .or. (dir .eq. SE)) then
       if (state%parallel%my_coords(2) .eq. 0) then
         yshift = state%global_grid%top(2)-state%global_grid%bottom(2)+&
         state%global_grid%resolution(2)
       endif
     endif

     !$OMP END SINGLE

     !now we want to copy data in
     !$OMP WORKSHARE
     buff(:,X_INDEX) = state%parcels%x(myindex) + xshift
     buff(:,Y_INDEX) = state%parcels%y(myindex) + yshift
     buff(:,Z_INDEX) = state%parcels%z(myindex)
     buff(:,P_INDEX) = state%parcels%p(myindex)
     buff(:,Q_INDEX) = state%parcels%q(myindex)
     buff(:,R_INDEX) = state%parcels%r(myindex)
     buff(:,DXDT_INDEX) = state%parcels%dxdt(myindex)
     buff(:,DYDT_INDEX) = state%parcels%dydt(myindex)
     buff(:,DZDT_INDEX) = state%parcels%dzdt(myindex)
     buff(:,DPDT_INDEX) = state%parcels%dpdt(myindex)
     buff(:,DQDT_INDEX) = state%parcels%dqdt(myindex)
     buff(:,DRDT_INDEX) = state%parcels%drdt(myindex)
     buff(:,H_INDEX) = state%parcels%h(myindex)
     buff(:,B_INDEX) = state%parcels%b(myindex)
     buff(:,VOL_INDEX) = state%parcels%vol(myindex)
     buff(:,STRETCH_INDEX) = state%parcels%stretch(myindex)
     buff(:,TAG_INDEX) = state%parcels%tag(myindex)
     !$OMP END WORKSHARE

     !$OMP SINGLE
     do q=1,state%parcels%qnum
       buff(:,Q_START_INDEX+(q-1)) = state%parcels%qvalues(q,myindex)
     enddo
     !$OMP END SINGLE

     q=state%parcels%qnum+1

     !$OMP SINGLE
     deallocate(myindex)
     !$OMP END SINGLE

   end subroutine






   subroutine send_buffers(state,dir,request,count)
     type(model_state_type), intent(inout) :: state
     integer, intent(in) :: dir
     integer, intent(out) :: request
     integer, intent(in) :: count
     !integer:: dest

     real(kind=DEFAULT_PRECISION), dimension(:,:), pointer :: buff

     call get_buffer_ptr(buff,dir)

     if (count .eq. 0) then !so as to not go out of bounds we send a dummy variable if the count is 0
       call MPI_ISend(count,&
                      count*(state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum),&
                      PRECISION_TYPE,&
                      destinations(dir),&
                      dir,&
                      state%parallel%monc_communicator,&
                      request,&
                      ierr)
     else
       call MPI_ISend(buff(1,1),&
                      count*(state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum),&
                      PRECISION_TYPE,&
                      destinations(dir),&
                      dir,&
                      state%parallel%monc_communicator,&
                      request,&
                      ierr)
     endif

    ! print *, "sent in direction", dir

    end subroutine





    subroutine check_for_message(state,dir,nrecv,src)
      type(model_state_type), intent(in) :: state
      integer, intent(out) :: dir, nrecv, src
      integer :: status(MPI_STATUS_SIZE)

      call MPI_probe(MPI_ANY_SOURCE,&
                     MPI_ANY_TAG,&
                     state%parallel%monc_communicator,&
                     status,&
                     ierr)

      src=status(MPI_SOURCE)
      dir=status(MPI_TAG)

      call MPI_get_count(status,PRECISION_TYPE,nrecv,ierr)

      nrecv=nrecv/(state%parcels%n_properties+state%parcels%n_rk+state%parcels%qnum)

      !print *, "message: src, tag, size=", src, dir, nrecv

    end subroutine




    subroutine unpack_parcels(state,buff,nrecv,lastparcel,n_initial)
      type(model_state_type), intent(inout) :: state
      real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: buff
      integer, intent(in) :: nrecv
      integer(kind=PARCEL_INTEGER), intent(inout) :: lastparcel !the last parcel we replaced +1
      integer(kind=PARCEL_INTEGER), intent(in) :: n_initial

      integer(kind=PARCEL_INTEGER), dimension(:), allocatable, save :: myindex
      integer(kind=PARCEL_INTEGER) :: istart, i, num
      integer :: q



      ! we need to create an index (myindex) saying where each parcel in the buffer will go
      !to do this we look through the main index to see where parcels have been removed and
      !arrange to put the buffer parcels in these places

      !$OMP SINGLE
      call timer_start(handle_unpack)
      allocate(myindex(nrecv))

      istart=lastparcel
      do num=1,nrecv
        do i=istart,state%parcels%maxparcels_local
          if (i .gt. n_initial .or. index(i) .ne. 0) then !this parcel location is free to be overwritten
          !if (index(i) .ne. 0) then !this parcel location is free to be
            myindex(num) = i
            if (i .gt. state%parcels%maxparcels_local) then
              print *, "Error: maximum parcel count reached (parcel haloswap)"
              error stop "Maxparcels reached"
            endif
            index(i) = 0
            istart=i+1
            exit
          endif
        enddo
      enddo

      lastparcel=istart

      !$OMP END SINGLE

      !we now copy from the buffer to the parcels
      !$OMP WORKSHARE
      state%parcels%x(myindex) = buff(:,X_INDEX)
      state%parcels%y(myindex) = buff(:,Y_INDEX)
      state%parcels%z(myindex) = buff(:,Z_INDEX)
      state%parcels%p(myindex) = buff(:,P_INDEX)
      state%parcels%q(myindex) = buff(:,Q_INDEX)
      state%parcels%r(myindex) = buff(:,R_INDEX)
      state%parcels%dxdt(myindex) = buff(:,DXDT_INDEX)
      state%parcels%dydt(myindex) = buff(:,DYDT_INDEX)
      state%parcels%dzdt(myindex) = buff(:,DZDT_INDEX)
      state%parcels%dpdt(myindex) = buff(:,DPDT_INDEX)
      state%parcels%dqdt(myindex) = buff(:,DQDT_INDEX)
      state%parcels%drdt(myindex) = buff(:,DRDT_INDEX)
      state%parcels%h(myindex) = buff(:,H_INDEX)
      state%parcels%b(myindex) = buff(:,B_INDEX)
      state%parcels%vol(myindex) = buff(:,VOL_INDEX)
      state%parcels%stretch(myindex) = buff(:,STRETCH_INDEX)
      state%parcels%tag(myindex) = buff(:,TAG_INDEX)
      !$OMP END WORKSHARE

      !$OMP SINGLE
      do q=1,state%parcels%qnum
        state%parcels%qvalues(q,myindex) =  buff(:,Q_START_INDEX+(q-1))
      enddo
      !$OMP END SINGLE

      q=state%parcels%qnum+1

      !$OMP SINGLE
      deallocate(myindex)
      call timer_stop(handle_unpack)
      !$OMP END SINGLE

    end subroutine





    subroutine backfill(state,index,ntofill)
      type(model_state_type), intent(inout) :: state
      integer, dimension(:), intent(inout) :: index
      integer(kind=PARCEL_INTEGER), intent(in) :: ntofill
      integer, save :: nswap

      integer(kind=PARCEL_INTEGER), dimension(:), allocatable, save :: from, to
      integer(kind=PARCEL_INTEGER) :: i, num
      integer :: q
      !$OMP SINGLE
      !count number of existing parcels in region of array that will be removed
      !this is the number of parcels that needs to be swapped
      nswap=0
      do i=state%parcels%numparcels_local-ntofill+1,state%parcels%numparcels_local
        if (index(i) .eq. 0) nswap=nswap+1
      enddo

      !$OMP END SINGLE


      if (nswap .gt. 0) then
        !$OMP SINGLE
        allocate(to(nswap))
        allocate(from(nswap))

        !generate "to" array - look for gaps in parcels from the beginning of the array
        num=1
        do i=1,state%parcels%numparcels_local-ntofill
          if (index(i) .ne. 0) then
            to(num) = i
            num=num+1
            if (num .gt. nswap) exit
          endif
        enddo

        !if (num .ne. nswap+1) error stop "We have less parcels to swap into than we intended"

        !generate "from" array - look for esisting parcels at the end of the array
        num=1
        do i=state%parcels%numparcels_local,state%parcels%numparcels_local-ntofill+1,-1
          if (index(i) .eq. 0) then
            from(num) = i
            num=num+1
            if (num .gt. nswap) exit
          endif
        enddo

        !$OMP END SINGLE

        !if (num .ne. nswap+1) error stop "We have less parcels to swap from than we intended"


        !now we can copy parcels at "from" to parcels at "to"
        !$OMP WORKSHARE
        state%parcels%x(to) = state%parcels%x(from)
        state%parcels%y(to) = state%parcels%y(from)
        state%parcels%z(to) = state%parcels%z(from)
        state%parcels%p(to) = state%parcels%p(from)
        state%parcels%q(to) = state%parcels%q(from)
        state%parcels%r(to) = state%parcels%r(from)
        state%parcels%dxdt(to) = state%parcels%dxdt(from)
        state%parcels%dydt(to) = state%parcels%dydt(from)
        state%parcels%dzdt(to) = state%parcels%dzdt(from)
        state%parcels%dpdt(to) = state%parcels%dpdt(from)
        state%parcels%dqdt(to) = state%parcels%dqdt(from)
        state%parcels%drdt(to) = state%parcels%drdt(from)
        state%parcels%h(to) = state%parcels%h(from)
        state%parcels%b(to) = state%parcels%b(from)
        state%parcels%vol(to) = state%parcels%vol(from)
        state%parcels%stretch(to) = state%parcels%stretch(from)
        state%parcels%tag(to) = state%parcels%tag(from)
        !$OMP END WORKSHARE

        !$OMP SINGLE
        do q=1,state%parcels%qnum
          state%parcels%qvalues(q,to) =  state%parcels%qvalues(q,from)
        enddo
        !$OMP END SINGLE

        !$OMP SINGLE
        deallocate(to)
        deallocate(from)
        !$OMP END SINGLE

      endif

    end subroutine



end module
