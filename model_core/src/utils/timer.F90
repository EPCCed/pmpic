!Allows parts of MONC/MPIC to be timed, and provides statistics at the end of the run
module timer_mod

  use MPI, only: MPI_Wtime, MPI_Gather, MPI_DOUBLE_PRECISION
  use state_mod, only: model_state_type

  implicit none

  integer, parameter :: max_entries = 100
  integer :: num !number of active entries
  integer :: idum

  type timing_data_type
    character (len=20) :: name
    integer :: n
    double precision :: tstart
    double precision :: t
    double precision :: dt
    double precision :: dt2 ! the sum of dt^2 (to calculate standard deviation)
    double precision :: max, min
    double precision :: mean
    integer :: unit = -1
    logical :: paused
    logical :: running
  end type

  type(timing_data_type), allocatable, dimension(:) :: timings

contains

  !initialise the timing system
  subroutine init_timing(state)
    type(model_state_type), intent(in) :: state
    allocate(timings(max_entries))
    num=0

    if (state%parallel%my_rank .eq. 0) then
      write(*,"('Self-timing initialised with a maximum of',i4,' tracable routines')") max_entries
    endif
  end subroutine


  !registers a routine with 'name' with the timing system, and returns a handle
  subroutine register_routine_for_timing(name,handle,state)
    character (len=*), intent(in) :: name
    integer, intent(out) :: handle
    type(model_state_type), intent(in) :: state

    integer :: length

    if (num .ge. max_entries) then
      error stop "maximum number of tracable routines reached"
    endif

    handle=num+1
    num=num+1

    length=len(name)

    if (length .ge. 20) then
      print *, "Warning, routine '",name,"' name will be truncated"
      length=19
    endif

    write(timings(handle)%name,*) name(1:length)
    timings(handle)%tstart=0.
    timings(handle)%dt=0.
    timings(handle)%t = 0.
    timings(handle)%dt2 = 0.
    timings(handle)%max=0.
    timings(handle)%min = 0.
    timings(handle)%paused = .false.
    timings(handle)%running = .false.
    timings(handle)%n = 0

    if (state%parallel%my_rank .eq. 0) then
      write(*,"('Registered routine `',a,'` for timing with handle =',i4)") &
       timings(handle)%name,handle
    endif
  end subroutine

  subroutine timer_start(handle)
    integer, intent(in) :: handle

    if (timings(handle)%running) error stop "Timer is already running"

    timings(handle)%tstart = MPI_Wtime()
    timings(handle)%running = .true.
    timings(handle)%dt = 0.
  end subroutine


  subroutine timer_pause(handle)
    integer, intent(in) :: handle

    if (timings(handle)%paused) error stop "Timer is already paused"

    timings(handle)%dt = timings(handle)%dt + (MPI_Wtime() - timings(handle)%tstart)
    timings(handle)%paused = .true.
  end subroutine

  subroutine timer_resume(handle)
    integer, intent(in) :: handle

    if (.not. timings(handle)%paused) error stop "Timer is not paused - cannot resume"

    timings(handle)%paused = .false.
    timings(handle)%tstart = MPI_Wtime()
  end subroutine


  subroutine timer_stop(handle)
    integer, intent(in) :: handle

    if (timings(handle)%paused) error stop "cannot stop as it is paused"

    timings(handle)%dt = timings(handle)%dt + (MPI_Wtime()-timings(handle)%tstart)
    timings(handle)%t = timings(handle)%t + timings(handle)%dt

    timings(handle)%dt2 = timings(handle)%dt2 + timings(handle)%dt*timings(handle)%dt



    if (timings(handle)%n .eq. 1) then
      timings(handle)%max = timings(handle)%dt
      timings(handle)%min = timings(handle)%dt
    else
      if (timings(handle)%dt .gt. timings(handle)%max) timings(handle)%max = timings(handle)%dt
      if (timings(handle)%dt .lt. timings(handle)%min) timings(handle)%min = timings(handle)%dt
    endif

    timings(handle)%running = .false.

    timings(handle)%dt = 0.
    timings(handle)%n = timings(handle)%n + 1

  end subroutine

  subroutine finalize_timing(state)
     type(model_state_type), intent(in) :: state
     integer :: n
     double precision :: mean, percentage, imbalance, stddev, mxtime
     double precision, allocatable, dimension(:) :: maxtime, meantime, load_imb, tdiff
     integer, allocatable, dimension(:) :: fastest, slowest
     integer, allocatable, dimension(:) :: index
     double precision, allocatable, dimension(:) :: values
     double precision, allocatable, dimension(:,:) :: rank_times
     integer :: i

     mxtime=timings(1)%t

     allocate(index(num), values(num))

     do n=1,num
       values(n) = timings(n)%t
     enddo

     call sort_values(values,index,num)

     !now print stats from rank 0

     if (state%parallel%my_rank .eq. 0) then

       print *, ""
       print *, "Summary of times for the master rank (sorted from longest to shortest time):"

         write(*,'(" ------------------------------------------------------------------------------ ")')
         write(*,'("|     Function       |   #    |  Total   |  % of  |  Mean   |  Max    | stddev |")')
         write(*,'("|       Name         | Calls  |  Time    |  Time  |Duration |Duration |    %   |")')
         write(*,'("|------------------------------------------------------------------------------|")')

       do i=1,num

         n=index(i)

         if (timings(n)%n .gt. 0) then
           mean = timings(n)%t/timings(n)%n
           stddev =sqrt( timings(n)%dt2 / timings(n)%n - mean*mean)
           !imbalance = (timings(n)%max - mean)/mean *100.
           imbalance = (stddev/mean)*100
         else
           mean=0.
           stddev=0.
           imbalance=0.
         endif
         percentage= timings(n)%t/mxtime*100

         !print *, timings(n)%max, mean, timings(n)%max-mean

         write(*,"('|',a,'|', i8,'|', f9.3,'s', '| ', f6.2,'%','|', f8.3,'s|',f8.3, 's|', f7.2,'%|')")  &
         timings(n)%name, timings(n)%n,timings(n)%t, percentage , mean, timings(n)%max ,imbalance
      enddo
       write(*,'(" ------------------------------------------------------------------------------  ")')
     endif


     if (state%parallel%processes .gt. 1) then
       !if (state%parallel%my_rank .eq. 0) print *, "We have more than one rank!"

       allocate(rank_times(state%parallel%processes, num))
       allocate(slowest(num), fastest(num), meantime(num), maxtime(num), tdiff(num), load_imb(num))

       if (state%parallel%my_rank .eq. 0) then

         print *, ""
         print *, "Load imbalances over processes (sorted by difference in time):"

         write(*,'(" ------------------------------------------------------------------------------ ")')
         write(*,'("|     Function       |   Mean   |   Max    |Difference |    %    | Slow | Fast |")')
         write(*,'("|       Name         |   Time   |   Time   |(max-mean) |Imbalance| Proc | Proc |")')
         write(*,'("|------------------------------------------------------------------------------|")')
       endif

       !get all the data
       do i=1,num

         !n=index(i)
         call MPI_Gather(sendbuf=timings(i)%t,&
                         sendcount=1,&
                         sendtype=MPI_DOUBLE_PRECISION,&
                         recvbuf=rank_times(1,i),&
                         recvcount=1,&
                         recvtype=MPI_DOUBLE_PRECISION,&
                         root=0,&
                         comm=state%parallel%monc_communicator,&
                         ierror=idum)



         meantime(i)=sum(rank_times(:,i))/state%parallel%processes

         if (meantime(i) .gt. 0) then

           maxtime(i)=maxval(rank_times(:,i))

           tdiff(i)=maxtime(i)-meantime(i)

           load_imb(i) = (maxtime(i)-meantime(i))/meantime(i) * 100

           slowest(i)=maxloc(rank_times(:,i), dim=1)-1
           fastest(i)=minloc(rank_times(:,i), dim=1)-1
         else
           load_imb(i) = 0.
         endif


       enddo

       values=tdiff
       call sort_values(values,index,num)

       do i=1,num

         n=index(i)

          if (meantime(n) .gt. 0) then
            if (state%parallel%my_rank .eq. 0) then
              write(*,"('|',a,'|', f9.3,'s|', f9.3,'s', '| 'f9.5,'s| ', f7.2,'%','|', i6,'|',i6, '|')")  &
              timings(n)%name, meantime(n), maxtime(n) ,tdiff(n), load_imb(n), slowest(n) ,fastest(n)
            endif
          else

            maxtime=0.
            tdiff=0.

            if (state%parallel%my_rank .eq. 0) then
              write(*,"('|',a,'|', f9.3,'s|', f9.3,'s', '| 'f9.3,'s|   N/A   |  N/A |  N/A |')")  &
              timings(n)%name, meantime(n), maxtime(n), tdiff(n)
            endif

          endif

       enddo

       if (state%parallel%my_rank .eq. 0) then
          write(*,'(" ------------------------------------------------------------------------------  ")')
          print *, ""
       endif

       deallocate(rank_times)



     endif


  end subroutine


  !return index of sorted values
  ! as the number is small (~100) we just use bubble sort as it's not worth the effort to use
  ! a fancier sort algorithm such as merge sort
  subroutine sort_values(values,index,n)
    double precision, intent(inout) :: values(:)
    integer, intent(inout) :: index(:)
    integer, intent(in) :: n !number of items in the array

    integer :: i
    integer :: check
    double precision :: valdum
    integer :: idxdum

    do i=1,n
      index(i) = i
    enddo

    check=0
    do while (check .lt. num-1)
      check = 0
      do i=1,n-1
        if (values(i) .lt. values(i+1)) then
          valdum = values(i)
          idxdum = index(i)

          values(i) = values(i+1)
          values(i+1) = valdum

          index(i) = index(i+1)
          index(i+1) = idxdum
        else
          check=check+1
        endif
      enddo
    enddo

    !do i=1,num
  !    print *, values(i), index(i)
  !  enddo

  end subroutine






end module
