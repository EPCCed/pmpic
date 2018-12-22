!sets up some parcels and grids then tests grid2par
module grid2partest_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER
  use prognostics_mod, only: prognostic_field_type
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use parcel_interpolation_mod, only: cache_parcel_interp_weights, grid2par, par2grid, x_coords, y_coords, z_coords
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use MPI
  use parcel_haloswap_mod, only: parcel_haloswap
  use timer_mod

  implicit none

  integer :: nx, ny, nz
  real(kind=DEFAULT_PRECISION) :: dx, dy, dz!, xmin, xmax, ymin, ymax, zmin,zmax

  integer :: ierr

  integer :: handle_g2p, handle_p2g

contains

  type(component_descriptor_type) function grid2partest_get_descriptor()
    grid2partest_get_descriptor%name="grid2partest"
    grid2partest_get_descriptor%version=0.1
    grid2partest_get_descriptor%initialisation=>initialisation_callback
    grid2partest_get_descriptor%timestep=>timestep_callback
    grid2partest_get_descriptor%finalisation=>finalisation_callback
  end function grid2partest_get_descriptor


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    if (current_state%parallel%my_rank .eq. 0) then
      print *, ""
      print*, "grid2partest initialisation"
    endif

    !set up grid

    current_state%rksteps = 1

    dx = current_state%global_grid%resolution(3)
    dy = current_state%global_grid%resolution(2)
    dz = current_state%global_grid%resolution(1)


    call setup_parcels(current_state)

    call setup_grid(current_state)

    call register_routine_for_timing("verify_grid2par",handle_g2p,current_state)
    call register_routine_for_timing("verify_par2grid",handle_p2g,current_state)


  end subroutine



  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    call MPI_Barrier(current_state%parallel%monc_communicator,ierr)

    if (current_state%parallel%my_rank .eq. 0) then
      print *, ""
      print*, "grid2partest:"
    endif

    !call parcel_haloswap(current_state)

    !cache interpolation weights for parcels
    call cache_parcel_interp_weights(current_state)

    if (current_state%parallel%my_rank .eq. 0) print*, "Testing in x direction"
    call grid2par(current_state,current_state%u,current_state%parcels%dxdt)
    call check_parcels(current_state,current_state%parcels%dxdt,current_state%parcels%x)
    call par2grid(current_state,current_state%parcels%p,current_state%p,current_state%parcels%vol)
    call check_grid(current_state,current_state%p,current_state%u)

    if (current_state%parallel%my_rank .eq. 0) print *, "Testing in y direction"
    call grid2par(current_state,current_state%v,current_state%parcels%dydt)
    call check_parcels(current_state,current_state%parcels%dydt,current_state%parcels%y)

    call par2grid(current_state,current_state%parcels%q,current_state%q,current_state%parcels%vol)
    call check_grid(current_state,current_state%q,current_state%v)

    if (current_state%parallel%my_rank .eq. 0) print *, "Testing in Z direction"
    call grid2par(current_state,current_state%w,current_state%parcels%dzdt)
    call check_parcels(current_state,current_state%parcels%dzdt,current_state%parcels%z)

    call par2grid(current_state,current_state%parcels%r,current_state%r,current_state%parcels%vol)
    call check_grid(current_state,current_state%r,current_state%w)

    call MPI_Barrier(current_state%parallel%monc_communicator,ierr)

    if (current_state%parallel%my_rank .eq. 0) then
      print *, "grid2partest finished"
      print *, ""
    endif

    call MPI_Barrier(current_state%parallel%monc_communicator,ierr)


  end subroutine

  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    !print*, "grid2partest finalisation"
    deallocate(current_state%u%data, current_state%v%data, current_state%w%data)


  end subroutine

  ! sets up p=x, q=y, r=z
  subroutine setup_parcels(state)
    type(model_state_type), intent(inout) :: state

    integer(kind=PARCEL_INTEGER) :: n

    do n=1,state%parcels%numparcels_local
      state%parcels%p(n)=state%parcels%x(n)
      state%parcels%q(n)=state%parcels%y(n)
      state%parcels%r(n)=state%parcels%z(n)
    enddo

    if (state%parallel%my_rank .eq. 0) print *, "grid2partest: initialised parcels"

  end subroutine


  !sets up u=x, v=y, w=z
  subroutine setup_grid(state)
    type(model_state_type), intent(inout) :: state
    integer :: i, j, k

    !number of cells in local grid (cells belonging to grid and halos)
    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)


    ! set values for each variable:
    ! u = x
    ! v = y
    ! w = z
    do i=1,nx
      do j=1,ny
        do k=1,nz
          state%u%data(k,j,i) = x_coords(i)
          state%v%data(k,j,i) = y_coords(j)
          state%w%data(k,j,i) = z_coords(k)
        enddo
      enddo
    enddo

    if (state%parallel%my_rank .eq. 0) print *, "grid2partest: initialised grids"

  end subroutine


  !compares interpolated values with reference values
  subroutine check_parcels(state, values, reference)
    type(model_state_type), intent(in) :: state
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: values, reference
    real(kind=DEFAULT_PRECISION), parameter :: tol=1.e-9
    real(kind=DEFAULT_PRECISION) :: diff
    integer(kind=PARCEL_INTEGER) :: n, nparcels

    call timer_start(handle_g2p)

    nparcels=state%parcels%numparcels_local

    do n=1,nparcels
      diff =abs(values(n)-reference(n))
      if (diff .gt. tol) then
        !answer will be wrong at up, down left and right extreme boundaries of grid because test profile is
        !linear so it can't deal with the periodicity of the grid
        if ((state%parcels%x(n) .gt. state%global_grid%bottom(3)+dx) .and. &
            (state%parcels%x(n) .lt. state%global_grid%top(3)-dx)) then
            if ((state%parcels%y(n) .gt. state%global_grid%bottom(2)+dy) .and. &
                (state%parcels%y(n) .lt. state%global_grid%top(2)-dy)) then
                    print *, n, values(n), reference(n), diff
                    print*, state%parallel%my_rank
                    error stop "Wrong answer grid2par"
            endif
        endif
      endif
    enddo

    call MPI_Barrier(state%parallel%monc_communicator,ierr)

    if (state%parallel%my_rank .eq. 0) print*, "grid2par result: verified"

    call timer_stop(handle_g2p)

  end subroutine

    !compares interpolated values with reference values
  subroutine check_grid(state, values, reference)
    type(model_state_type), intent(in) :: state
    type(prognostic_field_type), intent(in) :: values, reference
    real(kind=DEFAULT_PRECISION), parameter :: tol=1.e-9
    real(kind=DEFAULT_PRECISION) :: diff

    integer :: xhalo, yhalo, zhalo
    integer :: xstart, ystart, zstart
    integer :: xstop, ystop, zstop
    integer :: i, j, k

    call timer_start(handle_p2g)

    !again, answer will be incorrect on the edges of the global grid due to periodicity
    !so we don't want to check cells on the edge of the grid

    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)

    xhalo = state%local_grid%halo_size(3)
    yhalo = state%local_grid%halo_size(2)
    zhalo = state%local_grid%halo_size(1)

    xstart=1
    xstop=nx

    if (x_coords(1) .lt. state%global_grid%bottom(3)) xstart=1+xhalo+1
    if (x_coords(nx) .gt. state%global_grid%top(3)) xstop=nx-xhalo-1

    ystart=1
    ystop=ny

    if (y_coords(1) .lt. state%global_grid%bottom(2)) ystart=1+yhalo+1
    if (y_coords(ny) .gt. state%global_grid%top(2)) ystop=ny-yhalo-1

    zstart=2
    zstop=nz-1

    do i=xstart, xstop
      do j=ystart,ystop
        do k=zstart,zstop
          diff = abs(values%data(k,j,i) - reference%data(k,j,i))
          if (diff .gt. tol ) then
            print*, "i,j,k=", i, j, k
            print*, "x, y, z=", x_coords(i), y_coords(j), z_coords(k)
            print*, "reference=", reference%data(k,j,i)
            print*, "value=", values%data(k,j,i)
            print *, "rank=", state%parallel%my_rank
            error stop "par2grid PROBLEM"
          endif
        enddo
      enddo
    enddo

    call MPI_Barrier(state%parallel%monc_communicator,ierr)

    if (state%parallel%my_rank .eq. 0) print*, "par2grid result: verified"

    call timer_stop(handle_p2g)

  end subroutine








end module
