!sets up some parcels and grids then tests grid2par
module grid2partest_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use prognostics_mod, only: prognostic_field_type
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use parcel_interpolation_mod, only: cache_parcel_interp_weights, grid2par, par2grid, x_coords, y_coords, z_coords
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array

  implicit none

  integer :: nx, ny, nz
  real(kind=DEFAULT_PRECISION) :: dx, dy, dz!, xmin, xmax, ymin, ymax, zmin,zmax

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

    !print*, "grid2partest initialisation"

    !set up grid

    call setup_parcels(current_state)

    call setup_grid(current_state)


  end subroutine



  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    print *, ""
    print*, "grid2partest:"

    !cache interpolation weights for parcels
    call cache_parcel_interp_weights(current_state)

    print *, "Testing in x direction"
    call grid2par(current_state,current_state%u,current_state%parcels%u)
    call check_parcels(current_state,current_state%parcels%u,current_state%parcels%x)

    call par2grid(current_state,current_state%parcels%r,current_state%r)
    call check_grid(current_state,current_state%r,current_state%u)

    print *, "Testing in y direction"
    call grid2par(current_state,current_state%v,current_state%parcels%v)
    call check_parcels(current_state,current_state%parcels%v,current_state%parcels%y)

    call par2grid(current_state,current_state%parcels%s,current_state%s)
    call check_grid(current_state,current_state%s,current_state%v)

    print *, "Testing in Z direction"
    call grid2par(current_state,current_state%w,current_state%parcels%w)
    call check_parcels(current_state,current_state%parcels%w,current_state%parcels%z)

    call par2grid(current_state,current_state%parcels%t,current_state%t)
    call check_grid(current_state,current_state%t,current_state%w)

    print *, "grid2partest finished"
    print *, ""


  end subroutine

  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    !print*, "grid2partest finalisation"
    deallocate(current_state%u%data, current_state%v%data, current_state%w%data)


  end subroutine


  subroutine setup_parcels(state)
    type(model_state_type), intent(inout) :: state

    integer :: n

    do n=1,state%parcels%numparcels_local
      state%parcels%r(n)=state%parcels%x(n)
      state%parcels%s(n)=state%parcels%y(n)
      state%parcels%t(n)=state%parcels%z(n)
    enddo

  end subroutine



  subroutine setup_grid(state)
    type(model_state_type), intent(inout) :: state
    integer :: i, j, k

    !number of cells in local grid (cells belonging to grid and halos)
    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)


    allocate(state%u%data(nz,ny,nx))
    allocate(state%v%data(nz,ny,nx))
    allocate(state%w%data(nz,ny,nx))

    allocate(state%r%data(nz,ny,nx))
    allocate(state%s%data(nz,ny,nx))
    allocate(state%t%data(nz,ny,nx))

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

    ! xmin = x_coords(1+state%local_grid%halo_size(3))
    ! ymin = y_coords(1+state%local_grid%halo_size(2))
    ! zmin = z_coords(1+state%local_grid%halo_size(1))
    !
    ! xmax = x_coords(nx-state%local_grid%halo_size(3))
    ! ymax = y_coords(ny-state%local_grid%halo_size(2))
    ! zmax = z_coords(nz-state%local_grid%halo_size(1))


    print *, "grid2partest: initialised grids"

  end subroutine



  subroutine check_parcels(state, values, reference)
    type(model_state_type), intent(in) :: state
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: values, reference
    real(kind=DEFAULT_PRECISION), parameter :: tol=1.e-9
    real(kind=DEFAULT_PRECISION) :: diff
    integer :: n, nparcels

    nparcels=state%parcels%numparcels_local

    do n=1,nparcels
      diff =abs(values(n)-reference(n))
      if (diff .gt. tol) then
        print *, n, values(n), reference(n), diff
        error stop "Wrong answer"
      endif
    enddo

    print*, "grid2par result: verified"

  end subroutine


  subroutine check_grid(state, values, reference)
    type(model_state_type), intent(in) :: state
    type(prognostic_field_type), intent(in) :: values, reference
    real(kind=DEFAULT_PRECISION), parameter :: tol=1.e-9
    real(kind=DEFAULT_PRECISION) :: diff

    integer :: xhalo, yhalo, zhalo
    integer :: i, j, k

    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)

    xhalo = state%local_grid%halo_size(3)
    yhalo = state%local_grid%halo_size(2)
    zhalo = state%local_grid%halo_size(1)

    do i=1+xhalo+1,nx-xhalo-1
      do j=1+yhalo+1,ny-yhalo-1
        do k=1+zhalo+1,nz-zhalo-1
          diff = abs(values%data(k,j,i) - reference%data(k,j,i))
          if (diff .gt. tol ) then
            print*, "i,j,k=", i, j, k
            print*, "x, y, z=", x_coords(i), y_coords(j), z_coords(k)
            print*, "reference=", reference%data(k,j,i)
            print*, "value=", values%data(k,j,i)
            error stop "PROBLEM"
          endif
        enddo
      enddo
    enddo

    print*, "par2grid result: verified"

  end subroutine








end module
