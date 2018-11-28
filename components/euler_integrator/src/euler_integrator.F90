!A basic Euler integrator component
module euler_integrator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, dz
  use parcel_haloswap_mod, only: parcel_haloswap
  use MPI
  use timer_mod

  implicit none

  real(kind=DEFAULT_PRECISION) :: originaldt
  real(kind=DEFAULT_PRECISION), parameter :: cfl=1
  integer :: ierr
  integer :: handle

contains

  type(component_descriptor_type) function euler_integrator_get_descriptor()
    euler_integrator_get_descriptor%name="euler_integrator"
    euler_integrator_get_descriptor%version=0.1
    euler_integrator_get_descriptor%initialisation=>initialisation_callback
    euler_integrator_get_descriptor%timestep=>timestep_callback
    !euler_integrator_get_descriptor%finalisation=>finalisation_callback
  end function euler_integrator_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) print *, "In Euler Integrator Initialisation"

    originaldt= options_get_real(state%options_database,"dtm")
    state%dtmax = options_get_real(state%options_database, "dtmax")

    state%dtm = originaldt

    if (state%parallel%my_rank .eq. 0) print *, "Starting dt=",originaldt

    if (state%rksteps .ne. 0) then
      print *, "Error - another integrator component is present. ABORTING"
      stop
    endif

    call register_routine_for_timing("Euler_integrator",handle,state)
    state%rksteps = 1

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: umax, vmax, wmax, maxdt, maxdtglobal, dt
    integer(kind=PARCEL_INTEGER) :: nparcels, n

    call timer_start(handle)

    nparcels=state%parcels%numparcels_local


    dt = state%dtm
    if (state%parallel%my_rank .eq. 0) print*, "t=", state%time,"dt=",state%dtm

    !$OMP PARALLEL DO
    do n=1,nparcels
      state%parcels%x(n) = state%parcels%x(n) + state%parcels%dxdt(n) * dt
      state%parcels%y(n) = state%parcels%y(n) + state%parcels%dydt(n) * dt
      state%parcels%z(n) = state%parcels%z(n) + state%parcels%dzdt(n) * dt
      state%parcels%p(n) = state%parcels%p(n) + state%parcels%dpdt(n) * dt
      state%parcels%q(n) = state%parcels%q(n) + state%parcels%dqdt(n) * dt
      state%parcels%r(n) = state%parcels%r(n) + state%parcels%drdt(n) * dt
    enddo
    !$OMP END PARALLEL DO

    call timer_stop(handle)

    call parcel_haloswap(state)

  end subroutine



end module
