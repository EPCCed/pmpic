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

  implicit none

  real(kind=DEFAULT_PRECISION) :: originaldt
  real(kind=DEFAULT_PRECISION), parameter :: cfl=0.1
  integer :: ierr

contains

  type(component_descriptor_type) function euler_integrator_get_descriptor()
    euler_integrator_get_descriptor%name="euler_integrator"
    euler_integrator_get_descriptor%version=0.1
    euler_integrator_get_descriptor%initialisation=>initialisation_callback
    euler_integrator_get_descriptor%timestep=>timestep_callback
    euler_integrator_get_descriptor%finalisation=>finalisation_callback
  end function euler_integrator_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print *, "In Euler Integrator Initialisation"

    originaldt= options_get_real(state%options_database,"dtm")

    state%dtm = originaldt

    print *, "Starting dt=",originaldt

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: umax, vmax, wmax, maxdt, maxdtglobal, dt
    integer(kind=PARCEL_INTEGER) :: nparcels, n

    nparcels=state%parcels%numparcels_local

    !look at velocities - get maximum then from that determine constraint on dt
    umax=maxval(abs(state%parcels%dxdt(1:nparcels)))
    vmax=maxval(abs(state%parcels%dydt(1:nparcels)))
    wmax=maxval(abs(state%parcels%dzdt(1:nparcels)))

    if (umax .eq. 0.) umax=1.e-10
    if (vmax .eq. 0.) vmax=1.e-10
    if (wmax .eq. 0.) wmax=1.e-10

    maxdt = minval( (/ dx/umax, dy/vmax, dz/wmax /) )*cfl

    if (nparcels .eq. 0) maxdt = originaldt

    call MPI_Allreduce(sendbuf=maxdt,&
                       recvbuf=maxdtglobal,&
                       count=1,&
                       datatype=PRECISION_TYPE,&
                       op=MPI_MIN,&
                       comm=state%parallel%monc_communicator,&
                       ierror=ierr)


    dt = minval( (/ maxdtglobal, originaldt/) )

    if (state%parallel%my_rank .eq. 0) print*, "t=", state%time," maxdt=",maxdtglobal, " New dt=",dt
    !print*, "vmax=",vmax

    state%dtm = dt
    state%dtm_new = dt

    !$OMP PARALLEL DO
    do n=1,nparcels
      state%parcels%x(n) = state%parcels%x(n) + state%parcels%dxdt(n) * dt
      state%parcels%y(n) = state%parcels%y(n) + state%parcels%dydt(n) * dt
      state%parcels%z(n) = state%parcels%z(n) + state%parcels%dzdt(n) * dt
    enddo
    !$OMP END PARALLEL DO

    call parcel_haloswap(state)

  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print *, "Integrator finalisation - nothing to see here"
  end subroutine

end module
