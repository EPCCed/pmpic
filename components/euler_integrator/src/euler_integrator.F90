!A basic Euler integrator component
module euler_integrator_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, meandz

  implicit none

  real(kind=DEFAULT_PRECISION) :: originaldt
  real(kind=DEFAULT_PRECISION), parameter :: cfl=0.1

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

    originaldt = state%dtm

    print *, "Starting dt=",originaldt

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: umax, vmax, wmax, maxdt, dt
    integer :: nparcels, n

    nparcels=state%parcels%numparcels_local

    !look at velocities - get maximum then from that determine constraint on dt
    umax=maxval(abs(state%parcels%dxdt(1:nparcels)))
    vmax=maxval(abs(state%parcels%dydt(1:nparcels)))
    wmax=maxval(abs(state%parcels%dzdt(1:nparcels)))

    if (umax .eq. 0.) umax=1.e-10
    if (vmax .eq. 0.) vmax=1.e-10
    if (wmax .eq. 0.) wmax=1.e-10

    maxdt = minval( (/ dx/umax, dy/vmax, meandz/wmax /) )*cfl
    dt = minval( (/ maxdt, originaldt/) )

    print*, "t=", state%time," maxdt=",maxdt, " New dt=",dt
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

  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print *, "Integrator finalisation - nothing to see here"
  end subroutine

end module
