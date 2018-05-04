!Component that prescribes parcel velocities

!current profiles in use are:
!1) shear flow (u,v,w) = (0,tanh(x-x0),0)
!2) cylindrical flow (u,v,w) = r*exp(-(r-r0)^2/l0^2)(cos(theta), -sin(theta), 0)

module prescribed_parcel_velocity_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, meandz

  implicit none

  integer :: profile_type
  real(kind=DEFAULT_PRECISION) :: x0, y0, r0

contains

  type(component_descriptor_type) function prescribed_parcel_velocity_get_descriptor()
    prescribed_parcel_velocity_get_descriptor%name="prescribed_parcel_velocity"
    prescribed_parcel_velocity_get_descriptor%version=0.1
    prescribed_parcel_velocity_get_descriptor%initialisation=>initialisation_callback
    prescribed_parcel_velocity_get_descriptor%timestep=>timestep_callback
    prescribed_parcel_velocity_get_descriptor%finalisation=>finalisation_callback
  end function prescribed_parcel_velocity_get_descriptor



  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer :: n

    profile_type=options_get_integer(state%options_database,"velocity_profile")
    !x0=options_get_real(state%options_database,"x0")

    !centre of grids in y and x
    x0=(state%global_grid%top(3)-state%global_grid%bottom(3))/2. + state%global_grid%bottom(3)
    y0=(state%global_grid%top(2)-state%global_grid%bottom(2))/2. + state%global_grid%bottom(2)

    r0=(state%global_grid%top(3)-state%global_grid%bottom(3))/8.


    !print*, "Velocity option read in as", profile_type, x0

    !we now need to tag parcels according to the profile
    !at present we use a basic tag where:
    ! y<x0 : tag=0
    !y >= x0 : tag=1



    !$OMP PARALLEL do
    do n=1,state%parcels%numparcels_local
      state%parcels%tag(n) = 0
      if (state%parcels%y(n) .gt. y0) then
        state%parcels%tag(n) = state%parcels%tag(n) + 2
      endif
      if (state%parcels%x(n) .gt. x0) then
        state%parcels%tag(n) = state%parcels%tag(n) + 1
      endif
    enddo
    !$OMP END PARALLEL do

    print*, "Tagged parcels"


  end subroutine initialisation_callback


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer :: n
    real(kind=DEFAULT_PRECISION) :: x, y, z, v0, theta, r2

    if (profile_type .eq. 1) then

      !$OMP PARALLEL DO PRIVATE(x)
      do n=1,state%parcels%numparcels_local
        x = state%parcels%x(n)
        !y = state%parcels%y(n)
        !z = state%parcels%z(n)

        state%parcels%dxdt(n) = 0.
        state%parcels%dydt(n) = tanh((x-x0)/r0)
        state%parcels%dzdt(n) = 0.
      enddo
      !$OMP END PARALLEL DO

    else if (profile_type .eq. 2) then

      !$OMP PARALLEL DO PRIVATE(x,y,v0,theta,r2)
      do n=1,state%parcels%numparcels_local
        x = state%parcels%x(n)
        y = state%parcels%y(n)
        !z = state%parcels%z(n)

        theta=atan(y-y0,x-x0)

        r2 = ((x-x0)*(x-x0) + (y-y0)*(y-y0))

        v0=exp( -r2/r0/r0/2 )/x0*10

        state%parcels%dxdt(n) = v0*sqrt(r2)*sin(theta)
        state%parcels%dydt(n) = -v0*sqrt(r2)*cos(theta)
        state%parcels%dzdt(n) = 0.
      enddo
      !$OMP END PARALLEL DO

    else
      error stop "invalid velocity option"
    endif

    print*, "velocities set"


  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print*, "Finalising velocity - nothing to see here"

  end subroutine

end module
