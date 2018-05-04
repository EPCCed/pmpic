!A Haloswapping routine for parcels.
!At present because haloswapping isn't working properly all this does is
!apply periodic boundary conditions to the parcels in a single process
module parcelswap_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, meandz, maxx, maxy, minx, miny

  implicit none

contains

  type(component_descriptor_type) function parcelswap_get_descriptor()
    parcelswap_get_descriptor%name="parcelswap"
    parcelswap_get_descriptor%version=0.1
    parcelswap_get_descriptor%initialisation=>initialisation_callback
    parcelswap_get_descriptor%timestep=>timestep_callback
    parcelswap_get_descriptor%finalisation=>finalisation_callback
  end function parcelswap_get_descriptor

  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

  end subroutine initialisation_callback

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer :: n

    !$OMP PARALLEL DO
    do n=1,state%parcels%numparcels_local
      !basic periodic BCs
      if (state%parcels%x(n) .gt. maxx) then
        state%parcels%x(n) = state%parcels%x(n)-(maxx-minx)
      endif
      if (state%parcels%x(n) .lt. minx) then
        state%parcels%x(n) = state%parcels%x(n)+(maxx-minx)
      endif

      if (state%parcels%y(n) .gt. maxy) then
        state%parcels%y(n) = state%parcels%y(n)-(maxy-miny)
      endif
      if (state%parcels%y(n) .lt. miny) then
        state%parcels%y(n) = state%parcels%y(n)+(maxy-miny)
      endif
    enddo
    !$OMP END PARALLEL DO

  end subroutine timestep_callback

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

  end subroutine finalisation_callback

end module
