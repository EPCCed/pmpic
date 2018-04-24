module parcel_mod
  use datadefn_mod, only : DEFAULT_PRECISION

  implicit none

private

  type, public :: parcel_type

    integer :: numparcels_global !number of active parcels (globally)
    integer :: numparcels_local !number of active parcels (belonging to process)
    integer :: maxparcels_global !maximum number of parcels globally
    integer :: maxparcels_local !maximum number of parcels (beloning to process)

    real (kind=DEFAULT_PRECISION), allocatable, dimension(:) :: &
                x, y, z, & !positions
                r, s, t, & !vorticities
                u, v, w, & !velocities
                q, b, vol  !humidity, buoyancy, volume
  end type parcel_type

end module
