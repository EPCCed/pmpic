module parcel_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER

  implicit none


private

  type, public :: parcel_type

    integer(kind=PARCEL_INTEGER) :: numparcels_global !number of active parcels (globally)
    integer(kind=PARCEL_INTEGER) :: numparcels_local !number of active parcels (belonging to process)
    integer(kind=PARCEL_INTEGER) :: maxparcels_global !maximum number of parcels globally
    integer(kind=PARCEL_INTEGER) :: maxparcels_local !maximum number of parcels (beloning to process)
    integer :: n_properties = 17 ! number of parcel properties (not including qvalues)
    integer :: n_rk = 0 !number of active RK variables
    integer :: qnum=0 !number of q values per parcel

    real (kind=DEFAULT_PRECISION), allocatable, dimension(:) :: &
                x, y, z, & !positions
                p, q, r, & !vorticities
                dxdt, dydt, dzdt, & !velocities
                dpdt, dqdt, drdt, & !vorticity tendency
                h, b, vol, & !humidity, buoyancy, volume
                stretch, tag, & !stretch, a tag to tag parcels so they can be tracked
                btot !total humidity
    real (kind=DEFAULT_PRECISION), allocatable, dimension(:,:) :: qvalues !an arrayof various values the user may wish to add into the code

    !rk4 variables (if the rk4 integrator component isn't enabled we don't bother allocating these)
    real (kind=DEFAULT_PRECISION), allocatable, dimension(:) :: &
                xo, yo, zo, &
                xf, yf, zf, &
                po, qo, ro, &
                pf, qf, rf



  end type parcel_type

end module
