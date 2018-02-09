!module defining parcel types and subroutines
module parcel_mod
    implicit none

    integer, parameter :: maxparcels = 20000000
    integer :: nparcels

    !parcel structure description
    type parcel
        double precision :: x, y, z
        double precision :: u,v,w
        double precision :: r,s,t
        double precision :: bl,thetal,q,volume
        logical :: active
        double precision :: dummyvars(10)
    end type parcel

    !array containing parcels
    type(parcel), allocatable, dimension(:) :: parcels

    !old array descriptions
    double precision, allocatable, dimension(:) :: xpos, ypos, zpos
    double precision, allocatable, dimension(:) :: uvel, vvel, wvel
    double precision, allocatable, dimension(:) :: rvort, svort, tvort
    double precision, allocatable, dimension(:) :: bl, thetal, q,volume
    logical, allocatable, dimension(:) :: active
    double precision, allocatable, dimension(:,:) :: dummyvars


end module
