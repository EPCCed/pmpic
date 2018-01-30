!module containing global variables etc
module global_mod

    implicit none

    !parameters for the computational domain:

    !base grid resolution (grid can be integer multiples greater or smaller than this)
    integer :: nxbase, nybase, nzbase


    ! size of the domain in units
    double precision :: xmin, xmax, ymin, ymax, zmin, zmax

contains



end module
