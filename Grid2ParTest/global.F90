!module containing global variables etc
module global_mod

    use MPI
    use OMP_LIB

    implicit none

    !parameters for the computational domain:

    !base grid resolution (grid can be integer multiples greater or smaller than this)
    integer :: nxbase, nybase, nzbase


    ! size of the domain in units
    double precision :: xmin, xmax, ymin, ymax, zmin, zmax

    double precision :: tp2g, tg2p

    integer :: ierror

contains



end module
