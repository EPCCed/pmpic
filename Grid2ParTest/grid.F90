!module which contains types and subroutines related to the gridded variables
module grid_mod
    use global_mod
    implicit none

    type gridded_variable

        !number of cells in each direction
        integer :: nx, ny, nz;

        !grid spacing in each direction
        double precision :: dx, dy, dz;

        !the grid itselg
        double precision, allocatable, dimension(:,:,:) :: data;

    end type gridded_variable

    type(gridded_variable) :: ugrid, vgrid, wgrid, pgrid

contains


    !Initialises the gridded variable, and allocates the memory for the grid
    subroutine Initialise_Grid(grid,nx,ny,nz)
        implicit none
        type(gridded_variable), intent(inout) ::  grid
        integer, intent(in) :: nx, ny, nz

        allocate(grid%data(nx,ny,nz))

        grid%nx=nx
        grid%ny=ny
        grid%nz=nz

        grid%dx=dble(nx/(xmax-xmin))
        grid%dy=dble(ny/(ymax-ymin))
        grid%dz=dble(nz/(zmax-zmin))

    end subroutine

    !finalises the gridded variable (deallocates the data)
    subroutine Finalize_grid(grid)
        implicit none
        type(gridded_variable), intent(inout) :: grid

        deallocate(grid%data)

    end subroutine


    !from a set of coordinates (x,y,z) returns the lower left grid indices (i,j,k)
    subroutine Get_Grid_Coords(grid,x,y,z,i,j,k)
        type(gridded_variable), intent(in) :: grid
        double precision, intent(in) :: x,y,z
        integer, intent(out) :: i,j,k

        !transformed x,y,z coordinates
        double precision :: xp, yp, zp

        xp=x-xmin
        yp=y-ymin
        zp=z-zmin

        i=floor(xp/grid%dx)
        j=floor(yp/grid%dy)
        k=floor(zp/grid%dz)

    end subroutine


end module grid_mod
