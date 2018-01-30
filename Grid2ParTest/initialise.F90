!module containing subroutines to initialise grids and parcels

!this is mosty just a placeholder to create somewhat plausible test data - not for production runs

module initialise_mod
    use global_mod
    use grid_mod
    use parcel_mod

    implicit none

contains

    subroutine initialise_parcels(structure)
        implicit none
        logical, intent(in) :: structure

        if (structure) then
            allocate(parcels(maxparcels))
        else
            allocate(xpos(maxparcels))
            allocate(ypos(maxparcels))
            allocate(zpos(maxparcels))
            allocate(uvel(maxparcels))
            allocate(vvel(maxparcels))
            allocate(wvel(maxparcels))
        endif

        !assume we have 2 parcels in each direction in a cell, so 8 per cell
        nparcels=8*nxbase*nybase*nzbase

        print *, "There will be ",nparcels," parcels"

    end subroutine

    subroutine finalize_parcels(structure)
        implicit none
        logical, intent(in) :: structure

        if (structure) then
            deallocate(parcels)
        else
            deallocate(xpos, ypos, zpos)
            deallocate(uvel, vvel, wvel)
        endif
        print *, "Deallocating parcels"
    end subroutine


    subroutine set_grid(grid)
        implicit none
        type(gridded_variable), intent(inout) :: grid

        grid%data(:,:,:) = 2.d0

    end subroutine

end module
