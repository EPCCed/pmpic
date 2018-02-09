Module interpolation_mod
    use global_mod
    use grid_mod
    use parcel_mod

    implicit none

    INTEGER, PARAMETER :: VEL_X = 0
    INTEGER, PARAMETER :: VEL_Y = 1
    INTEGER, PARAMETER :: VEL_Z = 2

    INTEGER, PARAMETER :: VORT_X = 10
    INTEGER, PARAMETER :: VORT_Y = 11
    INTEGER, PARAMETER :: VORT_Z = 12

contains

    subroutine grid2par(grid, structure, variable)
        type(gridded_variable), intent(in) :: grid
        logical, intent(in) :: structure !flag to tell the subroutine if it's using arrays or structures
        integer, intent(in) :: variable ! flag to tell the subroutine which variable to interpolate into

        integer :: n
        integer :: i, j, k

        double precision :: delx, dely, delz
        double precision :: c000, c001, c010, c011, c100, c101, c110, c111
        double precision :: c00, c01, c10, c11
        double precision :: c0, c1
        double precision :: c

        double precision :: t2, t1

        t1=MPI_Wtime()

        !loop over each parcel
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP             FIRSTPRIVATE(structure,variable) &
!$OMP             SHARED(grid,parcels,xpos, ypos, zpos, uvel, vvel, wvel) &
!$OMP             schedule(guided)
        do n=1, nparcels

            !get lower left corner of the cell the particle is contained within:

            if (structure) then
                call Get_Grid_Coords(grid,parcels(n)%x,parcels(n)%y,parcels(n)%z,i,j,k, delx, dely, delz)
            else
                call Get_Grid_Coords(grid,xpos(n),ypos(n),zpos(n),i,j,k,delx, dely, delz)
            endif

            !get corners of cube around parcel

            c000 = grid%data(i,j,k)
            c001 = grid%data(i,j,k+1)
            c010 = grid%data(i,j+1,k)
            c011 = grid%data(i,j+1,k+1)
            c100 = grid%data(i+1,j,k)
            c101 = grid%data(i+1,j,k+1)
            c110 = grid%data(i+1,j+1,k)
            c111 = grid%data(i+1,j+1,k+1)

            !interpolate in x direction to produce square around the parcel

            c00 = c000*(1-delx) + c100*delx
            c01 = c001*(1-delx) + c101*delx
            c10 = c010*(1-delx) + c110*delx
            c11 = c011*(1-delx) + c111*delx

            !interpolate in y direction to produce line through parcel

            c0 = c00*(1-dely) + c10*dely
            c1 = c01*(1-dely) + c11*dely

            !interpolate to parcel position

            c = c0*(1-delz) + c1*delz

            !now update parcel's variables

            if (variable .eq. VEL_X) then
                if (structure) then
                    parcels(n)%u = c
                else
                    uvel(n) = c
                endif
            else if (variable .eq. VEL_Y) then
                if (structure) then
                    parcels(n)%v = c
                else
                    vvel(n) = c
                endif
            else if (variable .eq. VEL_Z) then
                if (structure) then
                    parcels(n)%w = c
                else
                    wvel(n) = c
                endif
            else
                STOP "Undefined variable"
            ENDIF

        enddo
!$OMP END PARALLEL DO

        t2=MPI_Wtime()

        print*, "grid2par time=", t2-t1

    end subroutine



    subroutine par2grid(grid,structure,variable)
        type(gridded_variable), intent(inout) :: grid
        logical, intent(in) :: structure
        integer, intent(in) :: variable

        integer :: n
        double precision, allocatable, dimension(:,:,:) :: weights, data
        double precision :: w
        double precision :: v
        double precision :: delx, dely, delz
        double precision :: w000, w001, w010, w011, w100, w101, w110, w111
        integer :: i,j,k

        double precision :: t2, t1

        allocate(weights(grid%nx, grid%ny, grid%nz))
        allocate(data(grid%nx,grid%ny,grid%nz))

        !zero the weights and the grid
!$OMP PARALLEL DO
        do n=1,grid%nz
            weights(:,:,n) = 0.0d0
            data(:,:,n) = 0.0d0
        enddo
!$OMP END PARALLEL DO

        t1=MPI_Wtime()

        !loop over each parcel and add its contribtion to the grid

!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP             shared(parcels, structure, volume, grid, rvort, svort,tvort, variable)&
!$OMP             SCHEDULE(GUIDED)&
!$OMP             reduction(+:weights,data)
        do n=1,nparcels

            !get lower left corner of the cell the particle is contained within:

            if (structure) then
                call Get_Grid_Coords(grid,parcels(n)%x,parcels(n)%y,parcels(n)%z,i,j,k, delx, dely, delz)
            else
                call Get_Grid_Coords(grid,xpos(n),ypos(n),zpos(n),i,j,k,delx, dely, delz)
            endif


            !retrieve value of parcel variable and volume

            if (VARIABLE .eq. VORT_X) then
                if (structure) then
                    w = parcels(n)%r
                else
                    w = rvort(n)
                endif
            else if (VARIABLE .eq. VORT_Y) then
                if (structure) then
                    w = parcels(n)%s
                else
                    w = svort(n)
                endif
            else if (VARIABLE .eq. VORT_Z) then
                if (structure) then
                    w = parcels(n)%t
                else
                    w = tvort(n)
                endif
            else
                STOP "Invalid variable"
            endif

            if (structure) then
                v = parcels(n)%volume
            else
                v = volume(n)
            endif

            !calculate weights on each vertex of cube and add that to grid subtotals

            w000 = (1-delx)*(1-dely)*(1-delz)*v
!!$OMP ATOMIC
            data(i,j,k) = data(i,j,k) + w000*w
!!$OMP ATOMIC
            weights(i,j,k) = weights(i,j,k) + w000



            w001 = (1-delx)*(1-dely)*(delz)*v
!!$OMP ATOMIC
            data(i,j,k+1) = data(i,j,k+1) + w001*w
!!$OMP ATOMIC
            weights(i,j,k+1) = weights(i,j,k+1) + w001



            w010 = (1-delx)*(dely)*(1-delz)*v
!!$OMP ATOMIC
            data(i,j+1,k) = data(i,j+1,k) + w010*w
!!$OMP ATOMIC
            weights(i,j+1,k) = weights(i,j+1,k) + w010



            w011 = (1-delx)*(dely)*(delz)*v
!!$OMP ATOMIC
            data(i,j+1,k+1) = data(i,j+1,k+1) + w011*w
!!$OMP ATOMIC
            weights(i,j+1,k+1) = weights(i,j+1,k+1) + w011



            w100 = (delx)*(1-dely)*(1-delz)*v
!!$OMP ATOMIC
            data(i+1,j,k) = data(i+1,j,k) + w100*w
!!$OMP ATOMIC
            weights(i+1,j,k) = weights(i+1,j,k) + w100



            w101 = (delx)*(1-dely)*(delz)*v
!!$OMP ATOMIC
            data(i+1,j,k+1) = data(i+1,j,k+1) + w101*w
!!$OMP ATOMIC
            weights(i+1,j,k+1) = weights(i+1,j,k+1) + w101



            w110 = (delx)*(dely)*(1-delz)*v
!!$OMP ATOMIC
            data(i+1,j+1,k) = data(i+1,j+1,k) + w110*w
!!$OMP ATOMIC
            weights(i+1,j+1,k) = weights(i+1,j+1,k) + w110



            w111 = (delx)*(dely)*(delz)*v
!!$OMP ATOMIC
            data(i+1,j+1,k+1) = data(i+1,j+1,k+1) + w111*w
!!$OMP ATOMIC
            weights(i+1,j+1,k+1) = weights(i+1,j+1,k+1) + w111


        enddo
!$OMP END PARALLEL DO

        !divide grid by weights to get the value of the gridded variable
!$OMP PARALLEL DO
        do n=1,grid%nz
            grid%data(:,:,n) = data(:,:,n)/weights(:,:,n)
        enddo
!$OMP END PARALLEL DO

        t2=MPI_Wtime()

        print *, "par2grid time=", t2-t1

        deallocate(weights)

    end subroutine

















end module
