module testInterp_mod
    use global_mod
    use interpolation_mod
    use grid_mod
    use parcel_mod

    implicit none

    double precision, parameter :: tol=1.d-9

contains

    !tests whether a grid has been interpolated onto correctly
    subroutine testgrid(grid, variable)
        type(gridded_variable), intent(in) :: grid
        integer, intent(in) :: variable !which variable we interpolated from

        integer :: i, j, k
        integer :: nx, ny, nz
        double precision :: val

        double precision :: t2, t1
        double precision :: ref

        !if variable = vort_x  then should have 1
        !if variable = vort_y then should be -1
        !if vairable = vort_z then should be 0

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz

        t1=MPI_Wtime()

        !check internal values only (will be incorrect on boundaries)
        do k=2,nz-1
            if(variable .eq. VORT_Z) ref = zmin + (k-1)*grid%dz
            do j=2,ny-1
                if (variable .eq. VORT_Y) ref = ymin + (j-1)*grid%dy
                do i=2,nx-1
                    if (variable .eq. vort_x) ref = xmin + (i-1)*grid%dx


                    val=grid%data(i,j,k)

                    if (variable .eq. VORT_X) then
                        if (abs(val-ref) .gt. tol) then
                            print*, "testgrid error vort_x"
                            print*, i,j,k
                            print*, "Error:", val, ref, abs(val-ref)
                            stop
                        endif
                    else if (variable .eq. VORT_Y) then
                        if (abs(val-ref) .gt. tol) then
                            print*, "testgrid error vort_y"
                            print*, i,j,k
                            print*, "Error:", val, ref, abs(val-ref)
                            stop
                        endif
                    else if (variable .eq. VORT_Z) then
                        if (abs(val-ref) .gt. tol) then
                            print*, "testgrid error vort_z"
                            print *, i,j,k
                            print*, "Error:", val, ref, abs(val-ref)
                            stop
                        endif
                    endif
                enddo
            enddo
        enddo

        t2=MPI_Wtime()

        if (debug) print *, "verified correctly in t=", t2-t1

    end subroutine

    !tests whether we have interpolated onto the parcel correctly
    subroutine testparcel(structure,variable,gridtype)
        logical, intent(in) :: structure
        integer, intent(in) :: variable
        integer, intent(in) :: gridtype !which type of grid was set up

        !if gridtype = 0 then vaue should be equal to 2.0
        !if gridtype = 1 then value should be equal to x coord of parcel
        !if gridtype = 2 then value should be equal to y coord of parcel
        !if gridtype = 3 then value should be equal to z coord of parcel

        integer :: n
        double precision :: val, ref

        double precision :: t2, t1


        t1=MPI_Wtime()
        do n=1,nparcels


            !retrieve value

            if (variable .eq. VEL_X) then
                if (structure) then
                    val = parcels(n)%u
                else
                    val = uvel(n)
                endif
            else if (variable .eq. VEL_Y) then
                if (structure) then
                    val = parcels(n)%v
                else
                    val = vvel(n)
                endif
            else if (variable .eq. VEL_Z) then
                if (structure) then
                    val = parcels(n)%w
                else
                    val = wvel(n)
                endif
            endif

            !retrieve 'reference' value

            if (gridtype .eq. 0) then
                ref=2.d0
            else if (gridtype .eq. 1) then
                if (structure) then
                    ref = parcels(n)%x
                else
                    ref=xpos(n)
                endif
            else if (gridtype .eq. 2) then
                if (structure) then
                    ref = parcels(n)%y
                else
                    ref=ypos(n)
                endif
            else if (gridtype .eq. 3) then
                if (structure) then
                    ref = parcels(n)%z
                else
                    ref=zpos(n)
                endif
            endif

            !compare refrence with value

            if (abs(val-ref) .gt. tol) then
                print*, "testparcel error"
                print*, "Error:", val, ref, abs(val-ref)
                print *, "Parcel number=",n
                if (structure) then
                    print*, "i,j,k=",parcels(n)%i, parcels(n)%k, parcels(n)%k
                else
                    print*, "i,j,k=",is(n), js(n), ks(n)
                endif
                stop
            endif

        enddo
        t2=MPI_Wtime()

        if (debug) print *, "verified correctly in t=", t2-t1

    end subroutine

end module testInterp_mod
