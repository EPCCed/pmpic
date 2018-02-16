!module containing subroutines to initialise grids and parcels

!this is mosty just a placeholder to create somewhat plausible test data - not for production runs

module initialise_mod
    use global_mod
    use grid_mod
    use parcel_mod

    implicit none

contains

!allocates memory for parcels and sets up their initial positions
    subroutine initialise_parcels(structure, shuffle)
        implicit none
        logical, intent(in) :: structure
        logical, intent(in) :: shuffle

        double precision :: dx, dy, dz
        double precision :: ipos, jpos, kpos
        integer :: i, j, k, num
        double precision :: t2, t1

        t1=MPI_Wtime()
        if (structure) then
            allocate(parcels(maxparcels))
        else
            allocate(xpos(maxparcels))
            allocate(ypos(maxparcels))
            allocate(zpos(maxparcels))
            allocate(uvel(maxparcels))
            allocate(vvel(maxparcels))
            allocate(wvel(maxparcels))
            allocate(bl(maxparcels))
            allocate(thetal(maxparcels))
            allocate(q(maxparcels))
            allocate(volume(maxparcels))
            allocate(active(maxparcels))
            allocate(dummyvars(10,maxparcels))
            allocate(rvort(maxparcels))
            allocate(svort(maxparcels))
            allocate(tvort(maxparcels))
        endif
        t2=MPI_Wtime()

        Print *, "Allocation time=",t2-t1

        !not a necessary step, but marks all parcels as inactive (empty)
        ! if (structure) then
        !     t1=MPI_Wtime()
        !     do i=1,maxparcels
        !         parcels(i)%active=.false.
        !     enddo
        !     t2=MPI_Wtime()
        !     print *, "Parcel structure initialisation time=",t2-t1
        ! endif

        !assume we have 2 parcels in each direction in a cell, so 8 per cell
        nparcels=8*nxbase*nybase*nzbase

        print *, "There will be ",nparcels," parcels"

        !Now place initial parcels:

        !spacing of the parcels
        dx=(xmax-xmin)/nxbase/2.
        dy=(ymax-ymin)/nybase/2.
        dz=(zmax-zmin)/nzbase/2.



        t1=MPI_Wtime()

!$OMP PARALLEL DO PRIVATE(i,j,k,kpos,jpos,ipos,num)
        do k=1,nzbase*2
            kpos=0.25*dz + (k-1)*dz
            do j=1,nybase*2
                jpos=0.25*dy + (j-1)*dy

                do i=1,nxbase*2
                    !calculate parcel number
                    num=(nxbase*2*nybase*2)*(k-1) + (j-1)*nxbase*2 + i

                    ipos=0.25*dx + (i-1)*dx

                    if (structure) then

                        !set the parcel's position
                        parcels(num)%x = ipos
                        parcels(num)%y = jpos
                        parcels(num)%z = kpos
                        parcels(num)%volume = 1.d0

                        parcels(num)%r = 1.d0
                        parcels(num)%s = -1.d0
                        parcels(num)%t = 0.d0

                        !set that parcel to be active
                        parcels(num)%active=.true.

                    else
                        xpos(num) = ipos
                        ypos(num) = jpos
                        zpos(num) = kpos
                        volume(num) = 1.d0
                        rvort(num)= 1.d0
                        svort(num)= -1.d0
                        tvort(num)= 0.d0
                        active(num) = .true.

                    endif

                enddo
            enddo
        enddo
!$OMP END PARALLEL DO

        t2=MPI_Wtime()

        print *, "Setup time", t2-t1

        if (shuffle) then
            print *, "Shuffling parcels:"
            t1=MPI_Wtime()
            call shuffle_parcels(nparcels,structure)
            t2=MPI_Wtime()
            print *, "Shuffling time=",t2-t1
        endif


        !print *, " Set up", num, "parcels"

    end subroutine




!deallcoates parcels
    subroutine finalise_parcels(structure)
        implicit none
        logical, intent(in) :: structure

        if (structure) then
            deallocate(parcels)
        else
            deallocate(xpos, ypos, zpos)
            deallocate(uvel, vvel, wvel)
            deallocate(bl,thetal,q,volume)
            deallocate(active)
            deallocate(dummyvars)
        endif
        print *, "Deallocating parcels"
    end subroutine


!sets up a grid variable
    subroutine set_grid(grid,type)
        implicit none
        type(gridded_variable), intent(inout) :: grid
        integer, intent(in) :: type !flag to set the contents of the grid:
        ! 0 - uniform value
        ! 1 - gradient in x
        ! 2 - gradient in y
        ! 3 - gradient in z
        integer :: i,j,k

        print *, "Setting up grid with type=",type


        if (type .eq. 0) then

            grid%data(:,:,:) = 2.d0

        else if (type .eq. 1) then

            do i=1,grid%nx
                grid%data(i,:,:) = i*grid%dx
            enddo

        else if (type .eq. 2) then

            do j=1,grid%ny
                grid%data(:,j,:) = j*grid%dy
            enddo

        else if (type .eq. 3) then

            do k=1,grid%nz
                grid%data(:,:,k) = k*grid%dz
            enddo

        else

            STOP "INVALID type value"

        endif



    end subroutine





    subroutine shuffle_parcels(n,structure)
  ! performs a Fisher-Yates shuffle (aka Knuth shuffle) of the indices
  ! This may be either beneficial or bad for shared-balance parallelism in some
  ! instances, as it can help prevent processors accessing the same grid variables
  ! multiple times in sequence
        integer, intent(in) :: n
        logical, intent(in) :: structure
        integer :: shuffle_index, rand_target
        double precision :: temp_var
        type(parcel) :: temp_par
        real :: random_out


        do shuffle_index = n, 2, -1
            call random_number(random_out)
            rand_target = int(random_out * shuffle_index) + 1

            if (structure) then
                temp_par = parcels(rand_target)
                parcels(rand_target)=parcels(shuffle_index)
                parcels(shuffle_index) = temp_par
            else
                temp_var=xpos(rand_target)
                xpos(rand_target)=xpos(shuffle_index)
                xpos(shuffle_index) = temp_var

                temp_var=ypos(rand_target)
                ypos(rand_target)=ypos(shuffle_index)
                ypos(shuffle_index) = temp_var

                temp_var=zpos(rand_target)
                zpos(rand_target)=zpos(shuffle_index)
                zpos(shuffle_index) = temp_var
            endif

        end do

    end subroutine


end module
