!reads in parcel options from config file and allocates memory
!also places uniformly placed parcels in cells 
module parcelsetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use parcel_interpolation_mod, only: initialise_parcel_interp, finalise_parcel_interp, x_coords, y_coords, z_coords

  implicit none

  integer :: maxparcels_global, maxparcels_local
  integer :: nprocs
  integer :: myrank
  integer :: n_per_dir

contains

  type(component_descriptor_type) function parcelsetup_get_descriptor()
    parcelsetup_get_descriptor%name="parcelsetup"
    parcelsetup_get_descriptor%version=0.1
    parcelsetup_get_descriptor%initialisation=>initialisation_callback
    parcelsetup_get_descriptor%finalisation=>finalisation_callback
  end function parcelsetup_get_descriptor


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state


    myrank=current_state%parallel%my_rank
    nprocs=current_state%parallel%processes

    if (myrank .eq. 0) print *, "Parcel Setup - initialise"

    !get options from config file
    call read_configuration(current_state)

    if (myrank .eq. 0) print *, "maxparcels_global=",maxparcels_global

    if (mod(maxparcels_global,nprocs) .ne. 0) then
      error stop "Error: maxparcels not divisible by number of processes"
    endif

    maxparcels_local=maxparcels_global/nprocs


    if (myrank .eq. 0) print *, "maxparcels_local=",maxparcels_local

    current_state%parcels%maxparcels_global=maxparcels_global
    current_state%parcels%maxparcels_local=maxparcels_local

    if (myrank .eq. 0) print *, "allocating parcels"

    allocate(current_state%parcels%x(maxparcels_local))
    allocate(current_state%parcels%y(maxparcels_local))
    allocate(current_state%parcels%z(maxparcels_local))
    allocate(current_state%parcels%p(maxparcels_local))
    allocate(current_state%parcels%q(maxparcels_local))
    allocate(current_state%parcels%r(maxparcels_local))
    allocate(current_state%parcels%dxdt(maxparcels_local))
    allocate(current_state%parcels%dydt(maxparcels_local))
    allocate(current_state%parcels%dzdt(maxparcels_local))
    allocate(current_state%parcels%dpdt(maxparcels_local))
    allocate(current_state%parcels%dqdt(maxparcels_local))
    allocate(current_state%parcels%drdt(maxparcels_local))
    allocate(current_state%parcels%h(maxparcels_local))
    allocate(current_state%parcels%b(maxparcels_local))
    allocate(current_state%parcels%vol(maxparcels_local))
    allocate(current_state%parcels%stretch(maxparcels_local))
    allocate(current_state%parcels%tag(maxparcels_local))

    !initialise parcel interpolation 'component' of model core
    call initialise_parcel_interp(current_state)


    call setup_parcels(current_state)


    if (myrank .eq. 0) print *, "parcel setup done"


  end subroutine initialisation_callback



  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    if (myrank .eq. 0) print *, "Parcel setup - finalise"
    if (myrank .eq. 0) print *, "deallocating parcels..."

    deallocate(current_state%parcels%x)
    deallocate(current_state%parcels%y)
    deallocate(current_state%parcels%z)
    deallocate(current_state%parcels%p)
    deallocate(current_state%parcels%q)
    deallocate(current_state%parcels%r)
    deallocate(current_state%parcels%dxdt)
    deallocate(current_state%parcels%dydt)
    deallocate(current_state%parcels%dzdt)
    deallocate(current_state%parcels%dpdt)
    deallocate(current_state%parcels%dqdt)
    deallocate(current_state%parcels%drdt)
    deallocate(current_state%parcels%h)
    deallocate(current_state%parcels%b)
    deallocate(current_state%parcels%vol)
    deallocate(current_state%parcels%stretch)
    deallocate(current_state%parcels%tag)

    call finalise_parcel_interp(current_state)

    if (myrank .eq. 0) print *, "done!"


  end subroutine finalisation_callback



  subroutine setup_parcels(state)
    type(model_state_type), intent(inout) :: state
    integer :: n_per_cell
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION) :: dx, dy, dz
    integer :: nnx, nny, nnz
    integer :: nparcels
    integer :: n, i, j, k
    integer :: xstart, xstop, ystart, ystop, zstart, zstop
    real(kind=DEFAULT_PRECISION) :: ddx, ddy, ddz, xs, ys, zs
    integer :: ii, jj, kk

    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)



    !calculate number of parcels in a cell
    n_per_cell=n_per_dir**3

    !calculate number of cells in each direction (total grid size minus halo cells)
    nnx=nx - 2*state%local_grid%halo_size(3) - 1
    nny=ny - 2*state%local_grid%halo_size(2) - 1
    nnz=nz - 2*state%local_grid%halo_size(1) - 1

    nparcels=n_per_cell*(nnx)*(nny)*(nnz)

    if (state%parallel%my_rank .eq. 0) print*, n_per_dir, "parcels per direction per cell"
    if (state%parallel%my_rank .eq. 0) print*, "So", n_per_cell,"Parcels per cell"
    if (state%parallel%my_rank .eq. 0) print*, "Setting up", nparcels,"parcels in", nnx*nny*nnz, "cells"

    !print*, "nnx, nny, nnz=", nnx, nny, nnz
    !print*, "maxparcels=",state%parcels%maxparcels_local,  " nparcels=", nparcels

    if (nparcels .gt. state%parcels%maxparcels_local) then
      error stop "Maxparcels is not big enough for the number of parcels per cell requested"
    endif

    state%parcels%numparcels_local = nparcels

    !print *, "dx, dy, dz=", dx, dy, dz

    !start and end indices of the bit of grid belonging to that process
    xstart=state%local_grid%halo_size(3)+1
    ystart=state%local_grid%halo_size(2)+1
    zstart=state%local_grid%halo_size(1)+1

    xstop=nx-state%local_grid%halo_size(3)-1
    ystop=ny-state%local_grid%halo_size(2)-1
    zstop=nz-state%local_grid%halo_size(1)-1

    !print*, "xstart, xstop=", xstart, xstop
    !print*, "ystart, ystop=", ystart, ystop
    !print*, "zstart, zstop=", zstart, zstop

    !loop over every cell and put n_per_cell parcels in it
    n=0
    do i=xstart,xstop
      dx = x_coords(i+1)-x_coords(i)
      ddx=dx/n_per_dir
      do j=ystart,ystop
        dy = y_coords(j+1)-y_coords(j)
        ddy=dy/n_per_dir
        do k=zstart,zstop
          dz = z_coords(k+1)-z_coords(k)
          ddz=dz/n_per_dir

          do ii=1,n_per_dir
            do jj=1,n_per_dir
              do kk=1,n_per_dir
                n=n+1
                xs=x_coords(i)+ddx/2+(ii-1)*ddx
                ys=y_coords(j)+ddy/2+(jj-1)*ddy
                zs=z_coords(k)+ddz/2+(kk-1)*ddz

                state%parcels%x(n) = xs
                state%parcels%y(n) = ys
                state%parcels%z(n) = zs
                state%parcels%vol(n) = 1.0
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

       ! do n=1,8
       !   print *, state%parcels%x(n), state%parcels%y(n), state%parcels%z(n)
       ! enddo

      !print*, maxval(state%parcels%x(1:n)), xmax
      !print*, maxval(state%parcels%y(1:n)), ymax
      !print*, maxval(state%parcels%z(1:n)), zmax


      ! print *, n, nparcels

      print*, "parcels initialised"


  end subroutine




  subroutine read_configuration(state)
    type(model_state_type), intent(inout) :: state

    maxparcels_global=options_get_integer(state%options_database,"max_parcels")
    n_per_dir=options_get_integer(state%options_database,"parcels_per_cell_dir")

  end subroutine read_configuration






end module
