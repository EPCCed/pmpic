module grid2partest_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use parcel_interpolation_mod, only: cache_parcel_interp_weights, grid2par, x_coords, y_coords, z_coords
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array

  implicit none

  integer :: nx, ny, nz
  real(kind=DEFAULT_PRECISION) :: dx, dy, dz, xmin, xmax, ymin, ymax, zmin,zmax

contains

  type(component_descriptor_type) function grid2partest_get_descriptor()
    grid2partest_get_descriptor%name="grid2partest"
    grid2partest_get_descriptor%version=0.1
    grid2partest_get_descriptor%initialisation=>initialisation_callback
    grid2partest_get_descriptor%timestep=>timestep_callback
    grid2partest_get_descriptor%finalisation=>finalisation_callback
  end function grid2partest_get_descriptor


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    print*, "grid2partest initialisation"

    !set up grid

    call setup_grid(current_state)

    !set up parcels

    call setup_parcels(current_state)

  end subroutine



  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    print*, "grid2partest timestep"

    call cache_parcel_interp_weights(current_state)

    call grid2par(current_state,current_state%u,current_state%parcels%u)
    call check_result(current_state,current_state%parcels%u,current_state%parcels%x)

    call grid2par(current_state,current_state%v,current_state%parcels%v)
    call check_result(current_state,current_state%parcels%v,current_state%parcels%y)

    call grid2par(current_state,current_state%w,current_state%parcels%w)
    call check_result(current_state,current_state%parcels%w,current_state%parcels%z)

  end subroutine

  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    print*, "grid2partest finalisation"

  end subroutine


  subroutine setup_grid(state)
    type(model_state_type), intent(inout) :: state
    integer :: i, j, k

    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)

    !print *, "setup_grid", nx, ny, nz
    !print*, "y_coords=", y_coords(1), y_coords(2),y_coords(ny-1),y_coords(ny)



    !get nx, ny, nz

    allocate(state%u%data(nz,ny,nx))
    allocate(state%v%data(nz,ny,nx))
    allocate(state%w%data(nz,ny,nx))

    do i=1,nx
      do j=1,ny
        do k=1,nz
          state%u%data(k,j,i) = x_coords(i)
          state%v%data(k,j,i) = y_coords(j)
          state%w%data(k,j,i) = z_coords(k)
        enddo
      enddo
    enddo

    !print *, x_coords(1), x_coords(2)
    !print *, state%u%data(1,1,1), state%u%data(1,1,2)

    xmin = x_coords(1+state%local_grid%halo_size(3))
    ymin = y_coords(1+state%local_grid%halo_size(2))
    zmin = z_coords(1+state%local_grid%halo_size(1))

    xmax = x_coords(nx-state%local_grid%halo_size(3))
    ymax = y_coords(ny-state%local_grid%halo_size(2))
    zmax = z_coords(nz-state%local_grid%halo_size(1))

    !print*, "xmin, xmax=", xmin, xmax
    !print*, "ymin, ymax=", ymin, ymax
    !print*, "zmin, zmax=", zmin, zmax



    print *, "initialised grids"

  end subroutine

  subroutine setup_parcels(state)
    type(model_state_type), intent(inout) :: state
    integer :: n_per_cell, n_per_dir
    integer :: nnx, nny, nnz
    integer :: nparcels
    integer :: n, i, j, k
    integer :: xstart, xstop, ystart, ystop, zstart, zstop
    real(kind=DEFAULT_PRECISION) :: ddx, ddy, ddz, xs, ys, zs
    integer :: ii, jj, kk

    !read in number of parcels per cell in a given direction

    n_per_dir=options_get_integer(state%options_database,"parcels_per_cell")

    !calculate number of parcels in a cell
    n_per_cell=n_per_dir**3

    !calculate number of cells in each direction (total grid size minus halo cells)
    nnx=nx - 2*state%local_grid%halo_size(3) - 1
    nny=ny - 2*state%local_grid%halo_size(2) - 1
    nnz=nz - 2*state%local_grid%halo_size(1) - 1

    nparcels=n_per_cell*(nnx)*(nny)*(nnz)

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
      do j=ystart,ystop
        dy = y_coords(j+1)-y_coords(j)
        do k=zstart,zstop
          dz = z_coords(k+1)-z_coords(k)

          ddx=dx/n_per_dir
          ddy=dy/n_per_dir
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

  subroutine check_result(state, values, reference)
    type(model_state_type), intent(in) :: state
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: values, reference
    real(kind=DEFAULT_PRECISION), parameter :: tol=1.e-9
    real(kind=DEFAULT_PRECISION) :: diff
    integer :: n, nparcels

    nparcels=state%parcels%numparcels_local

    do n=1,nparcels
      diff =abs(values(n)-reference(n))
      if (diff .gt. tol) then
        print *, n, values(n), reference(n), diff
        error stop "Wrong answer"
      endif
    enddo

    print*, "values verified"

  end subroutine



end module
