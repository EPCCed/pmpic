!contains grid2par and par2grid and helper subroutines
module parcel_interpolation_mod
  use state_mod, only: model_state_type
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use MPI, only: MPI_Barrier
  use omp_lib


  implicit none

  !flag to see if this is intiialised or not
  logical :: initialised = .false.

  !cached variables
  integer, allocatable, dimension(:) ::  is, js, ks
  real(kind=DEFAULT_PRECISION), allocatable, dimension(:) :: delxs, delys, delzs
  integer :: nparcels

  !cached grid variables
  integer :: nx, ny, nz, ndz
  real(kind=DEFAULT_PRECISION) ::xmin, xmax, ymin, ymax, zmin, zmax
  real(kind=DEFAULT_PRECISION) :: dx, dy, meandz
  real(kind=DEFAULT_PRECISION), allocatable, dimension(:) :: z, dz

  real(kind=default_precision), allocatable, dimension(:) :: x_coords, y_coords, z_coords


contains

  !initialise parcel interpolation data structures (this subroutine should be called once during parcel setup)
  !This contains variables for caching interpolation values as well as variables converting grid
  !indices to coordinates
  subroutine initialise_parcel_interp(state)
    type(model_state_type), intent(in) :: state
    integer :: n
    integer :: xstart, xstop, ystart, ystop, zstart, zstop
    real(kind=DEFAULT_PRECISION) :: dzdummy

    if (initialised) error stop "parcel interpolation routines are already initialised - cannot initialise"

    !allocate cache arrays
    n=state%parcels%maxparcels_local

    allocate(is(n))
    allocate(js(n))
    allocate(ks(n))
    allocate(delxs(n))
    allocate(delys(n))
    allocate(delzs(n))

    initialised = .true.

    !get dx, dy and dz
    dx=state%global_grid%resolution(3)
    dy=state%global_grid%resolution(2)
    dzdummy=state%global_grid%resolution(1)

    !set number of cells in local grids
    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)


    if (dx .lt. 1.E-3 .or. dy .lt. 1.e-3) then
      error stop "no grid resolutions defined"
    endif


    if (state%parallel%my_rank .eq. 0 ) print *, "dx=", dx, " dy=", dy, " dz=", dzdummy

    !get global indices of the first and last elements in the local arrays
    xstart = state%local_grid%start(3)-state%local_grid%halo_size(3)
    xstop = state%local_grid%end(3)+state%local_grid%halo_size(3)

    ystart = state%local_grid%start(2)-state%local_grid%halo_size(2)
    ystop = state%local_grid%end(2)+state%local_grid%halo_size(2)

    zstart = state%local_grid%start(1)-state%local_grid%halo_size(1)
    zstop = state%local_grid%end(1)+state%local_grid%halo_size(1)


    !try to set up z
    if (allocated(state%global_grid%configuration%vertical%z)) then
      !if it is already set up then take its values from the state
    !  nz=size(state%global_grid%configuration%vertical%z)
      print*, "nz=",nz
      allocate(z(nz))
      z=state%global_grid%configuration%vertical%z


      ndz=size(state%global_grid%configuration%vertical%dz)
      print*, "dzn=",ndz
      allocate(dz(ndz))
    else !else we set it up ourselves
      if (state%parallel%my_rank .eq. 0 ) print *, "Warning: no z grid defined. Creating uniform grid"

      !nz=zstop-zstart+1
      allocate(z(nz))
      ndz=nz-1
      allocate(dz(ndz))

      !set dz to dx since we have no other estimate
      dz(1:ndz)=dzdummy

      meandz=sum(dz)/(nz-1)


      z(1)=(zstart-1)*dz(1)
      do n=2,nz
        z(n) = z(n-1)+dz(n-1)
      enddo
    endif

    !print*, state%parallel%my_rank, xstart, xstop, ystart, ystop, zstart, zstop
    ! do n=1,nz
    !   print*, state%parallel%my_rank, n, z(n)
    ! enddo

    xmin = (xstart-1)*dx !Coordinate of first point in the x grid
    ymin = (ystart-1)*dy
    zmin = z(1)

    xmax = (xstop-1)*dx !coordinate of last point in the x grid
    ymax = (ystop-1)*dy
    zmax = z(nz)

    !allocate arrays that will hold the coordinates of x,y and z for each grid cell
    allocate(x_coords(nx), y_coords(ny))
    allocate(z_coords(nz))

    do n=1,nx
      x_coords(n) = xmin + (n-1)*dx
    enddo

    do n=1,ny
      y_coords(n) = ymin + (n-1)*dy
    enddo


    z_coords(:) = z(:)

    !print *, xmin, xmax, ymin, ymax, zmin, zmax

    call flush()

    ! if (state%parallel%my_rank .eq. 0) then
    ! do n=1,nz
    !    print*, state%parallel%my_rank, n, z(n)
    ! enddo
    ! endif

    !sanity check to see if grid is set up right
    !print *, state%parallel%my_rank, ymin+dy, ymax-2*dy

    call MPI_Barrier(state%parallel%monc_communicator,n)

    if (state%parallel%my_rank .eq. 0 ) print *, "parcel_interp initialised"

  end subroutine


 !finalise parcel interpolation data structures
  subroutine finalise_parcel_interp(state)
    type(model_state_type), intent(in) :: state

    if (.not. initialised) error stop "parcel interpolation routines are not initialised - cannot finalise"

    deallocate(is)
    deallocate(js)
    deallocate(ks)
    deallocate(delxs)
    deallocate(delys)
    deallocate(delzs)

    deallocate(z)
    deallocate(dz)

    deallocate(x_coords,y_coords,z_coords)

    initialised=.false.

    if (state%parallel%my_rank .eq. 0 ) print *, "parcel_interp finalised"

  end subroutine


  !cache weights for parcels (should be called once after parcel positions have changed before interpolating values)
  subroutine cache_parcel_interp_weights(state)
    type(model_state_type), intent(in) :: state
    integer :: n
    integer :: i, j, k
    integer :: nn
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz
    real(kind=DEFAULT_PRECISION) :: xp, yp, zp


    nparcels=state%parcels%numparcels_local

    !for each parcel we determine which cell it belongs to

    !maybe have one loop for x, one for y and one for z for vectorisation?

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,xp,yp,zp,i,j,k,delx,dely,delz)
    !$OMP DO
    do n=1,nparcels

      !get positions of parcels
      xp=state%parcels%x(n)
      yp=state%parcels%y(n)
      zp=state%parcels%z(n)

      !get the index of the lower left corner of the cell that the parcel is in
      i=floor(xp-xmin)/dx+1
      j=floor(yp-ymin)/dy+1
      do nn=1,nz-1 !as z may be a variable size grid we need to search through z to get the cell
        if ((zp .gt. z(nn)) .and. (zp .lt. z(nn+1))) then
          k=nn
          exit
        endif
      enddo

      if ((xp .lt. xmin) .or. (xp .gt. xmax)) error stop "x too big/small"
      if ((yp .lt. ymin) .or. (yp .gt. ymax)) error stop "y too big/small"
      if ((zp .lt. zmin) .or. (zp .gt. zmax)) error stop "z too big/small"

      !get the fractional position in the cell (from the lower left corner) that the parcel is at
      delx= ((xp-xmin)-(i-1)*dx)/dx
      dely= ((yp-ymin)-(j-1)*dy)/dy
      delz= ((zp-zmin)-(k-1)*dz(k))/dz(k)


      if (delx .gt. 1. .or. delx .lt. 0) error stop "delx wrong size"
      if (dely .gt. 1. .or. dely .lt. 0) error stop "dely wrong size"
      if (delz .gt. 1. .or. delz .lt. 0) error stop "delz wrong size"

      !cache these values
      is(n) = i
      js(n) = j
      ks(n) = k

      delxs(n)=delx
      delys(n)=dely
      delzs(n)=delz


    enddo
    !$OMP END DO

    !$OMP MASTER
    print *, "nthreads=", omp_get_num_threads()
    !$OMP END MASTER

    !$OMP END PARALLEL

    print *, "weights cached"

  end subroutine



  !interpolate gridded variable (grid) to parcel variable (var)
  subroutine grid2par(state,grid,var)
    type(model_state_type), intent(inout) :: state
    type(prognostic_field_type), intent(in) :: grid
    real(kind=DEFAULT_PRECISION), dimension(:) :: var

    integer :: n
    real(kind=DEFAULT_PRECISION) :: c000, c001, c010, c011, c100, c101, c110, c111
    real(kind=DEFAULT_PRECISION) :: c00, c01, c10, c11
    real(kind=DEFAULT_PRECISION) :: c0, c1
    real(kind=DEFAULT_PRECISION) :: c

    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nparcels,delxs,delys,delzs,is,js,ks,grid,var)
    do n=1,nparcels

      !retrieve cached values
      delx=delxs(n)
      dely=delys(n)
      delz=delzs(n)
      i=is(n)
      j=js(n)
      k=ks(n)

      !retrieve corners of grid cell the nth parcel is in

      c000 = grid%data(k,j,i)
      c001 = grid%data(k,j,i+1)
      c010 = grid%data(k,j+1,i)
      c011 = grid%data(k,j+1,i+1)
      c100 = grid%data(k+1,j,i)
      c101 = grid%data(k+1,j,i+1)
      c110 = grid%data(k+1,j+1,i)
      c111 = grid%data(k+1,j+1,i+1)

      !interpolate in z direction to produce square around parcel in y-x plane
      c00 = c000*(1-delz) + c100*delz
      c01 = c001*(1-delz) + c101*delz
      c10 = c010*(1-delz) + c110*delz
      c11 = c011*(1-delz) + c111*delz

      !interpolate in y direction to produce line through parcel along x direction

      c0 = c00*(1-dely) + c10*dely
      c1 = c01*(1-dely) + c11*dely

      !interpolate to parcel position

      c = c0*(1-delx) + c1*delx

      !now update parcel's variable

      var(n) = c

    enddo
    !$OMP END PARALLEL DO

  end subroutine



  subroutine par2grid(state,var,grid)

    implicit none

    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: var
    type(prognostic_field_type), intent(inout) :: grid

    integer :: n
    double precision, allocatable, dimension(:,:,:) :: weights, data
    double precision :: w
    double precision :: v
    double precision :: delx, dely, delz
    double precision :: w000, w001, w010, w011, w100, w101, w110, w111
    integer :: i,j,k
    integer :: xhalo, yhalo, zhalo



    allocate(weights(nz,ny,nx))
    allocate(data(nz,ny,nx))

    !zero the weights and the grid
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP              SHARED(nparcels,nx,ny,nz,delx,dely,delz,state,var,grid, weights,data)

 !$OMP DO
    do n=1,nx
        weights(:,:,n) = 0.0d0
        data(:,:,n) = 0.0d0
    enddo
!$OMP END DO

!$OMP DO REDUCTION(+:weights,data)
    do n=1,nparcels

      !get cached grid positions
      delx=delxs(n)
      dely=delys(n)
      delz=delzs(n)
      i=is(n)
      j=js(n)
      k=ks(n)

      !parcel's volume
      v = state%parcels%vol(n)

      w=var(n)

      !calculate weights on each vertex of cube and add that to grid subtotals

      ! if (n .lt. 10) then
      !   print *,n, i, j, k, v, w
      ! endif

      w000 = (1-delz)*(1-dely)*(1-delx)*v
      data(k,j,i) = data(k,j,i) + w000*w
      weights(k,j,i) = weights(k,j,i) + w000



      w001 = (1-delz)*(1-dely)*(delx)*v
      data(k,j,i+1) = data(k,j,i+1) + w001*w
      weights(k,j,i+1) = weights(k,j,i+1) + w001



      w010 = (1-delz)*(dely)*(1-delx)*v
      data(k,j+1,i) = data(k,j+1,i) + w010*w
      weights(k,j+1,i) = weights(k,j+1,i) + w010



      w011 = (1-delz)*(dely)*(delx)*v
      data(k,j+1,i+1) = data(k,j+1,i+1) + w011*w
      weights(k,j+1,i+1) = weights(k,j+1,i+1) + w011



      w100 = (delz)*(1-dely)*(1-delx)*v
      data(k+1,j,i) = data(k+1,j,i) + w100*w
      weights(k+1,j,i) = weights(k+1,j,i) + w100



      w101 = (delz)*(1-dely)*(delx)*v
      data(k+1,j,i+1) = data(k+1,j,i+1) + w101*w
      weights(k+1,j,i+1) = weights(k+1,j,i+1) + w101



      w110 = (delz)*(dely)*(1-delx)*v
      data(k+1,j+1,i) = data(k+1,j+1,i) + w110*w
      weights(k+1,j+1,i) = weights(k+1,j+1,i) + w110



      w111 = (delz)*(dely)*(delx)*v
      data(k+1,j+1,i+1) = data(k+1,j+1,i+1) + w111*w
      weights(k+1,j+1,i+1) = weights(k+1,j+1,i+1) + w111


    enddo
    !$OMP END DO

    !divide grid by weights to get the value of the gridded variable

    xhalo = state%local_grid%halo_size(3)
    yhalo = state%local_grid%halo_size(2)
    zhalo = state%local_grid%halo_size(1)

!$OMP DO
    do n=1+xhalo,nx-xhalo
        grid%data(1+zhalo:nz-zhalo,1+yhalo:ny-yhalo,n) = &
        data(1+zhalo:nz-zhalo,1+yhalo:ny-yhalo,n)/weights(1+zhalo:nz-zhalo,1+yhalo:ny-yhalo,n)
    enddo
!$OMP END DO
!$OMP END PARALLEL


    deallocate(weights)
    deallocate(data)

  end subroutine




end module
