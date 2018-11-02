!contains grid2par and par2grid and helper subroutines
module parcel_interpolation_mod
  use state_mod, only: model_state_type
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use omp_lib
  use MPIC_Haloswap_mod, only: MPIC_Haloswap_init, grid2par_haloswap, par2grid_haloswap

  use optionsdatabase_mod, only : options_get_integer

  use communication_types_mod, only : halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, initiate_nonblocking_halo_swap, complete_nonblocking_halo_swap, &
       blocking_halo_swap, get_single_field_per_halo_cell,copy_corner_to_buffer,copy_buffer_to_corner
  use mpi
  use timer_mod


  implicit none

  type(halo_communication_type), save :: halo_swap_state


  !flag to see if this is intiialised or not
  logical :: initialised = .false.

  !cached variables
  integer, allocatable, dimension(:) ::  is, js, ks
  real(kind=DEFAULT_PRECISION), allocatable, dimension(:) :: delxs, delys, delzs
  integer(kind=PARCEL_INTEGER) :: nparcels

  !cached grid variables
  integer :: nx, ny, nz, ndz
  real(kind=DEFAULT_PRECISION) ::xmin, xmax, ymin, ymax, zmin, zmax
  real(kind=DEFAULT_PRECISION) :: maxx, maxy, maxz, minx, miny, minz
  real(kind=DEFAULT_PRECISION) :: dx, dy, dz
  integer :: hx, hy, hz
  !real(kind=DEFAULT_PRECISION), allocatable, dimension(:) :: z, dz

  real(kind=default_precision), allocatable, dimension(:) :: x_coords, y_coords, z_coords

  integer :: grid2par_handle, par2grid_handle, cache_handle, haloswap_handle, sumswap_handle


contains

  !initialise parcel interpolation data structures (this subroutine should be called once during parcel setup)
  !This contains variables for caching interpolation values as well as variables converting grid
  !indices to coordinates
  subroutine initialise_parcel_interp(state)
    type(model_state_type), intent(inout) :: state
    integer(kind=PARCEL_INTEGER) :: n
    integer :: xstart, xstop, ystart, ystop, zstart, zstop
    real(kind=DEFAULT_PRECISION) :: dzdummy
    integer :: halo_depth

    call MPIC_Haloswap_init(state)

    if (state%parallel%my_rank .eq. 0 ) print *, "Initialising parcel_interp module"
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
    dz=state%global_grid%resolution(1)

    !set number of cells in local grids
    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)


    if (dx .lt. 1.E-3 .or. dy .lt. 1.e-3) then
      error stop "no grid resolutions defined"
    endif


    !if (state%parallel%my_rank .eq. 0 ) print *, "dx=", dx, " dy=", dy, " dz=", dz

    !get global indices of the first and last elements in the local arrays
    xstart = state%local_grid%start(3)-state%local_grid%halo_size(3)
    xstop = state%local_grid%end(3)+state%local_grid%halo_size(3)

    ystart = state%local_grid%start(2)-state%local_grid%halo_size(2)
    ystop = state%local_grid%end(2)+state%local_grid%halo_size(2)

    zstart = state%local_grid%start(1)-state%local_grid%halo_size(1)
    zstop = state%local_grid%end(1)+state%local_grid%halo_size(1)


    ! !try to set up z
    ! if (allocated(state%global_grid%configuration%vertical%z)) then
    !   !if it is already set up then take its values from the state
    ! !  nz=size(state%global_grid%configuration%vertical%z)
    !   print*, "nz=",nz
    !   allocate(z(nz))
    !   z=state%global_grid%configuration%vertical%z
    !
    !
    !   ndz=size(state%global_grid%configuration%vertical%dz)
    !   print*, "dzn=",ndz
    !   allocate(dz(ndz))
    ! else !else we set it up ourselves
    !   if (state%parallel%my_rank .eq. 0 ) print *, "Warning: no z grid defined. Creating uniform grid"
    !
    !   !nz=zstop-zstart+1
    !   allocate(z(nz))
    !   ndz=nz-1
    !   allocate(dz(ndz))
    !
    !   !set dz
    !   dz(1:ndz)=dzdummy
    !
    !   meandz=sum(dz)/(nz-1)
    !
    !
    !   z(1)=(zstart-1)*dz(1)
    !   do n=2,nz
    !     z(n) = z(n-1)+dz(n-1)
    !   enddo
    ! endif

    !print*, state%parallel%my_rank, xstart, xstop, ystart, ystop, zstart, zstop
    ! do n=1,nz
    !   print*, state%parallel%my_rank, n, z(n)
    ! enddo

    xmin = (xstart-1)*dx + state%global_grid%bottom(3) !Coordinate of first point in the x grid (inc halo cells)
    ymin = (ystart-1)*dy + state%global_grid%bottom(2)
    zmin = (zstart-1)*dz + state%global_grid%bottom(1)

     !print *, "xmin, ymin, zmin", xmin, ymin, zmin
    !


    !coordinate of first point belonging to that grid
    minx = (xstart-1+state%local_grid%halo_size(3))*dx + state%global_grid%bottom(2)
    miny = (ystart-1+state%local_grid%halo_size(2))*dy + state%global_grid%bottom(2)
    minz = (zstart-1+state%local_grid%halo_size(1))*dz + state%global_grid%bottom(1)

    xmax = (xstop-1)*dx !coordinate of last point in the x grid (inc halo cells)
    ymax = (ystop-1)*dy
    zmax = (zstop-1)*dz

    !coordinate of last point belonging to the grid
    maxx= (xstop-1-state%local_grid%halo_size(3))*dx
    maxy= (ystop-1-state%local_grid%halo_size(2))*dy
    maxz= (zstop-1-state%local_grid%halo_size(1))*dz

    !allocate arrays that will hold the coordinates of x,y and z for each grid cell
    allocate(x_coords(nx), y_coords(ny))
    allocate(z_coords(nz))

    do n=1,nx
      x_coords(n) = xmin + (n-1)*dx + state%global_grid%bottom(3)
    enddo

    do n=1,ny
      y_coords(n) = ymin + (n-1)*dy + state%global_grid%bottom(2)
    enddo

    do n=1,nz
      z_coords(n) = zmin + (n-1)*dz + state%global_grid%bottom(1)
    enddo


    !z_coords(:) = z(:)

    !print *, "parcel_interp setup:", xmin, xmax, ymin, ymax, zmin, zmax

    !call flush()

    ! if (state%parallel%my_rank .eq. 0) then
    ! do n=1,nz
    !    print*, state%parallel%my_rank, n, z(n)
    ! enddo
    ! endif

    !sanity check to see if grid is set up right
    !print *, state%parallel%my_rank, ymin+dy, ymax-2*dy

    !call MPI_Barrier(state%parallel%monc_communicator,n)

    !print *, state%parallel%my_rank, y_coords(ny-2:ny), y_coords(1:3)

   !  call MPI_Barrier(state%parallel%monc_communicator,n)
     !call MPI_Finalize(n)
   ! !
     !stop

     hx=state%local_grid%halo_size(3)
     hy=state%local_grid%halo_size(2)
     hz=state%local_grid%halo_size(1)

     call register_routine_for_timing("grid2par",grid2par_handle,state)
     call register_routine_for_timing("par2grid",par2grid_handle,state)
     call register_routine_for_timing("cache_weights",cache_handle,state)
     call register_routine_for_timing("grid_HSwp", haloswap_handle, state)
     call register_routine_for_timing("grid_HSwp_sum",sumswap_handle, state)

     halo_depth = options_get_integer(state%options_database, "halo_depth")
     call init_halo_communication(state, get_single_field_per_halo_cell, halo_swap_state, &
          halo_depth, .true.)




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

    !deallocate(z)
    !deallocate(dz)

    deallocate(x_coords,y_coords,z_coords)

    initialised=.false.

    if (state%parallel%my_rank .eq. 0 ) print *, "parcel_interp finalised"

  end subroutine


  !cache weights for parcels (should be called once after parcel positions have changed before interpolating values)
  subroutine cache_parcel_interp_weights(state)
    type(model_state_type), intent(in) :: state
    integer(kind=PARCEL_INTEGER) :: n
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz
    real(kind=DEFAULT_PRECISION) :: xp, yp, zp

    call timer_start(cache_handle)


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
      i=floor((xp-xmin)/dx)+1
      j=floor((yp-ymin)/dy)+1
      ! do nn=1,nz-1 !as z may be a variable size grid we need to search through z to get the cell
      !   if ((zp .gt. z(nn)) .and. (zp .lt. z(nn+1))) then
      !     k=nn
      !     exit
      !   endif
      ! enddo
      k=floor((zp-zmin)/dz)+1

      !if ((xp .lt. xmin) .or. (xp .gt. xmax)) error stop "x too big/small"
      !if ((yp .lt. ymin) .or. (yp .gt. ymax)) error stop "y too big/small"
      !if ((zp .lt. zmin) .or. (zp .gt. zmax)) error stop "z too big/small"

      if ((i .lt. 1) .or. (i .ge. nx)) error stop "parcel out of box (x direction)"
      if ((j .lt. 1) .or. (j .ge. ny)) error stop "parcel out of box (y direction)"
      if ((k .lt. 1) .or. (k .ge. nz)) then
        print *, n, xp, yp, zp, zmin,zmax
        error stop "parcel out of box (z direction)"
      endif

      !get the fractional position in the cell (from the lower left corner) that the parcel is at
      delx= ((xp-xmin)-(i-1)*dx)/dx
      dely= ((yp-ymin)-(j-1)*dy)/dy
      delz= ((zp-zmin)-(k-1)*dz)/dz


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
    !if (state%parallel%my_rank .eq. 0) print *, "nthreads=", omp_get_num_threads()
    !$OMP END MASTER

    !$OMP END PARALLEL

    call timer_stop(cache_handle)

    !print *, "weights cached"

  end subroutine



  !interpolate gridded variable (grid) to parcel variable (var)
  subroutine grid2par(state,grid,var)
    type(model_state_type), intent(inout) :: state
    type(prognostic_field_type), intent(inout) :: grid
    real(kind=DEFAULT_PRECISION), dimension(:) :: var

    integer(kind=PARCEL_INTEGER) :: n
    real(kind=DEFAULT_PRECISION) :: c000, c001, c010, c011, c100, c101, c110, c111
    real(kind=DEFAULT_PRECISION) :: c00, c01, c10, c11
    real(kind=DEFAULT_PRECISION) :: c0, c1
    real(kind=DEFAULT_PRECISION) :: c

    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz

    call timer_start(grid2par_handle)

    call grid2par_haloswap(state,grid%data)
    !call perform_halo_swap(state,grid%data,perform_sum=.false.)

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

    call timer_stop(grid2par_handle)

  end subroutine



  subroutine par2grid(state,var,grid)

    implicit none

    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: var
    type(prognostic_field_type), intent(inout) :: grid

    integer(kind=PARCEL_INTEGER) :: n
    double precision, allocatable, dimension(:,:,:) :: weights, data
    double precision :: w
    double precision :: v
    double precision :: delx, dely, delz
    double precision :: w000, w001, w010, w011, w100, w101, w110, w111
    integer :: i,j,k
    integer :: xhalo, yhalo, zhalo

    call timer_start(par2grid_handle)

    allocate(weights(nz,ny,nx))
    allocate(data(nz,ny,nx))

    !zero the weights and the grid
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP              SHARED(nparcels,nx,ny,nz,state,var,grid, weights,data) &
!$OMP              SHARED(is, js, ks, delxs, delys, delzs) &
!$OMP SHARED(x_coords, y_coords, z_coords)

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

    !$OMP BARRIER

    !$OMP MASTER
    !call perform_halo_swap(state,weights,perform_sum=.true.)
    !call perform_halo_swap(state,data,perform_sum=.true.)
    call par2grid_haloswap(state,weights)
    call par2grid_haloswap(state,data)
    !$OMP END MASTER

    !$OMP BARRIER

    !$OMP WORKSHARE
    grid%data(1,:,:) = 2*data(1,:,:)
    grid%data(nz,:,:) = 2*data(nz,:,:)
    !$OMP END WORKSHARE


!$OMP DO
    do n=1,nx
        !set areas with weight=0 to 1 to prevent divide by 0
        where (weights(:,:,n) .eq. 0.) weights(:,:,n) = 1.
        grid%data(:,:,n) = data(:,:,n)/weights(:,:,n)
    enddo
!$OMP END DO

!$OMP END PARALLEL




    deallocate(weights)
    deallocate(data)

    call timer_stop(par2grid_handle)

  end subroutine




  !halo swapping functionality beyond this point



  !main halo swap routine to be called
  ! state - current model state
  ! data - 3D array of grid data
  ! perform_sum is an optional argument (default .false.)
  ! if .false. : Perform a "normal" haloswap (grid cells -> halo cells) using monc haloswapper
  ! if .true. : Take halo cells and add their contents to the grid cells (halo -> grid += halo)
    subroutine perform_halo_swap(state,data,perform_sum)
      type(model_state_type), intent(inout), target :: state
      real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:), target, intent(inout) :: data
      logical, intent(in), optional :: perform_sum

      logical :: sum
      real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:) :: left_buf, right_buf, up_buf, down_buf
      integer :: left, right, up, down
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr

      type(field_data_wrapper_type) :: source_data

      if (present(perform_sum)) then
        sum=perform_sum
      else
        sum = .false.
      endif

      source_data%data=>data

      !if we are summing then we need to move halo values across to the adjacent grid and add these to the
      !values already there. This means we cannot use MONC's inbuilt haloswapping (because it moves grid to halo
      ! not halo to grid) so we have our own version
      if (sum) then
        call timer_start(sumswap_handle)

        !allocate the buffers
        allocate(left_buf(nz,ny,hx), right_buf(nz,ny,hx), up_buf(nz,hy,nx), down_buf(nz,hy,nx))

        !determine the ranks in each direction from the "neighbours" array sotred in the state
        ! (seems to be neighbours(3,2*halosize) corresponding tot he 3 directions, and an entry for each halo cell)
        ! we assume that halo size is 2 and each halo cell in a given direction belongs to the same rank.
        down=state%local_grid%neighbours(2,1) ! in -y direction
        up=state%local_grid%neighbours(2,3) !in +y firection
        left=state%local_grid%neighbours(3,1) ! in -x direction
        right=state%local_grid%neighbours(3,3) ! in +x direction

        !send right, receive left

        right_buf(:,:,:) = data(:,:,nx-hx+1:nx) !copy halo to buffer

        call MPI_Sendrecv(right_buf, &
                          nz*ny*hx, &
                          PRECISION_TYPE, &
                          right, &
                          0, &
                          left_buf,&
                          nz*ny*hx,&
                          PRECISION_TYPE,&
                          left,&
                          0,&
                          state%parallel%monc_communicator,&
                          status,&
                          ierr)

        data(:,:,hx+1:hx+2) = data(:,:,hx+1:hx+2) + left_buf(:,:,:) !add buffer to grid


        ! !send left, receive right
        !
        ! left_buf(:,:,:) = data(:,:,1:hx) !copy halo to buffer
        !
        ! call MPI_Sendrecv(left_buf, &
        !                   nz*ny*hx, &
        !                   PRECISION_TYPE, &
        !                   left, &
        !                   0, &
        !                   right_buf,&
        !                   nz*ny*hx,&
        !                   PRECISION_TYPE,&
        !                   right,&
        !                   0,&
        !                   state%parallel%monc_communicator,&
        !                   status,&
        !                   ierr)
        !
        ! data(:,:,nx-2*hx+1:nx-hx) = data(:,:,nx-2*hx+1:nx-hx) + right_buf(:,:,:)

        !send up recv down

        up_buf(:,:,:) = data(:,ny-hy+1:ny,:) !copy halo to buffer

        call MPI_Sendrecv(up_buf, &
                          nz*nx*hy, &
                          PRECISION_TYPE, &
                          up, &
                          0, &
                          down_buf,&
                          nz*nx*hy,&
                          PRECISION_TYPE,&
                          down,&
                          0,&
                          state%parallel%monc_communicator,&
                          status,&
                          ierr)

         data(:,hy+1:2*hy,:) = data(:,hy+1:2*hy,:) + down_buf(:,:,:)

         ! !send down recv up
         !
         ! down_buf(:,:,:) = data(:,1:hy,:) !copy halo to buffer
         !
         ! call MPI_Sendrecv(down_buf, &
         !                   nz*nx*hy, &
         !                   PRECISION_TYPE, &
         !                   down, &
         !                   0, &
         !                   up_buf,&
         !                   nz*nx*hy,&
         !                   PRECISION_TYPE,&
         !                   up,&
         !                   0,&
         !                   state%parallel%monc_communicator,&
         !                   status,&
         !                   ierr)
         !
         !  data(:,ny-2*hy+1:ny-hy,:) = data(:,ny-2*hy+1:ny-hy,:) + up_buf(:,:,:)






         !then swap the halos conventionally to ensure the halos also have the correct values (do we need this?)

         call blocking_halo_swap(state, halo_swap_state, grid2buff, &
                                local_copy,buff2halo,&
                               copy_corners_to_halo_buffer=corner2buff,&
                               copy_from_halo_buffer_to_corner=buff2corner,&
                               source_data=(/source_data/))





        deallocate(up_buf, down_buf, left_buf, right_buf)
        call timer_stop(sumswap_handle)

      else
        call timer_start(haloswap_handle)
        call blocking_halo_swap(state, halo_swap_state, grid2buff, &
                               local_copy,buff2halo,&
                              copy_corners_to_halo_buffer=corner2buff,&
                              copy_from_halo_buffer_to_corner=buff2corner,&
                              source_data=(/source_data/))
        call timer_stop(haloswap_handle)

      endif


    end subroutine





    ! routines for interfacing with model_core halo swapping functionality

    subroutine grid2buff(current_state, neighbour_description, dim, source_index, &
         pid_location, current_page, source_data)
      type(model_state_type), intent(inout) :: current_state
      integer, intent(in) :: dim, pid_location, source_index
      integer, intent(inout) :: current_page(:)
      type(neighbour_description_type), intent(inout) :: neighbour_description
      type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

      type(field_data_wrapper_type) :: selected_source

      selected_source=source_data(1)

      call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, selected_source%data, &
           dim, source_index, current_page(pid_location))

      current_page(pid_location)=current_page(pid_location)+1
    end subroutine !copy_my_data_to_halo_buffer


    subroutine corner2buff(current_state, neighbour_description, &
         dim, x_source_index, &
         y_source_index, pid_location, current_page, source_data)
      !import model_state_type, neighbour_description_type, field_data_wrapper_type
      type(model_state_type), intent(inout) :: current_state
      integer, intent(in) :: dim, pid_location, x_source_index, y_source_index
      integer, intent(inout) :: current_page(:)
      type(neighbour_description_type), intent(inout) :: neighbour_description
      type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

      type(field_data_wrapper_type) :: selected_source

      selected_source=source_data(1)

      call copy_corner_to_buffer(current_state%local_grid, neighbour_description%send_corner_buffer,&
      selected_source%data, dim, &
           x_source_index, y_source_index, current_page(pid_location))

      current_page(pid_location)=current_page(pid_location)+1

    end subroutine


    subroutine local_copy(current_state, halo_depth, involve_corners, source_data)
     type(model_state_type), intent(inout) :: current_state
     integer, intent(in) :: halo_depth
     logical, intent(in) :: involve_corners
     type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

     type(field_data_wrapper_type) :: selected_source

     selected_source=source_data(1)

     call perform_local_data_copy_for_field(selected_source%data, current_state%local_grid, &
          current_state%parallel%my_rank, halo_depth, involve_corners)
   end subroutine !perform_local_data_copy_for_my_data


    subroutine buff2halo(current_state, neighbour_description, dim, target_index, &
         neighbour_location, current_page, source_data)
      type(model_state_type), intent(inout) :: current_state
      integer, intent(in) :: dim, target_index, neighbour_location
      integer, intent(inout) :: current_page(:)
      type(neighbour_description_type), intent(inout) :: neighbour_description
      type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

      type(field_data_wrapper_type) :: selected_source

      selected_source=source_data(1)

      call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, selected_source%data, &
           dim, target_index, current_page(neighbour_location))

      current_page(neighbour_location)=current_page(neighbour_location)+1
    end subroutine !copy_halo_buffer_to_my_data


    subroutine buff2corner(current_state, neighbour_description,&
         corner_loc, x_target_index, &
         y_target_index, neighbour_location, current_page, source_data)
    !  import model_state_type, neighbour_description_type, field_data_wrapper_type
      type(model_state_type), intent(inout) :: current_state
      integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
      integer, intent(inout) :: current_page(:)
      type(neighbour_description_type), intent(inout) :: neighbour_description
      type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

      type(field_data_wrapper_type) :: selected_source

      selected_source=source_data(1)

      call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer,&
       selected_source%data, corner_loc, &
           x_target_index, y_target_index, current_page(neighbour_location))

      current_page(neighbour_location) = current_page(neighbour_location)+1

    end subroutine buff2corner





end module
