!reads in parcel options from config file and allocates memory
!also places uniformly placed parcels in cells
module plume_parcelsetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, MPI_PARCEL_INT, STRING_LENGTH
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array, options_get_string
  use parcel_interpolation_mod, only:  x_coords, y_coords, z_coords
  use MPI
  use parcel_haloswap_mod, only: initialise_parcel_haloswapping


  implicit none

  integer(kind=PARCEL_INTEGER) :: maxparcels_global, maxparcels_local
  integer :: nprocs
  integer :: myrank
  integer :: n_per_dir
  integer :: ierr


!nicked from parameters.f90 in mpic for now...
  !------------------------------------------------------------------------
 !Note: we take the characteristic scale height 1/lambda = 1, so adjust
 !ellz above as needed.

 !Dimensionless latent heat of condensation, L*h_0/(c_p*theta_l0),
 !where h_0 is the saturation specific humidity at ground level,
 !c_p is the specific heat at constant pressure, and theta_l0 is
 !the mean liquid-water potential temperature:
REAL(KIND=DEFAULT_PRECISION),parameter:: latent=0.125

 !This is obtained from taking L/c_p=2500, h_0=0.015 and theta_l0=300.

 !A characteristic variation of the liquid water potential temperature
 !over a unit scale height, divided by its mean value theta_l0:
REAL(KIND=DEFAULT_PRECISION),parameter:: fdthetal0=0.01

 !Take the mean bouyancy frequency g*fdthetal0*lambda = 1 (over a
 !characteristic scale height 1/lambda):
REAL(KIND=DEFAULT_PRECISION),parameter:: gravity=1./fdthetal0

 !The group of parameters g*L*h_0/(c_p*theta_l0):
REAL(KIND=DEFAULT_PRECISION),parameter:: glat=gravity*latent


integer, parameter :: n_per_cell_dir_plume = 4
integer, parameter :: n_per_cell_dir_bg = 2

contains

  type(component_descriptor_type) function plume_parcelsetup_get_descriptor()
    plume_parcelsetup_get_descriptor%name="plume_parcelsetup"
    plume_parcelsetup_get_descriptor%version=0.1
    plume_parcelsetup_get_descriptor%initialisation=>initialisation_callback
  end function plume_parcelsetup_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    real(kind=DEFAULT_PRECISION) :: lambda, rhb, z_c, mu, z_d, z_m, r_plume, e_values(3)
    real(kind=DEFAULT_PRECISION) :: h_pl, h_bg, z_b, dbdz, b_pl
    logical :: master
    real(kind=DEFAULT_PRECISION) :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION) :: x_c_pl, y_c_pl, z_c_pl !x, y and z centres of plume
    real(kind=DEFAULT_PRECISION) :: dx, dy, dz
    real(kind=DEFAULT_PRECISION) :: dxplume, dyplume, dzplume
    real(kind=DEFAULT_PRECISION) :: dxbg, dybg, dzbg
    real(kind=DEFAULT_PRECISION) :: vol
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: xp, yp, zp, x, y, z
    integer(kind=PARCEL_INTEGER) :: n, n_plume, n_bg



    call MPI_Barrier(state%parallel%monc_communicator,ierr)

    master = state%parallel%my_rank .eq. 0

    lambda=options_get_real(state%options_database,"lambda")

    if (master) write(*,"('Scale height= ',f7.2,'m')") lambda

    rhb=options_get_real(state%options_database,"H")

    if (rhb .gt. 1) then
      if (master) write(*,'("Error: Relative humidity fraction must be less than 1. The selected value is ", f6.3)'), rhb
      call MPI_Finalize(ierr)
      stop
    endif

    z_c=options_get_real(state%options_database,"z_c")

    h_pl=exp(-z_c)
    if (master) write(*,"('Humidity inside the plume is ',f6.3)") h_pl


    mu=options_get_real(state%options_database,"mu")

    if (mu .gt. 1 .or. mu .le. rhb) then
      if (master) write(*,'("Error: mu must be between H and 1. The selected value is ", f6.3)'), mu
      call MPI_Finalize(ierr)
      stop
    endif

    h_bg=mu*h_pl
    if (master) write(*,"('Background humidity is ',f6.3)") h_bg

    z_b=log(rhb/h_bg)
    if (master) write(*,"('Base of mixed layer is ',f6.3)") z_b


    z_d=options_get_real(state%options_database,"z_d")
    z_m=options_get_real(state%options_database,"z_m")

    dbdz=glat*(h_pl-exp(-z_m))/(z_m-z_d)
    if (master) write(*,"('The buoyancy frequency in the stratified zone is ',f6.3)") sqrt(dbdz)

    !Also obtain the plume liquid-water buoyancy (using also z_b):
    b_pl=dbdz*(z_d-z_b)
    if (master) write(*,'(a,f7.5)') '  The plume liquid water buoyancy b_pl = ',b_pl
    if (master) write(*,'(a,f7.5)') '  corresponding to (theta_l-theta_l0)/theta_l0 = ',b_pl/gravity



    r_plume=options_get_real(state%options_database,"r_plume")
    if (2.*r_plume .gt. z_b) then
      if (master) write(*,"('Error: Plume radius is too big. At most it can be ',f7.5)") z_b/2.
      call MPI_Finalize(ierr)
      stop
    endif

    call options_get_real_array(state%options_database,"e_values",e_values)



    !now we need to transform dimensionless lengths into dimensional using the scale length lambda

    z_c=z_c*lambda
    z_b=z_b*lambda
    z_d=z_d*lambda
    z_m=z_m*lambda
    r_plume=r_plume*lambda

    e_values=e_values/r_plume/r_plume

    if (master) then
      write(*,*) "Box layout:"
      write(*,*) "z_max=", state%global_grid%top(1)
      write(*,*) "z_m  =", z_m
      write(*,*) "z_d  =", z_d
      write(*,*) "z_c  =", z_c
      write(*,*) "z_b  =", z_b
      write(*,*) "zmin  =", state%global_grid%bottom(1)
      write(*,*) "top of plume    =", 2*r_plume
      write(*,*) "bottom of plume =", 0.
    endif

    call MPI_Barrier(state%parallel%monc_communicator,ierr)

    ! determine x, y and z limits of this rank's box
    nx=state%local_grid%size(3)
    ny=state%local_grid%size(2)
    nz=state%local_grid%size(1)

    xmin_local=x_coords(state%local_grid%halo_size(3)+1)
    xmax_local=x_coords(2*state%local_grid%halo_size(3)+nx-1)

    ymin_local=y_coords(state%local_grid%halo_size(2)+1)
    ymax_local=y_coords(2*state%local_grid%halo_size(2)+ny-1)

    zmin_local=z_coords(1)
    zmax_local=z_coords(nz)



    !determine centre of plume (defined as the centre of the global computational domain)

    x_c_pl = (state%global_grid%top(3)-state%global_grid%bottom(3))/2 +state%global_grid%bottom(3)
    y_c_pl = (state%global_grid%top(2)-state%global_grid%bottom(2))/2 +state%global_grid%bottom(2)
    z_c_pl = r_plume

    !get grid spacings
    dx = state%global_grid%resolution(3)
    dy = state%global_grid%resolution(2)
    dz = state%global_grid%resolution(1)

    !determine spacing of parcels in plume
    dxplume = dx / n_per_cell_dir_plume
    dyplume = dy / n_per_cell_dir_plume
    dzplume = dz / n_per_cell_dir_plume


    !create parcels in plume
    n=0

    open(unit=10+state%parallel%my_rank)


    vol=(dxplume*dyplume*dzplume)/(dx*dy*dz)


    if (xmax_local .gt. x_c_pl-r_plume .and. xmin_local .lt. x_c_pl+r_plume ) then
      if (ymax_local .gt. y_c_pl-r_plume .and. ymin_local .lt. y_c_pl+r_plume ) then

        x=xmin_local+dxplume/2.
        y=ymin_local+dyplume/2.
        z=zmin_local+dzplume/2.
        do i=1,(nx)*n_per_cell_dir_plume
          xp=x-x_c_pl
          do j=1,(ny)*n_per_cell_dir_plume
            yp=y-y_c_pl
            do k=1,ceiling(2.*r_plume/dzplume)
              zp=z-z_c_pl
              !print *, x, y, z
              if (zp*zp + xp*xp + yp*yp .le. r_plume*r_plume) then
                n=n+1
                state%parcels%x(n) = x
                state%parcels%y(n) = y
                state%parcels%z(n) = z
                state%parcels%b(n) = b_pl*(1. + e_values(1)*xp*yp + e_values(2)*xp*zp + e_values(3)*yp*zp)
                state%parcels%h(n) = h_pl
                state%parcels%vol(n) = vol
                write(10+state%parallel%my_rank,*) x, y, z, state%parcels%b(n)
              endif
              z=z+dzplume
            enddo
            z=zmin_local+dzplume/2.
            y=y+dyplume
          enddo
          y=ymin_local+dyplume/2.
          x=x+dxplume
        enddo

      endif
    endif

    n_plume=n

    !create background parcels

    dxbg = dx/n_per_cell_dir_bg
    dybg = dy/n_per_cell_dir_bg
    dzbg = dz/n_per_cell_dir_bg

    vol = dxbg*dybg*dzbg/(dx*dy*dz)

    x=xmin_local+dxbg/2.
    y=ymin_local+dybg/2.
    z=zmin_local+dzbg/2.
    do i=1,nx*n_per_cell_dir_bg
      do j=1,ny*n_per_cell_dir_bg
        do k=1,(nz-1)*n_per_cell_dir_bg
          !print *, x, y, z

            if (z*z + y*y + x*x .gt. r_plume*r_plume) then

              n=n+1
              state%parcels%x(n) = x
              state%parcels%y(n) = y
              state%parcels%z(n) = z

              if (z .lt. z_b) then
            ! Mixed layer:
                state%parcels%b(n)=0.
                state%parcels%h(n)=h_bg
              else
            ! Stratified layer
                state%parcels%b(n)=dbdz/lambda*(z-z_b)
                state%parcels%h(n)=rhb*exp(-z/lambda)
              endif


              state%parcels%vol(n) = vol
              write(10+state%parallel%my_rank,*) x, y, z, state%parcels%b(n)
            endif

          z=z+dzbg
        enddo
        z=zmin_local+dzbg/2.
        y=y+dybg
      enddo
      y=ymin_local+dybg/2.
      x=x+dxbg
    enddo

    n_bg = n-n_plume


    close(10+state%parallel%my_rank)

    state%parcels%numparcels_local=n

    !update global parcel count
    call MPI_Allreduce(sendbuf=state%parcels%numparcels_local,&
                       recvbuf=state%parcels%numparcels_global,&
                       count=1,&
                       datatype=MPI_PARCEL_INT,&
                       op=MPI_SUM,&
                       comm=state%parallel%monc_communicator,&
                       ierror=ierr)



     call MPI_Allreduce(sendbuf=n_bg,&
                        recvbuf=n,&
                        count=1,&
                        datatype=MPI_PARCEL_INT,&
                        op=MPI_SUM,&
                        comm=state%parallel%monc_communicator,&
                        ierror=ierr)


    if (master) print *, "Total number of parcels=", state%parcels%numparcels_global
    if (master) print *, "background=", n
    if (master) print *, "plume=", state%parcels%numparcels_global-n






  end subroutine initialisation_callback









end module
