!Creates an initial condition of a spherical thermal from Clark and Grabowski papers
module thermal_wojtek_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, MPI_PARCEL_INT, STRING_LENGTH
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array, options_get_string
  use parcel_interpolation_mod, only:  x_coords, y_coords, z_coords
  use MPI
  use science_constants_mod, only : G,rlvap,cp,thref0,pref0,surf_temp_base,surf_pres_base,&
  inv_exn_factor,r_over_cp,pi,lapse_ref,lapse_exp,p_ref_prefactor
  use saturation_mod, only: qsaturation
  use q_indices_mod, only: get_q_index,standard_q_names
  
  implicit none

  integer(kind=PARCEL_INTEGER) :: maxparcels_global, maxparcels_local
  integer :: nprocs
  integer :: myrank
  integer :: n_per_dir
  integer :: ierr

  integer, parameter :: n_per_cell_dir_thermal = 4
  integer, parameter :: n_per_cell_dir_bg = 2

contains

  type(component_descriptor_type) function thermal_wojtek_get_descriptor()
    thermal_wojtek_get_descriptor%name="thermal_wojtek"
    thermal_wojtek_get_descriptor%version=0.1
    thermal_wojtek_get_descriptor%initialisation=>initialisation_callback
  end function thermal_wojtek_get_descriptor

  ! Reference pressure, given surface temperature, surface pressure and a constant lapse rate
  pure real(kind=DEFAULT_PRECISION) function p_ref(z)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in)  :: z
    p_ref=p_ref_prefactor*((surf_temp_base-lapse_ref*z)**lapse_exp)
    return
  end function p_ref
  
  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    real(kind=DEFAULT_PRECISION) :: r_thermal,r_edge,rh_thermal,z_thermal,rh_env,dlnthdz, e_values(3)
    logical :: master
    real(kind=DEFAULT_PRECISION) :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION) :: x_c_thermal, y_c_thermal, z_c_thermal !x, y and z centres of thermal
    real(kind=DEFAULT_PRECISION) :: dx, dy, dz
    real(kind=DEFAULT_PRECISION) :: dxthermal, dythermal, dzthermal
    real(kind=DEFAULT_PRECISION) :: dxbg, dybg, dzbg
    real(kind=DEFAULT_PRECISION) :: vol
    integer :: i, j, k, iqc, iqv
    real(kind=DEFAULT_PRECISION) :: xp, yp, zp, x, y, z
    real(kind=DEFAULT_PRECISION) :: p,exn,theta,theta_surf,r_outside
    integer(kind=PARCEL_INTEGER) :: n, n_thermal, n_bg
    iqv=get_q_index(standard_q_names%VAPOUR, 'saturation_adjust')
    iqc=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'saturation_adjust')
    
    call MPI_Barrier(state%parallel%monc_communicator,ierr)

    master = state%parallel%my_rank .eq. 0

    r_thermal=options_get_real(state%options_database,"r_thermal")
    r_edge=options_get_real(state%options_database,"r_edge")
    rh_thermal=options_get_real(state%options_database,"rh_thermal")
    z_c_thermal=options_get_real(state%options_database,"z_c_thermal")
    rh_env=options_get_real(state%options_database,"rh_env")
    dlnthdz=options_get_real(state%options_database,"dlnthdz")
    call options_get_real_array(state%options_database,"e_values",e_values)

    e_values=e_values/r_thermal/r_thermal

    if (master) write(*,'(a,f11.5)') 'r_thermal = ',r_thermal
    if (master) write(*,'(a,f11.5)') 'r_edge = ',r_edge
    if (master) write(*,'(a,f11.5)') 'rh_thermal = ',rh_thermal
    if (master) write(*,'(a,f11.5)') 'z_c_thermal = ',z_c_thermal
    if (master) write(*,'(a,f11.5)') 'rh_env = ',rh_env
    if (master) write(*,'(a,f11.5)') 'dlnthdz = ',dlnthdz
    if (master) write(*,'(a,f11.5)') 'e_values(1) = ',e_values(1)
    if (master) write(*,'(a,f11.5)') 'e_values(2) = ',e_values(2)
    if (master) write(*,'(a,f11.5)') 'e_values(3) = ',e_values(3)

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



    !determine centre of thermal (defined as the centre of the global computational domain)

    !get grid spacings
    dx = state%global_grid%resolution(3)
    dy = state%global_grid%resolution(2)
    dz = state%global_grid%resolution(1)

    x_c_thermal = (state%global_grid%top(3)+dx-state%global_grid%bottom(3))/2 +state%global_grid%bottom(3)
    y_c_thermal = (state%global_grid%top(2)+dy-state%global_grid%bottom(2))/2 +state%global_grid%bottom(2)

    !determine spacing of parcels in thermal
    dxthermal = dx / n_per_cell_dir_thermal
    dythermal = dy / n_per_cell_dir_thermal
    dzthermal = dz / n_per_cell_dir_thermal


    !create parcels in thermal
    n=0

    !open(unit=10+state%parallel%my_rank)


    vol=(dxthermal*dythermal*dzthermal)
    theta_surf=surf_temp_base/(inv_exn_factor*surf_pres_base**r_over_cp)

    if (xmax_local .gt. x_c_thermal-r_thermal-r_edge .and. xmin_local .lt. x_c_thermal+r_thermal+r_edge) then
      if (ymax_local .gt. y_c_thermal-r_thermal-r_edge .and. ymin_local .lt. y_c_thermal+r_thermal+r_edge) then

        x=xmin_local+dxthermal/2.
        y=ymin_local+dythermal/2.
        z=zmin_local+dzthermal/2.
        do i=1,(nx)*n_per_cell_dir_thermal
          xp=x-x_c_thermal
          do j=1,(ny)*n_per_cell_dir_thermal
            yp=y-y_c_thermal
            do k=1,(nz-1)*n_per_cell_dir_thermal
              zp=z-z_c_thermal
              !print *, x, y, z
              if (zp*zp + xp*xp + yp*yp .le. r_thermal*r_thermal) then
                n=n+1
                if (n .gt. state%parcels%maxparcels_local) then
                  print *, "Error! Maxparcels reached in thermal_parcelsetup"
                  error stop "Maxparcels reached"
                endif
                state%parcels%x(n) = x
                state%parcels%y(n) = y
                state%parcels%z(n) = z
                p=p_ref(z)
                exn=inv_exn_factor*p**r_over_cp
                theta=theta_surf*exp(dlnthdz*z)
                state%parcels%b(n)=theta-thref0   
                state%parcels%qvalues(iqv,n) =  qsaturation(theta*exn,0.01_DEFAULT_PRECISION*p)*&
                (rh_env+(rh_thermal-rh_env)*(1. + e_values(1)*xp*yp + e_values(2)*xp*zp + e_values(3)*yp*zp))
                state%parcels%qvalues(iqc,n)=0.0_DEFAULT_PRECISION
                state%parcels%vol(n) = vol
!                write(10+state%parallel%my_rank,*) x, y, z, state%parcels%b(n)
              else if (zp*zp + xp*xp + yp*yp .le. (r_thermal+r_edge)*(r_thermal+r_edge)) then
                n=n+1
                if (n .gt. state%parcels%maxparcels_local) then
                  print *, "Error! Maxparcels reached in thermal_parcelsetup"
                  error stop "Maxparcels reached"
                endif
                state%parcels%x(n) = x
                state%parcels%y(n) = y
                state%parcels%z(n) = z
                p=p_ref(z)
                exn=inv_exn_factor*p**r_over_cp
                theta=theta_surf*exp(dlnthdz*z)
                state%parcels%b(n)=theta-thref0   
                r_outside=sqrt(zp*zp + xp*xp + yp*yp)-r_thermal
                state%parcels%qvalues(iqv,n) =  qsaturation(theta*exn,0.01_DEFAULT_PRECISION*p)*&
                (rh_env+(rh_thermal-rh_env)*(1. + e_values(1)*xp*yp + e_values(2)*xp*zp + e_values(3)*yp*zp)*&
                cos(2.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)*r_outside/r_edge)**2)
                state%parcels%qvalues(iqc,n)=0.0_DEFAULT_PRECISION
                state%parcels%vol(n) = vol
!                write(10+state%parallel%my_rank,*) x, y, z, state%parcels%b(n)
              end if
              z=z+dzthermal
            enddo
            z=zmin_local+dzthermal/2.
            y=y+dythermal
          enddo
          y=ymin_local+dythermal/2.
          x=x+dxthermal
        enddo

      endif
    endif

    n_thermal=n

    !create background parcels

    dxbg = dx/n_per_cell_dir_bg
    dybg = dy/n_per_cell_dir_bg
    dzbg = dz/n_per_cell_dir_bg

    vol = dxbg*dybg*dzbg

    x=xmin_local+dxbg/2.
    y=ymin_local+dybg/2.
    z=zmin_local+dzbg/2.
    do i=1,nx*n_per_cell_dir_bg
      xp=x-x_c_thermal
      do j=1,ny*n_per_cell_dir_bg
        yp=y-y_c_thermal
        do k=1,(nz-1)*n_per_cell_dir_bg
          zp=z-z_c_thermal
          !print *, x, y, z

            if (zp*zp + yp*yp + xp*xp .gt. (r_thermal+r_edge)*(r_thermal+r_edge)) then

              n=n+1
              if (n .gt. state%parcels%maxparcels_local) then
                print *, "Error! Maxparcels reached in thermal_parcelsetup"
                error stop "Maxparcels reached"
              endif
              state%parcels%x(n) = x
              state%parcels%y(n) = y
              state%parcels%z(n) = z

              p=p_ref(z)
              exn=inv_exn_factor*p**r_over_cp
              theta=theta_surf*exp(dlnthdz*z)
              state%parcels%b(n)=theta-thref0   
              state%parcels%qvalues(iqv,n)=qsaturation(theta*exn,0.01_DEFAULT_PRECISION*p)*rh_env
              state%parcels%qvalues(iqc,n)=0.0_DEFAULT_PRECISION
              state%parcels%vol(n) = vol
                
              !write(10+state%parallel%my_rank,*) x, y, z, state%parcels%b(n)
            endif

          z=z+dzbg
        enddo
        z=zmin_local+dzbg/2.
        y=y+dybg
      enddo
      y=ymin_local+dybg/2.
      x=x+dxbg
    enddo

    n_bg = n-n_thermal


    !close(10+state%parallel%my_rank)

    state%parcels%numparcels_local=n

    !update global parcel count
    call MPI_Allreduce(state%parcels%numparcels_local,&
                       state%parcels%numparcels_global,&
                       1,&
                       MPI_PARCEL_INT,&
                       MPI_SUM,&
                       state%parallel%monc_communicator,&
                       ierr)



     call MPI_Allreduce(n_bg,&
                        n,&
                        1,&
                        MPI_PARCEL_INT,&
                        MPI_SUM,&
                        state%parallel%monc_communicator,&
                        ierr)


    if (master) print *, "Total number of parcels=", state%parcels%numparcels_global
    if (master) print *, "background=", n
    if (master) print *, "thermal=", state%parcels%numparcels_global-n

  end subroutine initialisation_callback

end module
