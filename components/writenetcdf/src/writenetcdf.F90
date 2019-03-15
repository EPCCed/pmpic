!writes grid information
!just dumps x, y, z and the tag
module writenetcdf_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, LONG_INTEGER, DOUBLE_PRECISION, STRING_LENGTH
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer,options_get_logical,options_get_string,options_get_real
  use prognostics_mod, only : prognostic_field_type
  use conversions_mod, only : conv_to_string
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use parcel_interpolation_mod, only : x_coords, y_coords, z_coords, par2grid, cache_parcel_interp_weights
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use netcdf, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       NF90_COLLECTIVE, nf90_def_var, nf90_var_par_access, nf90_def_var_fill, nf90_put_att, nf90_create, nf90_put_var, &
       nf90_def_dim, nf90_enddef, nf90_close, nf90_inq_dimid, nf90_inq_varid,&
       nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_noerr, nf90_strerror,NF90_UNLIMITED,nf90_def_var_deflate
  use logging_mod, only : LOG_ERROR, log_log
  use mpi, only : MPI_INFO_NULL

  implicit none

  character(len=*), parameter :: CHECKPOINT_TITLE = "MONC checkpoint file" !< Title of the NetCDF file

  integer :: num, ppersteps, gpersteps, pwritten, gwritten
 
  integer :: parcel_netcdf_deflate_level
  
  integer :: handlep, handleg, ierr

  real(kind=DEFAULT_PRECISION) :: dtparcels, dtgrids, tnextgrids, tnextparcels

  CHARACTER(len=4) :: mode
                          
  character(len=*), parameter ::  U_KEY = "u", & 
                                  V_KEY = "v", &
                                  W_KEY = "w", &      
                                  P_KEY = "p", &      
                                  Q_KEY = "q", &      
                                  R_KEY = "r", &
                                  B_KEY = "b", &
                                  H_KEY = "h", &
                                  HG_KEY = "hg", &
                                  HGLIQ_KEY = "hgliq", &
                                  VOL_KEY = "vol", &
                                  STRETCH_KEY = "stetch", &
                                  TAG_KEY = "tag", &
                                  QVALUES_KEY = "qvalues", &
                                  X_KEY = "x", &
                                  Y_KEY = "y", &
                                  Z_KEY = "z", &
                                  DXDT_KEY = "dxdt", &
                                  DYDT_KEY = "dydt", &
                                  DZDT_KEY = "dzdt", &
                                  DPDT_KEY = "dpdt", &
                                  DQDT_KEY = "dqdt", &
                                  DRDT_KEY = "drdt", &
                                  TIMESTEP="timestep", &        
                                  TIME_KEY="time",&
                                  NPARCEL_DIM_KEY="parcel",&
                                  QVAL_DIM_KEY="qvalue",&
                                  DTM_KEY="dtm",&
                                  CREATED_ATTRIBUTE_KEY="created",&
                                  TITLE_ATTRIBUTE_KEY="title",&
                                  NQVALUE_KEY="nqvalue",&
                                  NPARCEL_KEY="nparcel"

  public write_ncgrid_file
  
  
contains

  type(component_descriptor_type) function writenetcdf_get_descriptor()
    writenetcdf_get_descriptor%name="writenetcdf"
    writenetcdf_get_descriptor%version=0.1
    writenetcdf_get_descriptor%initialisation=>initialisation_callback
    writenetcdf_get_descriptor%timestep=>timestep_callback
    writenetcdf_get_descriptor%finalisation=>finalisation_callback
  end function writenetcdf_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) print *, "Writer initialisation"

    !determine the writing mode we want. Write a fixed number of timesteps ("steps") or a
    ! fixed time interval ("time")

    mode = options_get_string(state%options_database,"writenetcdf_mode")

    if (state%parallel%my_rank .eq. 0) then
      print *, "Selected writing mode: ",mode
    endif

    parcel_netcdf_deflate_level = options_get_real(state%options_database,"parcel_netcdf_deflate_level")
      
    if (mode .eq. "time") then
      dtparcels = options_get_real(state%options_database,"parcel_time_write_netcdf_frequency")
      dtgrids = options_get_real(state%options_database,"grid_time_write_netcdf_frequency")
      tnextgrids=0.
      tnextparcels=0.
    else if (mode .eq. "step") then
      ppersteps=options_get_integer(state%options_database,"parcel_step_write_netcdf_frequency")
      gpersteps=options_get_integer(state%options_database,"grid_step_write_netcdf_frequency")
    else if (mode .eq. "none") then
      if (state%parallel%my_rank .eq. 0) print *, "No grid/parcel files will be produced"
    else
      if (state%parallel%my_rank .eq. 0) then
        print *, "Error: mode '",mode,"' not recognised."
        print *, "The valid options are 'none', 'time' or 'step'"
        print *, 'Aborting'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    num=state%iterations
    pwritten=0
    gwritten=0

    call register_routine_for_timing("write_netcdf_parcels",handlep,state)
    call register_routine_for_timing("write_netcdf_grids",handleg,state)

  end subroutine


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    character (len=18) :: filename
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels
    integer :: i,j,k

    num=state%iterations

    if (mode .eq. "step") then

      if (ppersteps .ne. 0) then
        if(mod(num,ppersteps) .eq. 0) then
          call timer_start(handlep)
          call write_parcels_to_file(state)
          pwritten=pwritten+1
          call timer_stop(handlep)
        endif
      endif

      if (gpersteps .ne. 0) then
        if (mod(num,gpersteps) .eq. 0) then
          call timer_start(handleg)
          call write_grids_to_file(state)
          gwritten=gwritten+1
          call timer_stop(handleg)
        endif
      endif

    else if (mode .eq. "time") then
      if (dtgrids .ne. 0.) then
        if (state%time .ge. tnextgrids) then
          call timer_start(handleg)
          call write_grids_to_file(state)
          gwritten=gwritten+1
          call timer_stop(handleg)
          tnextgrids = tnextgrids + dtgrids
        endif
      endif

      if (dtparcels .ne. 0.) then
        if (state%time .ge. tnextparcels) then
          call timer_start(handlep)
          call write_parcels_to_file(state)
          pwritten=pwritten+1
          call timer_stop(handlep)
          tnextparcels = tnextparcels + dtparcels
        endif
      endif

    endif

    num=num+1
  end subroutine
    
  subroutine write_grids_to_file(state)
    type(model_state_type), intent(inout), target :: state
    character (len=23) :: filename
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels
    integer :: i,j,k
    
    call cache_parcel_interp_weights(state)
    call par2grid(state,state%parcels%b,state%b)
    call par2grid(state,state%parcels%p,state%p)
    call par2grid(state,state%parcels%q,state%q)
    call par2grid(state,state%parcels%r,state%r)

    ! obtain the humidity and liquid humidity
    call par2grid(state,state%parcels%h,state%hg)
    !$OMP PARALLEL DO
    do i=1,size(state%hgliq%data,3)
      do j=1,size(state%hgliq%data,2)
        do k=1,size(state%hgliq%data,1)
          state%hgliq%data(k,j,i) = max(0.,state%hg%data(k,j,i) - exp(-z_coords(k)))
          state%b%data(k,j,i) = state%b%data(k,j,i) + 12.5*state%hgliq%data(k,j,i)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    write(filename,"(A,I4.4,A3)") "grids_",num,".nc"
    call write_ncgrid_file(state, filename)
    gwritten=gwritten+1
  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) then
      print *, "Written", gwritten, "grid dumps"
    endif

  end subroutine

  !> Will write out the current model state_mod into a NetCDF checkpoint file
  !! @param currentState The current model state_mod
  !! @param filename The filename of the NetCDF file that will be written
  subroutine write_ncgrid_file(current_state, filename)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: filename

    integer :: ncid,u_id,v_id,w_id,p_id,q_id,r_id,b_id,hg_id,hgliq_id,vol_id,time_dim_id,&
               x_dim_id,y_dim_id,z_dim_id,timestep_id,time_id,dtm_id,x_id,y_id,z_id,i

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_arr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: y_arr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: x_arr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: time_arr

    call check_status(nf90_create(filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
         comm = current_state%parallel%monc_communicator, info = MPI_INFO_NULL))

    call write_out_global_attributes(ncid)

!    ! define dimensions
    call check_status(nf90_def_dim(ncid, TIME_KEY, 1, time_dim_id))
    call check_status(nf90_def_dim(ncid, Z_KEY, current_state%global_grid%size(Z_INDEX), z_dim_id))
    call check_status(nf90_def_dim(ncid, Y_KEY, current_state%global_grid%size(Y_INDEX), y_dim_id))
    call check_status(nf90_def_dim(ncid, X_KEY, current_state%global_grid%size(X_INDEX), x_dim_id))
    call define_1d_variable(ncid, time_dim_id, field_name=TIME_KEY, field_id=time_id,field_units="-")
    call define_1d_variable(ncid, z_dim_id, field_name=Z_KEY, field_id=z_id,field_units="-")
    call define_1d_variable(ncid, y_dim_id, field_name=Y_KEY, field_id=y_id,field_units="-")
    call define_1d_variable(ncid, x_dim_id, field_name=X_KEY, field_id=x_id,field_units="-")
    allocate(time_arr(1))
    allocate(z_arr(current_state%global_grid%size(Z_INDEX)))
    allocate(y_arr(current_state%global_grid%size(Y_INDEX)))
    allocate(x_arr(current_state%global_grid%size(X_INDEX)))
    do i=1,current_state%global_grid%size(Z_INDEX)
        z_arr(i)=current_state%global_grid%bottom(Z_INDEX)+current_state%global_grid%resolution(Z_INDEX)*(i-1.0)
    end do
    do i=1,current_state%global_grid%size(Y_INDEX)
        y_arr(i)=current_state%global_grid%bottom(Y_INDEX)+current_state%global_grid%resolution(Y_INDEX)*(i-0.5)
    end do
    do i=1,current_state%global_grid%size(X_INDEX)
        x_arr(i)=current_state%global_grid%bottom(X_INDEX)+current_state%global_grid%resolution(X_INDEX)*(i-0.5)
    end do
    time_arr(1)=current_state%time+current_state%dtm
    call check_status(nf90_put_var(ncid, time_id, time_arr))
    call check_status(nf90_put_var(ncid, z_id, z_arr))
    call check_status(nf90_put_var(ncid, y_id, y_arr))
    call check_status(nf90_put_var(ncid, x_id, x_arr))

!    !define prognostic variables
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=U_KEY, field_id=u_id,field_units="m s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=V_KEY, field_id=v_id,field_units="m s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=W_KEY, field_id=w_id,field_units="m s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=P_KEY, field_id=p_id,field_units="s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=Q_KEY, field_id=q_id,field_units="s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=R_KEY, field_id=r_id,field_units="s-1")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=B_KEY, field_id=b_id,field_units="m s-2")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, field_name=HG_KEY, field_id=hg_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, &
         field_name=HGLIQ_KEY, field_id=hgliq_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, time_dim_id, &
         field_name=VOL_KEY, field_id=vol_id,field_units="-")
    
    call define_misc_variables(ncid, timestep_id, dtm_id)

    call check_status(nf90_enddef(ncid))
    
    !write prognostic variables
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%u, u_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%v, v_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%w, w_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%p, p_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%q, q_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%r, r_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%b, b_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%hg, hg_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%hgliq, hgliq_id)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%vol, vol_id)

    if (current_state%parallel%my_rank==0) then
      call write_out_misc_variables(current_state, ncid, timestep_id, dtm_id)
    end if
    
    call check_status(nf90_close(ncid))

    deallocate(time_arr,z_arr,y_arr,x_arr)
        
  end subroutine write_ncgrid_file

  subroutine write_parcels_to_file(state)
    type(model_state_type), intent(inout), target :: state
    character (len=23) :: filename
    character (len=24) :: fnamedummy
    character (len=20) :: fieldnamedummy
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels
    integer :: nqvalues    
    integer :: ncid,p_id,q_id,r_id,b_id,h_id,vol_id,tag_id,stretch_id,time_dim_id,&
               x_id,y_id,z_id,timestep_id,time_id,dtm_id,i,&
               dxdt_id,dydt_id,dzdt_id,dpdt_id,dqdt_id,drdt_id,nqvalues_dim_id,&
               nparcels_dim_id,nparcel_id,nqvalue_id,&
               x_start_id,y_start_id,z_start_id,&
               x_end_id,y_end_id,z_end_id,dummy_id,&
               nparcels_id

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: nparcel_arr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: time_arr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: qvalues_id

    proc=state%parallel%my_rank
    nparcels=state%parcels%numparcels_local
    nqvalues=state%parcels%qnum

    write(filename,"(A8,i5.5,A1,I5.5,A3)") "parcels_", proc,"_",num,".nc"
    write(fnamedummy,"(A,I5.5,A3)") "parcels_[rank]_", num,".nc"

    if (proc .eq. 0) print *, "Writing parcels to '",fnamedummy,"'"

    !file contains:
    ! t
    ! xmin, xmax, ymin, ymax, zmin, zmax
    ! n_parcels
    ! data
    
    call check_status(nf90_create(filename, NF90_NETCDF4, ncid))
         
    call write_out_global_attributes(ncid)

!    ! define dimensions
    call check_status(nf90_def_dim(ncid, TIME_KEY, 1, time_dim_id))
    call check_status(nf90_def_dim(ncid, NPARCEL_KEY, NF90_UNLIMITED, nparcels_dim_id))
    call define_1d_variable(ncid, time_dim_id, field_name=TIME_KEY, field_id=time_id,field_units="-")
    call define_1d_variable(ncid, nparcels_dim_id, field_name=NPARCEL_DIM_KEY, field_id=nparcel_id,field_units="-")
    allocate(time_arr(1))
    allocate(nparcel_arr(nparcels))
    allocate(qvalues_id(nqvalues))

    do i=1,nparcels
        nparcel_arr(i)=i
    end do
    
    time_arr(1)=state%time+state%dtm
    call check_status(nf90_put_var(ncid, time_id, time_arr))
    call check_status(nf90_put_var(ncid, nparcel_id, nparcel_arr))
    
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=X_KEY, field_id=x_id,field_units="m")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=Y_KEY, field_id=y_id,field_units="m")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=Z_KEY, field_id=z_id,field_units="m")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=P_KEY, field_id=p_id,field_units="s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=Q_KEY, field_id=q_id,field_units="s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=R_KEY, field_id=r_id,field_units="s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DXDT_KEY, field_id=dxdt_id,field_units="m s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DYDT_KEY, field_id=dydt_id,field_units="m s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DZDT_KEY, field_id=dzdt_id,field_units="m s-1")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DPDT_KEY, field_id=dpdt_id,field_units="s-2")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DQDT_KEY, field_id=dqdt_id,field_units="s-2")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=DRDT_KEY, field_id=drdt_id,field_units="s-2")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=H_KEY, field_id=h_id,field_units="-")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=B_KEY, field_id=b_id,field_units="m s-2")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=VOL_KEY, field_id=vol_id,field_units="-")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=STRETCH_KEY, field_id=stretch_id,field_units="-")
    call define_parcel_variable(ncid, nparcels_dim_id, field_name=TAG_KEY, field_id=tag_id,field_units="-")
    
    do i=1,nqvalues
      fieldnamedummy=QVALUES_KEY//"_"//trim(conv_to_string(i))
      call define_parcel_variable(ncid, nparcels_dim_id, field_name=fieldnamedummy, field_id=dummy_id,field_units="-")
      qvalues_id(i)=dummy_id
    end do

    ! STILL EXPORT THESE AS WELL        
    call check_status(nf90_def_var(ncid, "x_start", NF90_DOUBLE, x_start_id))
    call check_status(nf90_def_var(ncid, "y_start", NF90_DOUBLE, y_start_id))
    call check_status(nf90_def_var(ncid, "z_start", NF90_DOUBLE, z_start_id))
    call check_status(nf90_def_var(ncid, "x_end", NF90_DOUBLE, x_end_id))
    call check_status(nf90_def_var(ncid, "y_end", NF90_DOUBLE, y_end_id))
    call check_status(nf90_def_var(ncid, "z_end", NF90_DOUBLE, z_end_id))

    call check_status(nf90_put_var(ncid, x_start_id, state%local_grid%local_domain_start_index(3)))
    call check_status(nf90_put_var(ncid, x_end_id, state%local_grid%local_domain_end_index(3)+1))
    call check_status(nf90_put_var(ncid, y_start_id, state%local_grid%local_domain_start_index(2)))
    call check_status(nf90_put_var(ncid, y_end_id, state%local_grid%local_domain_end_index(2)+1))
    call check_status(nf90_put_var(ncid, z_start_id, state%local_grid%local_domain_start_index(1)))
    call check_status(nf90_put_var(ncid, z_end_id, state%local_grid%local_domain_end_index(1)))

    call define_misc_variables(ncid, timestep_id, dtm_id)

    call check_status(nf90_enddef(ncid))

    call check_status(nf90_put_var(ncid, x_id,state%parcels%x(1:nparcels)))
    call check_status(nf90_put_var(ncid, y_id,state%parcels%y(1:nparcels)))
    call check_status(nf90_put_var(ncid, z_id,state%parcels%z(1:nparcels)))
    call check_status(nf90_put_var(ncid, p_id,state%parcels%p(1:nparcels)))
    call check_status(nf90_put_var(ncid, q_id,state%parcels%q(1:nparcels)))
    call check_status(nf90_put_var(ncid, r_id,state%parcels%r(1:nparcels)))
    call check_status(nf90_put_var(ncid, dxdt_id, state%parcels%dxdt(1:nparcels)))
    call check_status(nf90_put_var(ncid, dydt_id, state%parcels%dydt(1:nparcels)))
    call check_status(nf90_put_var(ncid, dzdt_id, state%parcels%dzdt(1:nparcels)))
    call check_status(nf90_put_var(ncid, dpdt_id, state%parcels%dpdt(1:nparcels)))
    call check_status(nf90_put_var(ncid, dqdt_id, state%parcels%dqdt(1:nparcels)))
    call check_status(nf90_put_var(ncid, drdt_id, state%parcels%drdt(1:nparcels)))
    call check_status(nf90_put_var(ncid, h_id, state%parcels%h(1:nparcels)))
    call check_status(nf90_put_var(ncid, b_id, state%parcels%b(1:nparcels)))
    call check_status(nf90_put_var(ncid, vol_id, state%parcels%vol(1:nparcels)))
    call check_status(nf90_put_var(ncid, stretch_id, state%parcels%stretch(1:nparcels)))
    call check_status(nf90_put_var(ncid, tag_id, state%parcels%tag(1:nparcels)))
    
    do i=1,nqvalues
       dummy_id=qvalues_id(i)
       call check_status(nf90_put_var(ncid, dummy_id, state%parcels%qvalues(i,1:nparcels)))
    end do
    
    call write_out_misc_variables(state, ncid, timestep_id, dtm_id)
    
    call check_status(nf90_close(ncid))

    deallocate(time_arr,nparcel_arr,qvalues_id)

    
  end subroutine
  
  !> Writes out global attributes into the checkpoint
  !! @param ncid NetCDF file id
  subroutine write_out_global_attributes(ncid)
    integer, intent(in) :: ncid

    integer :: date_values(8)

    call date_and_time(values=date_values)

    call check_status(nf90_put_att(ncid, NF90_GLOBAL, TITLE_ATTRIBUTE_KEY, CHECKPOINT_TITLE))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, CREATED_ATTRIBUTE_KEY, trim(conv_to_string(date_values(3)))//"/"//&
         trim(conv_to_string(date_values(2)))//"/"//trim(conv_to_string(date_values(1)))//" "//trim(conv_to_string(&
         date_values(5)))// ":"//trim(conv_to_string(date_values(6)))//":"//trim(conv_to_string(date_values(7)))))
  end subroutine write_out_global_attributes


  !> Will write out a single velocity field to the checkpoint file. If there are multiple processes then will determine
  !! the bounds, otherwise for serial just dump data field
  !! @param ncid The NetCDF file id
  !! @param field The model prognostic field to write out
  !! @param variable_id The NetCDF variable dimension id
  subroutine write_out_velocity_field(ncid, local_grid, field, variable_id)
    integer, intent(in) :: ncid, variable_id
    type(prognostic_field_type), intent(in) :: field
    type(local_grid_type), intent(inout) :: local_grid

    integer :: start(4), count(4), i, map(3)

    do i=1,3
      if (i==1) then
        map(i)=1
      else
        map(i)=map(i-1)*local_grid%size(i-1)
      end if
      start(i) = local_grid%start(i)
      count(i) = local_grid%size(i)
    end do
    start(4)=1
    count(4)=1

    call check_status(nf90_put_var(ncid, variable_id, field%data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
         local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
         local_grid%local_domain_end_index(X_INDEX)), start=start, count=count))

  end subroutine write_out_velocity_field

  !> Defines misc variables in the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable
  subroutine define_misc_variables(ncid, timestep_id, dtm_id)
    integer, intent(in) :: ncid
    integer, intent(out) :: timestep_id, dtm_id

    call check_status(nf90_def_var(ncid, TIMESTEP, NF90_INT, timestep_id))
    call check_status(nf90_def_var(ncid, DTM_KEY, NF90_DOUBLE, dtm_id))

  end subroutine define_misc_variables

  !> Will dump out (write) misc model data to the checkpoint
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable id
  subroutine write_out_misc_variables(current_state, ncid, timestep_id, dtm_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, timestep_id, dtm_id

    call check_status(nf90_put_var(ncid, timestep_id, current_state%timestep))
    call check_status(nf90_put_var(ncid, dtm_id, current_state%dtm))

  end subroutine write_out_misc_variables
  
  subroutine define_3d_variable(ncid, dimone, dimtwo, dimthree, dimt, field_name, field_id, field_units)
    integer, intent(in) :: ncid, dimone, dimtwo, dimthree, dimt
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: field_units

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(4))
    dimids = (/ dimone, dimtwo, dimthree, dimt /)

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
    call check_status(nf90_var_par_access(ncid, field_id, NF90_COLLECTIVE))
    call check_status(nf90_put_att(ncid, field_id, "units", field_units))

  end subroutine define_3d_variable

  subroutine define_1d_variable(ncid, dimone, field_name, field_id, field_units)
    integer, intent(in) :: ncid, dimone
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: field_units

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(1))
    dimids = (/ dimone /)

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
    call check_status(nf90_put_att(ncid, field_id, "units", field_units))

  end subroutine define_1d_variable

  subroutine define_parcel_variable(ncid, dimone, field_name, field_id, field_units)
    integer, intent(in) :: ncid, dimone
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: field_units

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(1))
    dimids = (/ dimone /)

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    call check_status(nf90_def_var_deflate(ncid, field_id, 0, 1, parcel_netcdf_deflate_level))
    call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
    call check_status(nf90_put_att(ncid, field_id, "units", field_units))

  end subroutine define_parcel_variable
  
  subroutine define_parcel_q_variable(ncid, dimone, dimtwo, field_name, field_id, field_units)
    integer, intent(in) :: ncid, dimone, dimtwo
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: field_units

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(2))
    dimids = (/ dimone, dimtwo /)

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    call check_status(nf90_def_var_deflate(ncid, field_id, 0, 1, parcel_netcdf_deflate_level))
    call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
    call check_status(nf90_put_att(ncid, field_id, "units", field_units))

  end subroutine define_parcel_q_variable
   
  !> Will check a NetCDF status and write to log_log error any decoded statuses. Can be used to decode
  !! whether a dimension or variable exists within the NetCDF data file
  !! @param status The NetCDF status flag
  !! @param foundFlag Whether the field has been found or not
  subroutine check_status(status, found_flag)
    integer, intent(in) :: status
    logical, intent(out), optional :: found_flag

    if (present(found_flag)) then
      found_flag = status /= nf90_ebaddim .and. status /= nf90_enotatt .and. status /= nf90_enotvar
      if (.not. found_flag) return
    end if

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status

end module
