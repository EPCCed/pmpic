!writes grid information
!just dumps x, y, z and the tag
module writegrids_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, LONG_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer,options_get_logical
  use prognostics_mod, only : prognostic_field_type
  use conversions_mod, only : conv_to_string
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use parcel_interpolation_mod, only : x_coords, y_coords, z_coords, par2grid, cache_parcel_interp_weights
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use netcdf, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       NF90_COLLECTIVE, nf90_def_var, nf90_var_par_access, nf90_def_var_fill, nf90_put_att, nf90_create, nf90_put_var, &
       nf90_def_dim, nf90_enddef, nf90_close, nf90_inq_dimid, nf90_inq_varid,&
       nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_noerr, nf90_strerror
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use logging_mod, only : LOG_ERROR, log_log
  use mpi, only : MPI_INFO_NULL

  implicit none

  character(len=*), parameter :: CHECKPOINT_TITLE = "MONC checkpoint file" !< Title of the NetCDF file

  integer :: num, gpersteps, gwritten

  integer :: handleg
                          
  character(len=*), parameter ::  U_KEY = "u", & 
                                  V_KEY = "v", &
                                  W_KEY = "w", &      
                                  P_KEY = "p", &      
                                  Q_KEY = "q", &      
                                  R_KEY = "r", &
                                  B_KEY = "b", &
                                  HG_KEY = "hg", &
                                  HGLIQ_KEY = "hgliq", &
                                  VOL_KEY = "vol", &
                                  X_KEY = "x", &
                                  Y_KEY = "y", &
                                  Z_KEY = "z", &
                                  TIMESTEP="timestep", &        
                                  TIME_KEY="time",&
                                  DTM_KEY="dtm",&
                                  CREATED_ATTRIBUTE_KEY="created",&
                                  TITLE_ATTRIBUTE_KEY="title"

  public write_checkpoint_file
  
  
contains

  type(component_descriptor_type) function writegrids_get_descriptor()
    writegrids_get_descriptor%name="writegrids"
    writegrids_get_descriptor%version=0.1
    writegrids_get_descriptor%initialisation=>initialisation_callback
    writegrids_get_descriptor%timestep=>timestep_callback
    writegrids_get_descriptor%finalisation=>finalisation_callback
  end function writegrids_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) print *, "Writer initialisation"

    gpersteps=options_get_integer(state%options_database,"grid_dump_frequency")

    num=state%iterations
    gwritten=0

    call register_routine_for_timing("write_grids",handleg,state)

  end subroutine


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    character (len=18) :: filename
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels
    integer :: i,j,k

    num=state%iterations

    if (gpersteps .ne. 0) then
      if (mod(num,gpersteps) .eq. 0) then


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


        call timer_start(handleg)
        write(filename,"(A,i3.3,A1,I4.4,A4)") "grids_", num, ".dat"
        call write_checkpoint_file(state, filename)
        gwritten=gwritten+1
        call timer_stop(handleg)
      endif
    endif

    num=num+1
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
  subroutine write_checkpoint_file(current_state, filename)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: filename

    integer :: ncid,u_id,v_id,w_id,p_id,q_id,r_id,b_id,hg_id,hgliq_id,vol_id,&
               x_dim_id,y_dim_id,z_dim_id,timestep_id,time_id,dtm_id


    call check_status(nf90_create(filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
         comm = current_state%parallel%monc_communicator, info = MPI_INFO_NULL))
    call write_out_global_attributes(ncid)

    ! define dimensions
    call check_status(nf90_def_dim(ncid, Z_KEY, current_state%global_grid%size(Z_INDEX), z_dim_id))
    call check_status(nf90_def_dim(ncid, Y_KEY, current_state%global_grid%size(Y_INDEX), y_dim_id))
    call check_status(nf90_def_dim(ncid, X_KEY, current_state%global_grid%size(X_INDEX), x_dim_id))
    
    !define prognostic variables
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=U_KEY, field_id=u_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=V_KEY, field_id=v_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=W_KEY, field_id=w_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=P_KEY, field_id=p_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=Q_KEY, field_id=q_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=R_KEY, field_id=r_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=B_KEY, field_id=b_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=HG_KEY, field_id=hg_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=HGLIQ_KEY, field_id=hgliq_id,field_units="-")
    call define_3d_variable(ncid, z_dim_id, y_dim_id, x_dim_id, field_name=VOL_KEY, field_id=vol_id,field_units="-")
    
    call define_misc_variables(ncid, timestep_id, time_id, dtm_id)

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
      call write_out_misc_variables(current_state, ncid, timestep_id, time_id,dtm_id)
    end if
    
    call check_status(nf90_close(ncid))
    
  end subroutine write_checkpoint_file

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

    call check_status(nf90_put_var(ncid, variable_id, field%data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
         local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
         local_grid%local_domain_end_index(X_INDEX)), start=start, count=count))

  end subroutine write_out_velocity_field

  !> Defines misc variables in the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable
  subroutine define_misc_variables(ncid, timestep_id, time_id, dtm_id)
    integer, intent(in) :: ncid
    integer, intent(out) :: timestep_id, time_id, dtm_id

    call check_status(nf90_def_var(ncid, TIMESTEP, NF90_INT, timestep_id))
    call check_status(nf90_def_var(ncid, TIME_KEY, NF90_DOUBLE, time_id))
    call check_status(nf90_def_var(ncid, DTM_KEY, NF90_DOUBLE, dtm_id))

  end subroutine define_misc_variables

  !> Will dump out (write) misc model data to the checkpoint
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable id
  subroutine write_out_misc_variables(current_state, ncid, timestep_id, time_id, dtm_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, timestep_id, time_id, dtm_id

    call check_status(nf90_put_var(ncid, timestep_id, current_state%timestep))
    call check_status(nf90_put_var(ncid, time_id, current_state%time+current_state%dtm))
    call check_status(nf90_put_var(ncid, dtm_id, current_state%dtm))

  end subroutine write_out_misc_variables
  
  subroutine define_3d_variable(ncid, dimone, dimtwo, dimthree, field_name, field_id, field_units)
    integer, intent(in) :: ncid, dimone, dimtwo, dimthree
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: field_units

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(3))
    dimids = (/ dimone, dimtwo, dimthree /)

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
    call check_status(nf90_var_par_access(ncid, field_id, NF90_COLLECTIVE))
    call check_status(nf90_put_att(ncid, field_id, "units", field_units))

  end subroutine define_3d_variable

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
