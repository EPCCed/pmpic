!very basic parcel and grid writing routine
module writebinary_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, LONG_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer,options_get_logical,options_get_string,options_get_real
  use parcel_interpolation_mod, only : x_coords, y_coords, z_coords, par2grid, cache_parcel_interp_weights
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use q_indices_mod, only: get_q_index,standard_q_names

  implicit none

  integer :: num, ppersteps, gpersteps, pwritten, gwritten

  integer :: handlep, handleg, ierr

  real(kind=DEFAULT_PRECISION) :: dtparcels, dtgrids, tnextgrids, tnextparcels

  CHARACTER(len=4) :: mode

contains

  type(component_descriptor_type) function writebinary_get_descriptor()
    writebinary_get_descriptor%name="writebinary"
    writebinary_get_descriptor%version=0.1
    writebinary_get_descriptor%initialisation=>initialisation_callback
    writebinary_get_descriptor%timestep=>timestep_callback
    writebinary_get_descriptor%finalisation=>finalisation_callback
  end function writebinary_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state


    if (state%parallel%my_rank .eq. 0) print *, "Writer initialisation"

    !determine the writing mode we want. Write a fixed number of timesteps ("steps") or a
    ! fixed time interval ("time")

    mode = options_get_string(state%options_database,"writebinary_mode")

    if (state%parallel%my_rank .eq. 0) then
      print *, "Selected writing mode: ",mode
    endif

    if (mode .eq. "time") then
      dtparcels = options_get_real(state%options_database,"parcel_time_write_frequency")
      dtgrids = options_get_real(state%options_database,"grid_time_write_frequency")
      tnextgrids=0.
      tnextparcels=0.
    else if (mode .eq. "step") then
      ppersteps=options_get_integer(state%options_database,"parcel_step_write_frequency")
      gpersteps=options_get_integer(state%options_database,"grid_step_write_frequency")
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

    call register_routine_for_timing("write_parcels",handlep,state)
    call register_routine_for_timing("write_grids",handleg,state)



  end subroutine


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    character (len=23) :: filename
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

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) then
      print *, "Written", pwritten, "parcel dumps"
      print *, "Written", gwritten, "grid dumps"
    endif


  end subroutine



  subroutine write_parcels_to_file(state)
    type(model_state_type), intent(inout), target :: state
    character (len=23) :: filename
    character (len=24) :: fnamedummy
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels


    proc=state%parallel%my_rank
    nparcels=state%parcels%numparcels_local

    write(filename,"(A8,i5.5,A1,I5.5,A4)") "parcels_", proc,"_", num, ".dat"
    write(fnamedummy,"(A,I5.5,A4)") "parcels_[rank]_", num, ".dat"

    if (proc .eq. 0) print *, "Writing parcels to '",fnamedummy,"'"

    !file contains:
    ! t
    ! xmin, xmax, ymin, ymax, zmin, zmax
    ! n_parcels
    ! data
    open(unit=10,file=filename,form="unformatted",access="stream",status="replace")
    write(10) state%time

    write(10) x_coords(state%local_grid%local_domain_start_index(3))
    write(10) x_coords(state%local_grid%local_domain_end_index(3)+1)
    write(10) y_coords(state%local_grid%local_domain_start_index(2))
    write(10) y_coords(state%local_grid%local_domain_end_index(2)+1)
    write(10) z_coords(state%local_grid%local_domain_start_index(1))
    write(10) z_coords(state%local_grid%local_domain_end_index(1))

    write(10) nparcels

    write(10) state%parcels%x(1:nparcels)
    write(10) state%parcels%y(1:nparcels)
    write(10) state%parcels%z(1:nparcels)

    write(10) state%parcels%p(1:nparcels)
    write(10) state%parcels%q(1:nparcels)
    write(10) state%parcels%r(1:nparcels)

    write(10) state%parcels%dxdt(1:nparcels)
    write(10) state%parcels%dydt(1:nparcels)
    write(10) state%parcels%dzdt(1:nparcels)

    write(10) state%parcels%dpdt(1:nparcels)
    write(10) state%parcels%dqdt(1:nparcels)
    write(10) state%parcels%drdt(1:nparcels)

    write(10) state%parcels%h(1:nparcels)
    write(10) state%parcels%b(1:nparcels)
    write(10) state%parcels%vol(1:nparcels)

    write(10) state%parcels%stretch(1:nparcels)

    write(10) state%parcels%tag(1:nparcels)
    write(10) state%parcels%qvalues(:,1:nparcels)

    close(10)

  end subroutine

  subroutine write_grids_to_file(state)
    type(model_state_type), intent(inout) :: state
    character (len=21) :: filename
    character (len=22) :: fnamedummy
    integer :: proc
    integer :: i1, i2, j1, j2, k1, k2
    integer :: i, j, k, iqv, iqc
    real(kind=DEFAULT_PRECISION) :: x1, x2, y1, y2, z1, z2

    iqv=get_q_index(standard_q_names%VAPOUR, 'saturation_adjust')
    iqc=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'saturation_adjust')
    
    call cache_parcel_interp_weights(state)
    call par2grid(state,state%parcels%btot,state%btot)
    call par2grid(state,state%parcels%b,state%b)
    call par2grid(state,state%parcels%p,state%p)
    call par2grid(state,state%parcels%q,state%q)
    call par2grid(state,state%parcels%r,state%r)
       
    ! obtain the humidity and liquid humidity
    call par2grid(state,state%parcels%qvalues(iqv,:),state%hg)
    call par2grid(state,state%parcels%qvalues(iqc,:),state%hgliq)

    i1=state%local_grid%local_domain_start_index(1)
    i2=state%local_grid%local_domain_end_index(1)

    j1=state%local_grid%local_domain_start_index(2)
    j2=state%local_grid%local_domain_end_index(2)

    k1=state%local_grid%local_domain_start_index(3)
    k2=state%local_grid%local_domain_end_index(3)

    x1=x_coords(k1)
    x2=x_coords(k2+1)

    y1=y_coords(j1)
    y2=y_coords(j2+1)

    z1=z_coords(i1)
    z2=z_coords(i2)

    proc=state%parallel%my_rank

    write(filename,"(A,i5.5,A1,I5.5,A4)") "grids_", proc,"_", num, ".dat"
    write(fnamedummy,"(A,I5.5,A4)") "grids_[rank]_", num, ".dat"

    if (proc .eq. 0) print *, "Writing grids to '",fnamedummy,"'"

    ! print *, "nx, ny, nz=", state%local_grid%size(3), state%local_grid%size(2), state%local_grid%size(1)
    ! print *, "xmin, xmax=", x_coords(state%local_grid%local_domain_start_index(3)), &
    ! x_coords(state%local_grid%local_domain_end_index(3))
    ! print *, "ymin, ymax=", y_coords(state%local_grid%local_domain_start_index(2)), &
    ! y_coords(state%local_grid%local_domain_end_index(2))
    ! print *, "zmin, zmax=", z_coords(state%local_grid%local_domain_start_index(1)), &
    ! z_coords(state%local_grid%local_domain_end_index(1))


    !file contains:
    ! t
    ! xmin, xmax, ymin, ymax, zmin, zmax
    ! nx, ny, nz
    ! data
    open(unit=10,file=filename,form="unformatted",access="stream", status="replace")
    write(10) state%time
    write(10) x_coords(state%local_grid%local_domain_start_index(3))
    write(10) x_coords(state%local_grid%local_domain_end_index(3))
    write(10) y_coords(state%local_grid%local_domain_start_index(2))
    write(10) y_coords(state%local_grid%local_domain_end_index(2))
    write(10) z_coords(state%local_grid%local_domain_start_index(1))
    write(10) z_coords(state%local_grid%local_domain_end_index(1))
    write(10) state%local_grid%size(3), state%local_grid%size(2), state%local_grid%size(1)
    write(10) state%u%data(i1:i2, j1:j2, k1:k2)
    write(10) state%v%data(i1:i2, j1:j2, k1:k2)
    write(10) state%w%data(i1:i2, j1:j2, k1:k2)
    write(10) state%p%data(i1:i2, j1:j2, k1:k2)
    write(10) state%q%data(i1:i2, j1:j2, k1:k2)
    write(10) state%r%data(i1:i2, j1:j2, k1:k2)
    write(10) state%b%data(i1:i2, j1:j2, k1:k2)
    write(10) state%btot%data(i1:i2, j1:j2, k1:k2)
    write(10) state%hg%data(i1:i2, j1:j2, k1:k2)
    write(10) state%hgliq%data(i1:i2, j1:j2, k1:k2)
    write(10) state%vol%data(i1:i2, j1:j2, k1:k2)
    close(10)

  end subroutine


end module
