!reads in parcel options from config file and allocates memory
!also places uniformly placed parcels in cells
module parcelsetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, MPI_PARCEL_INT, STRING_LENGTH
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array, options_get_string
  use parcel_interpolation_mod, only: initialise_parcel_interp, finalise_parcel_interp, x_coords, y_coords, z_coords
  use MPI
  use parcel_haloswap_mod, only: initialise_parcel_haloswapping

  use basicsetup_mod
  use readfromfile_mod


  implicit none

  integer(kind=PARCEL_INTEGER) :: maxparcels_global, maxparcels_local
  integer :: nprocs
  integer :: myrank
  integer :: n_per_dir
  integer :: ierr

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

    if (myrank .eq. 0) then
       print *, ""
       print *, "Parcel Setup:"
    endif

    !get options from config file
    maxparcels_global=options_get_integer(current_state%options_database,"max_parcels")

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

    !initialise parcel haloswapping
    call initialise_parcel_haloswapping(current_state)

    !set up the parcels. Depending on the config file, we either read in from
    !existing files, or set up new parcels according to a custom subroutine
    call setup_parcels(current_state)


  !  if (myrank .eq. 0) print *, "parcel setup done"


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



  !set up the parcels. Depending on the config file, we either read in from
  !existing files, or set up new parcels according to a custom subroutine
  subroutine setup_parcels(state)
    type(model_state_type), intent(inout) :: state
    logical :: restart
    character (len=STRING_LENGTH) :: setup_routine

    restart=options_get_logical(state%options_database,"restart")


    if (restart) then
      call read_parcels_from_file(state)
    else

      setup_routine=options_get_string(state%options_database,"initialisation_routine")

      if (setup_routine .eq. "basic") then
        call basicsetup(state)
      else
        print *, "Selected initialisation routine '",trim(setup_routine),"' not valid"
        call MPI_Finalize(ierr)
        error stop "Select a valid initilisation routine"
      endif

    endif


  end subroutine





  subroutine read_configuration(state)
    type(model_state_type), intent(inout) :: state

    maxparcels_global=options_get_integer(state%options_database,"max_parcels")
  !  n_per_dir=options_get_integer(state%options_database,"parcels_per_cell_dir")

  end subroutine read_configuration




end module
