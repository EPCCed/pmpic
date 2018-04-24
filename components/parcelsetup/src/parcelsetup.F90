module parcelsetup_mod
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array

  implicit none

  integer :: maxparcels_global, maxparcels_local
  integer :: nprocs
  integer :: myrank

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
    allocate(current_state%parcels%r(maxparcels_local))
    allocate(current_state%parcels%s(maxparcels_local))
    allocate(current_state%parcels%t(maxparcels_local))
    allocate(current_state%parcels%u(maxparcels_local))
    allocate(current_state%parcels%v(maxparcels_local))
    allocate(current_state%parcels%w(maxparcels_local))
    allocate(current_state%parcels%q(maxparcels_local))
    allocate(current_state%parcels%b(maxparcels_local))
    allocate(current_state%parcels%vol(maxparcels_local))


    if (myrank .eq. 0) print *, "parcel setup done"


  end subroutine initialisation_callback



  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    if (myrank .eq. 0) print *, "Parcel setup - finalise"
    if (myrank .eq. 0) print *, "deallocating parcels..."

    deallocate(current_state%parcels%x)
    deallocate(current_state%parcels%y)
    deallocate(current_state%parcels%z)
    deallocate(current_state%parcels%r)
    deallocate(current_state%parcels%s)
    deallocate(current_state%parcels%t)
    deallocate(current_state%parcels%u)
    deallocate(current_state%parcels%v)
    deallocate(current_state%parcels%w)
    deallocate(current_state%parcels%q)
    deallocate(current_state%parcels%b)
    deallocate(current_state%parcels%vol)

    if (myrank .eq. 0) print *, "done!"


  end subroutine finalisation_callback




  subroutine read_configuration(state)
    type(model_state_type), intent(inout) :: state

    maxparcels_global=options_get_integer(state%options_database,"max_parcels")

  end subroutine read_configuration






end module
