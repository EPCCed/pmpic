!> Performs halo swapping. In the parallel case this is between neighbouring processes
!! and in the serial case it still needs to wrap the halos around for the boundary conditions. This module
!! determines the policy of halo swapping (i.e. the fields to communicate) and the halo communication
!! module is used to provide the actual mechanism.
module haloswapper_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_DEBUG, log_get_logging_level, log_log
  use conversions_mod, only : conv_to_string
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, initiate_nonblocking_halo_swap, complete_nonblocking_halo_swap, &
       blocking_halo_swap, get_single_field_per_halo_cell
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  use optionsdatabase_mod, only : options_get_integer
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable, target :: my_data

  type(halo_communication_type), save :: halo_swap_state

  public haloswapper_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function haloswapper_get_descriptor()
    haloswapper_get_descriptor%name="halo_swapper"
    haloswapper_get_descriptor%version=0.1
    haloswapper_get_descriptor%initialisation=>initialisation_callback
    haloswapper_get_descriptor%timestep=>timestep_callback
    haloswapper_get_descriptor%finalisation=>finalisation_callback
  end function haloswapper_get_descriptor

  !> Initialisation callback hook which will set up the halo swapping state and cache some
  !!precalculated data for fast(er) halo swapping at each timestep
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: halo_depth
    integer :: nx, ny, nz

    ! get halo_depth and pass it to the halo_swapping routines
    halo_depth = options_get_integer(current_state%options_database, "halo_depth")
    call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, &
         halo_depth, .false.)

         print*, "In haloswapper init"
         print*, current_state%local_grid%size(Z_INDEX), &
         current_state%local_grid%size(Y_INDEX),&
         current_state%local_grid%size(X_INDEX)

         nx=current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)
         ny=current_state%local_grid%size(Y_INDEX)+2*current_state%local_grid%halo_size(Y_INDEX)
         nz=current_state%local_grid%size(Z_INDEX)+2*current_state%local_grid%halo_size(Z_INDEX)


    allocate(my_data(nz, ny, &
        nx))
  end subroutine initialisation_callback

  !> Timestep callback hook which performs the halo swapping for each prognostic field
  !!
  !! In parallel this is performed with MPI communication calls and wrapping around. In serial
  !! still need to wrap data around
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(field_data_wrapper_type) :: source_data
    source_data%data=>my_data

    call initiate_nonblocking_halo_swap(current_state, halo_swap_state, &
         copy_my_data_to_halo_buffer, source_data=(/source_data/))

    ! Do something

    call complete_nonblocking_halo_swap(current_state, halo_swap_state, perform_local_data_copy_for_my_data, &
               copy_halo_buffer_to_my_data, source_data=(/source_data/))
  end subroutine timestep_callback

  !> The finalisation callback hook which will clean up and free the memory associated with the
  !! halo swapping
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(halo_swap_state)
    deallocate(my_data)
  end subroutine finalisation_callback

   subroutine perform_local_data_copy_for_my_data(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    type(field_data_wrapper_type) :: selected_source

    selected_source=source_data(1)

    call perform_local_data_copy_for_field(selected_source%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_my_data

  subroutine copy_my_data_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
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
  end subroutine copy_my_data_to_halo_buffer

  subroutine copy_halo_buffer_to_my_data(current_state, neighbour_description, dim, target_index, &
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
  end subroutine copy_halo_buffer_to_my_data
end module haloswapper_mod
