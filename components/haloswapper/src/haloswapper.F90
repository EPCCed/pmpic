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
       blocking_halo_swap, get_single_field_per_halo_cell,copy_corner_to_buffer,copy_buffer_to_corner
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE, mpi_barrier
  use optionsdatabase_mod, only : options_get_integer
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: my_data, append_data

  type(halo_communication_type), save :: halo_swap_state

  integer :: nx, ny, nz
  integer :: hx, hy, hz

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
    integer :: rank

    ! get halo_depth and pass it to the halo_swapping routines
    halo_depth = options_get_integer(current_state%options_database, "halo_depth")
    call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, &
         halo_depth, .true.)

    print*, "In haloswapper init"
    print*, current_state%local_grid%size(Z_INDEX), &
    current_state%local_grid%size(Y_INDEX),&
    current_state%local_grid%size(X_INDEX)

    nx=current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)
    ny=current_state%local_grid%size(Y_INDEX)+2*current_state%local_grid%halo_size(Y_INDEX)
    nz=current_state%local_grid%size(Z_INDEX)+2*current_state%local_grid%halo_size(Z_INDEX)

    !array to do normal halo swapping with - set to be rank in centre and 0 on halos
    allocate(my_data(nz, ny, nx))

    !array to do addition halo swapping - set to be rank for the whole array
    allocate(append_data(nz,ny,nx))

    rank=current_state%parallel%my_rank

    append_data(:,:,:) = rank



    hx=current_state%local_grid%halo_size(X_INDEX)
    hy=current_state%local_grid%halo_size(Y_INDEX)
    hz=current_state%local_grid%halo_size(Z_INDEX)

    my_data(:,:,:) = -1
    my_data(hz+1:nz-hz,hy+1:ny-hy,hx+1:nx-hx) = rank


  end subroutine initialisation_callback

  !> Timestep callback hook which performs the halo swapping for each prognostic field
  !!
  !! In parallel this is performed with MPI communication calls and wrapping around. In serial
  !! still need to wrap data around
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(field_data_wrapper_type) :: source_data

    integer :: rank, comm, ierr

    ! call initiate_nonblocking_halo_swap(current_state, halo_swap_state, &
    !      copy_my_data_to_halo_buffer, source_data=(/source_data/))

    ! Do something

    ! call complete_nonblocking_halo_swap(current_state, halo_swap_state, perform_local_data_copy_for_my_data, &
    !            copy_halo_buffer_to_my_data, source_data=(/source_data/))

    ! call blocking_halo_swap(current_state, halo_swap_state, halo2buff, &
    !                         local_copy,buff2halo,&
    !                        copy_corners_to_halo_buffer=corner2buff,&
    !                        copy_from_halo_buffer_to_corner=buff2corner,&
    !                        source_data=(/source_data/))


   !do a normal halo swap for my_data (halos are replaced by the neighbouring processes' values)
    call perform_halo_swap(current_state,my_data,perform_sum=.false.)


    rank=current_state%parallel%my_rank
    comm=current_state%parallel%monc_communicator
    call MPI_Barrier(comm,ierr)

    call flush()

    if (rank .eq. 0) print *, "Swap:"



    call MPI_Barrier(comm,ierr)

    write(*,"(i2.1,a,f3.0,f3.0, f3.0,f3.0)") rank, " halo    ", &
                        my_data(nz/2,ny/2,1), my_data(nz/2,ny/2,nx),&
                        my_data(nz/2,1,nx/2),my_data(nz/2,ny,nx/2)

    write(*,"(i2.1,a,f3.0,f3.0, f3.0,f3.0)") rank, " Corners ", &
                        my_data(nz/2,ny,nx), my_data(nz/2,1,nx),&
                        my_data(nz/2,1,1),my_data(nz/2,ny,1)


    call MPI_Barrier(comm, ierr)

    call flush()



    ! call blocking_halo_swap(current_state, halo_swap_state, halo2buff, &
    !                         local_copy,buff2halo,&
    !                        copy_corners_to_halo_buffer=corner2buff,&
    !                        copy_from_halo_buffer_to_corner=buff2corner,&
    !                        source_data=(/source_data/))

    !perform a summing halo swap - halos are set to halo + neighbouring processes data
    call perform_halo_swap(current_state,append_data,perform_sum=.true.)

    call MPI_Barrier(comm,ierr)

    if (rank .eq. 0) print *, "Add:"

    call MPI_Barrier(comm,ierr)

    write(*,"(i2.1,a,f3.0,f3.0, f3.0,f3.0)") rank, " halo    ", &
                        append_data(nz/2,ny/2,1), append_data(nz/2,ny/2,nx),&
                        append_data(nz/2,1,nx/2),append_data(nz/2,ny,nx/2)

    write(*,"(i2.1,a,f3.0,f3.0, f3.0,f3.0)") rank, " Corners ", &
                        append_data(nz/2,ny,nx), append_data(nz/2,1,nx),&
                        append_data(nz/2,1,1),append_data(nz/2,ny,1)


    call MPI_Barrier(comm, ierr)

    ! call MPI_Finalize(ierr)
    ! stop "finished"

  end subroutine timestep_callback

  !> The finalisation callback hook which will clean up and free the memory associated with the
  !! halo swapping
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(halo_swap_state)
    deallocate(my_data)
    deallocate(append_data)
  end subroutine finalisation_callback


!wrapper for intrnal halo swapping routines
! if the optional perform_sum argument is supplued (and true) then the halo values will have the
! neighbouring processes values added onto them. To do this we make a copy of the original halo
! values and add these back on after the swap
  subroutine perform_halo_swap(state,data,perform_sum)
    type(model_state_type), intent(inout), target :: state
    real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:), target, intent(inout) :: data
    logical, intent(in), optional :: perform_sum

    logical :: sum
    real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:,:) :: left, right, up, down


    type(field_data_wrapper_type) :: source_data

    if (present(perform_sum)) then
      sum=perform_sum
    else
      sum = .false.
    endif

    !if we are summing then we have to make a copy of our halos
    if (sum) then

      allocate(left(nz,ny,hx), right(nz,ny,hx), up(nz,hy,nx-2*hx), down(nz,hy,nx-2*hx))

      !copy data into our caches
      left(:,:,:) = data(:,:,1:hx)
      right(:,:,:) = data(:,:,nx-hx+1:nx)
      down(:,:,:) = data(:,1:hy,hx+1:nx-hx)
      up(:,:,:) = data(:,ny-hy+1:ny,hx+1:nx-hx)
    endif

    source_data%data=>data

    call blocking_halo_swap(state, halo_swap_state, halo2buff, &
                            local_copy,buff2halo,&
                           copy_corners_to_halo_buffer=corner2buff,&
                           copy_from_halo_buffer_to_corner=buff2corner,&
                           source_data=(/source_data/))

    if (sum) then
      !add our cache back onto halos

      data(:,:,1:hx) = data(:,:,1:hx) + left(:,:,:)
      data(:,:,nx-hx+1:nx) = data(:,:,nx-hx+1:nx) + right(:,:,:)
      data(:,1:hy,hx+1:nx-hx) = data(:,1:hy,hx+1:nx-hx)+ down(:,:,:)
      data(:,ny-hy+1:ny,hx+1:nx-hx) = data(:,ny-hy+1:ny,hx+1:nx-hx) +up(:,:,:)

      deallocate(up, down, left, right)

    endif

  end subroutine





  ! routines for interfacing with model_core halo swapping functionality

  subroutine halo2buff(current_state, neighbour_description, dim, source_index, &
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


end module haloswapper_mod
