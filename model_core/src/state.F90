!> The model state which represents the current state of a run
module state_mod
  use collections_mod, only : hashmap_type
  use grids_mod, only : global_grid_type, local_grid_type
  use communication_types_mod, only : halo_communication_type
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> The constants defining the reason why the model has terminated
  integer, parameter, public :: TIME_TERMINATION_REASON=0, TIMESTEP_TERMINATION_REASON=1, MESSAGE_TERMINATION_REASON=2, &
       WALLTIME_TERMINATION_REASON=3

  !> Information about the parallel aspects of the system
  type, public :: parallel_state_type
     integer :: processes, & !> Total number of processes
          my_rank,&           !> My process rank in the system
          neighbour_comm,&     !> Neighbour communicator
          monc_communicator=-1, io_communicator=-1, corresponding_io_server_process
     integer, dimension(3) :: &
          my_coords,&         !> My process coordinates in each dimension
          dim_sizes           !> Number of processes in each dimension
     logical, dimension(3,2) :: wrapped_around
     procedure(), nopass, pointer :: decomposition_procedure => null() !> The decomposition procedure to use
  end type parallel_state_type
  
  !> The ModelState which represents the current state of a run
  !!
  !! This state is provided to each callback and may be used and modified as required by
  !! the callbacks. Apart from this state, there should be no other state (global) variables
  !! declared. This allows us to simply persist and retrieve the ModelState when suspending
  !! and reactivating MONC.
  type, public :: model_state_type
    logical :: continue_timestep=.true., initialised=.false., continuation_run=.false.
    !logical :: use_viscosity_and_diffusion=.true., &
    !   use_surface_boundary_conditions=.true., backscatter=.true.

    type(hashmap_type) :: options_database
    type(global_grid_type) :: global_grid
    type(local_grid_type) :: local_grid
    type(parallel_state_type) :: parallel
    real(kind=DEFAULT_PRECISION) :: time=.0_DEFAULT_PRECISION,& ! Model time in seconds
            dtm,& ! Modeltimestep (s)
            absolute_new_dtm, &
            timestep_runtime,&
            dtm_new
    integer :: timestep=1, start_timestep=1,  column_global_x, column_global_y, column_local_x, column_local_y,  termination_reason
    logical :: first_timestep_column, last_timestep_column, halo_column, first_nonhalo_timestep_column, update_dtm
    double precision :: model_start_wtime
  end type model_state_type
end module state_mod
