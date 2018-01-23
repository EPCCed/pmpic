!> General purpose debugger. By changing the priority and other logic we can plug it in
!! whereever we want in the run to dump out information
module debugger_mod
  use grids_mod
  use monc_component_mod
  use state_mod
  use logging_mod
  use prognostics_mod
  use conversions_mod
  implicit none

#ifndef TEST_MODE
  private
#endif

public debugger_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function debugger_get_descriptor()
    debugger_get_descriptor%name="debugger"
    debugger_get_descriptor%version=0.1
    debugger_get_descriptor%initialisation=>init_callback
    debugger_get_descriptor%timestep=>timestep_callback
  end function debugger_get_descriptor

  !> Called on MONC initialisation
  !! @param current_state The current model stat
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call log_log(LOG_WARN, "Debugger is active - disable this for production runs")
  end subroutine init_callback  

  !> Produces debugging information on each timestep
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
   
  end subroutine timestep_callback  
end module debugger_mod
