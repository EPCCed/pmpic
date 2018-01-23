!> Dry boundary layer test case, which represents test case 1 in the LEM
module basic_plume_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only :model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif


  public basic_plume_get_descriptor
contains

  type(component_descriptor_type) function basic_plume_get_descriptor()
    basic_plume_get_descriptor%name="basic_plume"
    basic_plume_get_descriptor%version=0.1
    basic_plume_get_descriptor%initialisation=>initialisation_callback
  end function basic_plume_get_descriptor

  !> Sets up the field values for this test case
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

   ! Initialisation code goes here


  end subroutine initialisation_callback
end module basic_plume_mod
