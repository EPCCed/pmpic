!writes grid information
!just dumps x, y, z and the tag
module writegrids_mod
  use monc_component_mod, only: component_descriptor_type
  implicit none
  public writegrids_get_descriptor
  
contains

  type(component_descriptor_type) function writegrids_get_descriptor()
    writegrids_get_descriptor%name="writegrids"
    writegrids_get_descriptor%version=0.1
  end function writegrids_get_descriptor

end module
