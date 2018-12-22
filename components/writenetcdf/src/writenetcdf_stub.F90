!writes grid information
!just dumps x, y, z and the tag
module writenetdf_mod
  use monc_component_mod, only: component_descriptor_type
  implicit none
  public writenetdf_get_descriptor
  
contains

  type(component_descriptor_type) function writenetdf_get_descriptor()
    writenetdf_get_descriptor%name="writenetdf"
    writenetdf_get_descriptor%version=0.1
  end function writenetdf_get_descriptor

end module
