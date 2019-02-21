!writes grid information
!just dumps x, y, z and the tag
module writenetcdf_mod
  use monc_component_mod, only: component_descriptor_type
  implicit none
  public writenetcdf_get_descriptor
  
contains

  type(component_descriptor_type) function writenetcdf_get_descriptor()
    writenetcdf_get_descriptor%name="writenetcdf"
    writenetcdf_get_descriptor%version=0.1
  end function writenetcdf_get_descriptor

end module
