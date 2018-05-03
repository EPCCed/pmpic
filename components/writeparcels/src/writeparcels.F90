!very basic parcel writing routine for basic debugging
!just dumps x, y, z and the tag
module writeparcels_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer

  implicit none

  integer :: num, persteps, written

contains

  type(component_descriptor_type) function writeparcels_get_descriptor()
    writeparcels_get_descriptor%name="writeparcels"
    writeparcels_get_descriptor%version=0.1
    writeparcels_get_descriptor%initialisation=>initialisation_callback
    writeparcels_get_descriptor%timestep=>timestep_callback
    writeparcels_get_descriptor%finalisation=>finalisation_callback
  end function writeparcels_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print *, "Writer initialisation"

    persteps=options_get_integer(state%options_database,"dump_frequency")

    num=0
    written=0

  end subroutine


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    character (len=20) :: filename
    integer :: proc, nparcels

    if(mod(num,persteps) .eq. 0) then

      proc=state%parallel%my_rank
      nparcels=state%parcels%numparcels_local

      write(filename,"(A8,i3.3,A1,I4.4,A4)") "parcels_", proc,"_", num, ".dat"

      print *, "Writing parcels to '",filename,"'"

      open(unit=10,file=filename,form="unformatted")
      write(10) state%time
      write(10) nparcels
      write(10) state%parcels%x(1:nparcels)
      write(10) state%parcels%y(1:nparcels)
      write(10) state%parcels%z(1:nparcels)
      write(10) state%parcels%tag(1:nparcels)
      close(10)

      written=written+1
    endif

    num=num+1
  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    print *, "Written", written, "files"

  end subroutine

end module
