module fileunits_mod

implicit none

integer, parameter :: minunit=20
integer, parameter :: maxunit=200

contains

  !returns the number of a free (i.e. unopened) unit number
  integer function get_free_file_unit()
    integer :: i
    logical :: isopen
    get_free_file_unit=-1
    do i=minunit,maxunit
      inquire(unit=i, opened=isopen)
      if (.not. isopen) then
        get_free_file_unit=i
        exit
      endif
    enddo

    if (get_free_file_unit .eq. -1) error stop "No free unit number available"
  end function

  ! checks whether a file exists or not
  logical function file_exists(fname)
    character (len=*), intent(in) :: fname

    inquire(file=fname,exist=file_exists)

  end function


end module fileunits_mod
