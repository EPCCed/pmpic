!very basic parcel writing routine for basic debugging
!just dumps x, y, z and the tag
module writeparcels_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, LONG_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer,options_get_logical
  use parcel_interpolation_mod, only : x_coords, y_coords, z_coords

  implicit none

  integer :: num, ppersteps, gpersteps, pwritten, gwritten

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

    if (state%parallel%my_rank .eq. 0) print *, "Writer initialisation"

    ppersteps=options_get_integer(state%options_database,"parcel_dump_frequency")
    gpersteps=options_get_integer(state%options_database,"grid_dump_frequency")

    if (options_get_logical(state%options_database,"restart")) then
      print *, "RESTART WRITEFILES"
      num=options_get_integer(state%options_database,"restart_num")
      !stop
    else
      num=0
    endif
    pwritten=0
    gwritten=0

  end subroutine


  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    character (len=20) :: filename
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels

    if(mod(num,ppersteps) .eq. 0) then
      call write_parcels_to_file(state)
      pwritten=pwritten+1
    endif

    if (mod(num,gpersteps) .eq. 0 .and. num .ne. 0) then
      call write_grids_to_file(state)
      gwritten=gwritten+1
    endif

    num=num+1
  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) then
      print *, "Written", pwritten, "parcel dumps"
      print *, "Written", gwritten, "grid dumps"
    endif


  end subroutine



  subroutine write_parcels_to_file(state)
    type(model_state_type), intent(inout), target :: state
    character (len=20) :: filename
    character (len=23) :: fnamedummy
    integer :: proc
    integer(kind=PARCEL_INTEGER) :: nparcels


    proc=state%parallel%my_rank
    nparcels=state%parcels%numparcels_local

    write(filename,"(A8,i3.3,A1,I4.4,A4)") "parcels_", proc,"_", num, ".dat"
    write(fnamedummy,"(A,I4.4,A4)") "parcels_[rank]_", num, ".dat"

    if (proc .eq. 0) print *, "Writing parcels to '",fnamedummy,"'"

    !file contains:
    ! t
    ! xmin, xmax, ymin, ymax, zmin, zmax
    ! n_parcels
    ! data
    open(unit=10,file=filename,form="unformatted",access="stream")
    write(10) state%time

    write(10) x_coords(state%local_grid%local_domain_start_index(3))
    write(10) x_coords(state%local_grid%local_domain_end_index(3)+1)
    write(10) y_coords(state%local_grid%local_domain_start_index(2))
    write(10) y_coords(state%local_grid%local_domain_end_index(2)+1)
    write(10) z_coords(state%local_grid%local_domain_start_index(1))
    write(10) z_coords(state%local_grid%local_domain_end_index(1))

    write(10) nparcels

    write(10) state%parcels%x(1:nparcels)
    write(10) state%parcels%y(1:nparcels)
    write(10) state%parcels%z(1:nparcels)

    write(10) state%parcels%p(1:nparcels)
    write(10) state%parcels%q(1:nparcels)
    write(10) state%parcels%r(1:nparcels)

    write(10) state%parcels%dxdt(1:nparcels)
    write(10) state%parcels%dydt(1:nparcels)
    write(10) state%parcels%dzdt(1:nparcels)

    write(10) state%parcels%dpdt(1:nparcels)
    write(10) state%parcels%dqdt(1:nparcels)
    write(10) state%parcels%drdt(1:nparcels)

    write(10) state%parcels%h(1:nparcels)
    write(10) state%parcels%b(1:nparcels)
    write(10) state%parcels%vol(1:nparcels)

    write(10) state%parcels%stretch(1:nparcels)

    write(10) state%parcels%tag(1:nparcels)
    write(10) state%parcels%qvalues(:,1:nparcels)

    close(10)

  end subroutine

  subroutine write_grids_to_file(state)
    type(model_state_type), intent(in) :: state
    character (len=18) :: filename
    character (len=21) :: fnamedummy
    integer :: proc
    integer :: i1, i2, j1, j2, k1, k2
    real(kind=DEFAULT_PRECISION) :: x1, x2, y1, y2, z1, z2

    i1=state%local_grid%local_domain_start_index(1)
    i2=state%local_grid%local_domain_end_index(1)

    j1=state%local_grid%local_domain_start_index(2)
    j2=state%local_grid%local_domain_end_index(2)

    k1=state%local_grid%local_domain_start_index(3)
    k2=state%local_grid%local_domain_end_index(3)

    x1=x_coords(k1)
    x2=x_coords(k2+1)

    y1=y_coords(j1)
    y2=y_coords(j2+1)

    z1=z_coords(i1)
    z2=z_coords(i2)

    proc=state%parallel%my_rank

    write(filename,"(A,i3.3,A1,I4.4,A4)") "grids_", proc,"_", num, ".dat"
    write(fnamedummy,"(A,I4.4,A4)") "grids_[rank]_", num, ".dat"

    if (proc .eq. 0) print *, "Writing grids to '",fnamedummy,"'"

    ! print *, "nx, ny, nz=", state%local_grid%size(3), state%local_grid%size(2), state%local_grid%size(1)
    ! print *, "xmin, xmax=", x_coords(state%local_grid%local_domain_start_index(3)), &
    ! x_coords(state%local_grid%local_domain_end_index(3))
    ! print *, "ymin, ymax=", y_coords(state%local_grid%local_domain_start_index(2)), &
    ! y_coords(state%local_grid%local_domain_end_index(2))
    ! print *, "zmin, zmax=", z_coords(state%local_grid%local_domain_start_index(1)), &
    ! z_coords(state%local_grid%local_domain_end_index(1))


    !file contains:
    ! t
    ! xmin, xmax, ymin, ymax, zmin, zmax
    ! nx, ny, nz
    ! data
    open(unit=10,file=filename,form="unformatted",access="stream")
    write(10) state%time
    write(10) x_coords(state%local_grid%local_domain_start_index(3))
    write(10) x_coords(state%local_grid%local_domain_end_index(3))
    write(10) y_coords(state%local_grid%local_domain_start_index(2))
    write(10) y_coords(state%local_grid%local_domain_end_index(2))
    write(10) z_coords(state%local_grid%local_domain_start_index(1))
    write(10) z_coords(state%local_grid%local_domain_end_index(1))
    write(10) state%local_grid%size(3), state%local_grid%size(2), state%local_grid%size(1)
    write(10) state%u%data(i1:i2, j1:j2, k1:k2)
    write(10) state%v%data(i1:i2, j1:j2, k1:k2)
    write(10) state%w%data(i1:i2, j1:j2, k1:k2)
    write(10) state%p%data(i1:i2, j1:j2, k1:k2)
    write(10) state%q%data(i1:i2, j1:j2, k1:k2)
    write(10) state%r%data(i1:i2, j1:j2, k1:k2)
    write(10) state%b%data(i1:i2, j1:j2, k1:k2)
    write(10) state%hg%data(i1:i2, j1:j2, k1:k2)
    write(10) state%hgliq%data(i1:i2, j1:j2, k1:k2)
    write(10) state%vol%data(i1:i2, j1:j2, k1:k2)
    close(10)

  end subroutine


end module
