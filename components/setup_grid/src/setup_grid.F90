!Component to set up prognostic grids for Moist Parcel In Cell code
!This should be called before the parcelsetup component

module setup_grid_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_ERROR, log_log


  implicit none

  integer :: nx, ny, nz
  real(kind=DEFAULT_PRECISION) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(kind=DEFAULT_PRECISION) :: dx, dy, dz

contains

  type(component_descriptor_type) function setup_grid_get_descriptor()
    setup_grid_get_descriptor%name="setup_grid"
    setup_grid_get_descriptor%version=0.1
    setup_grid_get_descriptor%initialisation=>initialisation_callback
  end function setup_grid_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    call read_configuration(state)

    call create_grid(state,state%global_grid)
    call decompose_grid(state)
    call allocate_prognostics(state)

  end subroutine

  subroutine read_configuration(state)
    type(model_state_type), intent(inout) :: state

    nx=options_get_integer(state%options_database, "nx")
    ny=options_get_integer(state%options_database, "ny")
    nz=options_get_integer(state%options_database, "nz")

    xmin=options_get_real(state%options_database,"xmin")
    ymin=options_get_real(state%options_database,"ymin")
    zmin=options_get_real(state%options_database,"zmin")
    xmax=options_get_real(state%options_database,"xmax")
    ymax=options_get_real(state%options_database,"ymax")
    zmax=options_get_real(state%options_database,"zmax")

  end subroutine

  subroutine create_grid(state, global_grid)
    type(model_state_type), intent(inout) :: state
    type(global_grid_type), intent(inout) :: global_grid

    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dz=(zmax-zmin)/(nz-1)

    global_grid%resolution(X_INDEX) = dx
    global_grid%resolution(Y_INDEX) = dy
    global_grid%resolution(Z_INDEX) = dz

    global_grid%bottom(X_INDEX) = xmin
    global_grid%bottom(Y_INDEX) = ymin
    global_grid%bottom(Z_INDEX) = zmin

    global_grid%top(X_INDEX) = xmax
    global_grid%top(Y_INDEX) = ymax
    global_grid%top(Z_INDEX) = zmax

    global_grid%size(X_INDEX) = nx
    global_grid%size(Y_INDEX) = ny
    global_grid%size(Z_INDEX) = nz

    global_grid%active(:) = .true.

    global_grid%dimensions = 3

  end subroutine

  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  subroutine allocate_prognostics(state)
    type(model_state_type), intent(inout) :: state
    integer :: nnx, nny, nnz

    !get number of cells in local grid
    nnx=state%local_grid%size(X_INDEX) + 2*state%local_grid%halo_size(X_INDEX)
    nny=state%local_grid%size(Y_INDEX) + 2*state%local_grid%halo_size(Y_INDEX)
    nnz=state%local_grid%size(Z_INDEX) + 2*state%local_grid%halo_size(Z_INDEX)


    print*, "Local grid sizes (+halos)=", nnx, nny, nnz

    !allocate local grids

    allocate(state%u%data(nnz,nny,nnx))
    allocate(state%v%data(nnz,nny,nnx))
    allocate(state%w%data(nnz,nny,nnx))

    allocate(state%p%data(nnz,nny,nnx))
    allocate(state%q%data(nnz,nny,nnx))
    allocate(state%r%data(nnz,nny,nnx))

    allocate(state%b%data(nnz,nny,nnx))
    allocate(state%hg%data(nnz,nny,nnx))
    allocate(state%hgliq%data(nnz,nny,nnx))

    allocate(state%vol%data(nnz,nny,nnx))

  end subroutine

end module
