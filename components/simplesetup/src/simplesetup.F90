module simplesetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use logging_mod, only : LOG_ERROR, log_log
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_integer_array, options_get_real_array
  use q_indices_mod, only: get_q_index, standard_q_names
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: x_size, y_size, z_size
  real(kind=DEFAULT_PRECISION) :: zztop, dxx, dyy
  logical :: enable_theta=.false.
  public simplesetup_get_descriptor
contains

  type(component_descriptor_type) function simplesetup_get_descriptor()
    simplesetup_get_descriptor%name="simplesetup"
    simplesetup_get_descriptor%version=0.1
    simplesetup_get_descriptor%initialisation=>initialisation_callback
  end function simplesetup_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    call read_configuration(current_state)
    if (.not. current_state%initialised) then
      current_state%dtm=options_get_real(current_state%options_database, "dtm")
      current_state%dtm_new=current_state%dtm
      call create_grid(current_state, current_state%global_grid)
      call decompose_grid(current_state)
      current_state%initialised=.true.

    end if
  end subroutine initialisation_callback

  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  subroutine create_grid(current_state, specific_grid)
    type(model_state_type), intent(inout) :: current_state
    type(global_grid_type), intent(inout) :: specific_grid

    integer, parameter :: KGD_SIZE=20
    integer :: number_kgd, i, kgd(KGD_SIZE)
    real(kind=DEFAULT_PRECISION) :: hgd(KGD_SIZE)

    kgd=-1

    call options_get_integer_array(current_state%options_database, "kgd", kgd)
    call options_get_real_array(current_state%options_database, "hgd", hgd)  

    if (kgd(1)==1)then
      if (hgd(1)/=0.0_DEFAULT_PRECISION)then
        call log_log(LOG_ERROR, "Lowest level is assumed to lie at the surface, check hgd(1)")
      else
        kgd(1:KGD_SIZE-1) = kgd(2:)
        hgd(1:KGD_SIZE-1) = hgd(2:)
      end if
    end if

    do i=1,size(kgd)
      if (kgd(i) == -1) exit      
    end do
    number_kgd=i-1

    if (number_kgd .gt. 0) then
      allocate(current_state%global_grid%configuration%vertical%kgd(number_kgd), &
           current_state%global_grid%configuration%vertical%hgd(number_kgd))
      current_state%global_grid%configuration%vertical%kgd=kgd(1:number_kgd)
      current_state%global_grid%configuration%vertical%hgd=hgd(1:number_kgd)
    end if

    specific_grid%bottom(Z_INDEX) = 0
    specific_grid%bottom(Y_INDEX) = 0
    specific_grid%bottom(X_INDEX) = 0

    specific_grid%top(Z_INDEX) = zztop
    specific_grid%top(Y_INDEX) = dyy * y_size
    specific_grid%top(X_INDEX) = dxx * x_size

    specific_grid%resolution(Z_INDEX) = zztop / z_size
    specific_grid%resolution(Y_INDEX) = dyy
    specific_grid%resolution(X_INDEX) = dxx

    specific_grid%size(Z_INDEX) = z_size
    specific_grid%size(Y_INDEX) = y_size
    specific_grid%size(X_INDEX) = x_size

    specific_grid%active(Z_INDEX) = .true.
    specific_grid%active(Y_INDEX) = .true.
    specific_grid%active(X_INDEX) = .true.

    specific_grid%dimensions = 3
  end subroutine create_grid

  subroutine read_configuration(current_state)
    type(model_state_type), intent(inout), target :: current_state

    x_size=options_get_integer(current_state%options_database, "x_size")
    y_size=options_get_integer(current_state%options_database, "y_size")
    z_size=options_get_integer(current_state%options_database, "z_size")
    dxx=options_get_real(current_state%options_database, "dxx")
    dyy=options_get_real(current_state%options_database, "dyy")
    zztop=options_get_real(current_state%options_database, "zztop")

  end subroutine read_configuration
end module simplesetup_mod
