!> Scientific constant values used throughout simulations. Each has a default value and this can be overridden by the
!! configuration supplied by the user
module science_constants_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use optionsdatabase_mod, only : options_get_string, options_get_real
  use state_mod, only : model_state_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION) :: smallp=1.0e-14, von_karman_constant, z0, z0th, alphah, betam, betah, &
       gammam, gammah, pi, surface_vapour_mixing_ratio, cp , rlvap, rlvap_over_cp, r, r_over_cp, G,&
       convective_limit, thref0, pref0, ratio_mol_wts,rlargep,q0,l_condense,lapse_ref,surf_temp_base,&
       surf_pres_base,lapse_exp,p_ref_prefactor,cp_over_rlvap,inv_exn_factor

  real(kind=DEFAULT_PRECISION) :: seconds_in_a_day=86400.0
  public smallp, von_karman_constant, z0, z0th, alphah, betam, betah, gammam, gammah, pi, cp, &
       rlvap, rlvap_over_cp, r, r_over_cp, G, convective_limit, ratio_mol_wts, rlargep, &
       initialise_science_constants, &
       seconds_in_a_day,thref0, pref0, q0,l_condense,lapse_ref,surf_temp_base,surf_pres_base,&
       lapse_exp,p_ref_prefactor,cp_over_rlvap,inv_exn_factor
contains

  !> Initialises the scientific constants to read in any values that are overridden in the configuration
  !! @param current_state The current model state
  subroutine initialise_science_constants(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    !von_karman_constant=options_get_real(current_state%options_database, "von_karman_constant")
    !z0=options_get_real(current_state%options_database, "z0")
    !z0th=options_get_real(current_state%options_database, "z0th")
    !alphah=options_get_real(current_state%options_database, "alphah")
    !betam=options_get_real(current_state%options_database, "betam")
    !betah=options_get_real(current_state%options_database, "betah")
    !gammam=options_get_real(current_state%options_database, "gammam")
    !gammah=options_get_real(current_state%options_database, "gammah")
    !pi=options_get_real(current_state%options_database, "pi")
    cp=options_get_real(current_state%options_database, "cp")
    rlvap=options_get_real(current_state%options_database, "rlvap")
    r=options_get_real(current_state%options_database, "r")
    G=options_get_real(current_state%options_database, "G")
    convective_limit=options_get_real(current_state%options_database, "convective_limit")
    thref0=options_get_real(current_state%options_database, "thref0")
    pref0=options_get_real(current_state%options_database, "pref0")
    q0=options_get_real(current_state%options_database, "q0")
    l_condense=options_get_real(current_state%options_database, "l_condense")
    lapse_ref=options_get_real(current_state%options_database, "lapse_ref")
    surf_temp_base=options_get_real(current_state%options_database, "surf_temp_base")
    surf_pres_base=options_get_real(current_state%options_database, "surf_pres_base")
    ratio_mol_wts=options_get_real(current_state%options_database, "ratio_mol_wts")
    !rlargep=options_get_real(current_state%options_database, "rlargep")
    lapse_exp=G/(r*lapse_ref)
    p_ref_prefactor=surf_pres_base/(surf_temp_base**lapse_exp)
    rlvap_over_cp=rlvap/cp
    cp_over_rlvap=cp/rlvap
    r_over_cp=r/cp
    inv_exn_factor=(1.0/pref0)**r_over_cp
  end subroutine initialise_science_constants  
end module science_constants_mod
