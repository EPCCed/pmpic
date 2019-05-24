!Component for adding large-scale forcings (cooling/drying etc)
module large_scale_forcing_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use state_mod, only: model_state_type
  use parcel_interpolation_mod, only:  x_coords, y_coords, z_coords
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_real, options_get_logical
  use MPI
  use timer_mod

  implicit none
  real(kind=DEFAULT_PRECISION) :: pi=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION) 
  real(kind=DEFAULT_PRECISION) :: tchar = 142.8571428571_DEFAULT_PRECISION
  real(kind=DEFAULT_PRECISION) :: fdtheta = 3.0_DEFAULT_PRECISION
  real(kind=DEFAULT_PRECISION) :: dryrate
  real(kind=DEFAULT_PRECISION) :: coolrate
  integer :: ierr
  integer :: handle
  integer :: iteration=0
  integer :: nx, ny, nz
  integer :: myrank
  logical :: l_cooling

contains

  type(component_descriptor_type) function large_scale_forcing_get_descriptor()
    large_scale_forcing_get_descriptor%name="large_scale_forcing"
    large_scale_forcing_get_descriptor%version=0.1
    large_scale_forcing_get_descriptor%initialisation=>initialisation_callback
    large_scale_forcing_get_descriptor%timestep=>timestep_callback
    large_scale_forcing_get_descriptor%finalisation=>finalisation_callback
  end function large_scale_forcing_get_descriptor
  
  pure real(kind=DEFAULT_PRECISION) function coolprof(z)
    real(kind=DEFAULT_PRECISION), intent(in) :: z
    coolprof = -coolrate*tchar/(fdtheta*86400.0_DEFAULT_PRECISION) 
  end function coolprof
  
  pure real(kind=DEFAULT_PRECISION) function dryprof(z,dt)
    real(kind=DEFAULT_PRECISION), intent(in) :: z, dt
    dryprof = exp(-dryrate*dt)
  end function dryprof

  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    
    myrank=state%parallel%my_rank
    

    if (myrank .eq. 0) print *, "In Large-scale Initialisation"
    
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)
    
    l_cooling = options_get_logical(state%options_database, "cooling_enabled")   
    coolrate = options_get_real(state%options_database, "coolrate")
    dryrate = options_get_real(state%options_database, "dryrate")
        
    call register_routine_for_timing("large_scale_forcing",handle,state)

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: dt, dz, x, y, z, x_c, y_c 
    integer :: start_x, start_y, end_x, end_y
    integer(kind=PARCEL_INTEGER) :: nparcels, n
    integer :: i, j

    !only do this if we're on the final step of the rk integrator
    if (mod(iteration,state%rksteps) == 0) then
      call timer_start(handle)

      nparcels=state%parcels%numparcels_local

      dt = state%dtm
      dz = state%global_grid%resolution(1)
 
      if (l_cooling) then
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, z)
        !$OMP DO 
        do n=1,nparcels
          z = state%parcels%z(n)
        
          state%parcels%b(n) = state%parcels%b(n)+coolprof(z)*dt
          state%parcels%h(n) = state%parcels%h(n)*dryprof(z,dt)
      
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
      endif
      call timer_stop(handle)
  endif
  
  iteration = iteration + 1
  
  end subroutine

subroutine finalisation_callback(state)
  type(model_state_type), intent(inout), target :: state
  
  if (myrank .eq. 0) print *, "Large-scale forcing finalisation "
  
  if (myrank .eq. 0) print *, "Done!"
  
end subroutine finalisation_callback

end module
