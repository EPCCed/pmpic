!A component for generating surface fluxes from an analytic profile
module surface_fluxes_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_real
  use MPI
  use timer_mod

  implicit none
  real(kind=DEFAULT_PRECISION) :: pi=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION) 
  integer :: ierr
  integer :: handle
  integer :: iteration=0

contains

  type(component_descriptor_type) function surface_fluxes_get_descriptor()
    surface_fluxes_get_descriptor%name="surface_fluxes"
    surface_fluxes_get_descriptor%version=0.1
    surface_fluxes_get_descriptor%initialisation=>initialisation_callback
    surface_fluxes_get_descriptor%timestep=>timestep_callback
    !surface_fluxes_get_descriptor%finalisation=>finalisation_callback
  end function surface_fluxes_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) print *, "In Surface Flux Initialisation"
 
    call register_routine_for_timing("surface_fluxes",handle,state)

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: dt, dx, dy, dz, x, y, z, x_c, y_c 
    real(kind=DEFAULT_PRECISION) :: x_shift, y_shift, fluxfac, AmpB, AmpH
    real(kind=DEFAULT_PRECISION) :: eps, width, del, y_transform, prefac
    integer(kind=PARCEL_INTEGER) :: nparcels, n
    

    !only do this if we're on the final step of the rk integrator
    if (mod(iteration,state%rksteps) == 0) then
      call timer_start(handle)

      nparcels=state%parcels%numparcels_local

      dt = state%dtm
      AmpB = 0.01_DEFAULT_PRECISION
      AmpH = 0.00005_DEFAULT_PRECISION
      eps = 0.2_DEFAULT_PRECISION
      width = 1.0_DEFAULT_PRECISION
      del = 0.2_DEFAULT_PRECISION
      
      dx = state%global_grid%resolution(3)
      dy = state%global_grid%resolution(2)
      dz = state%global_grid%resolution(1)

      x_c = (state%global_grid%top(3)+dx-state%global_grid%bottom(3))/2 +state%global_grid%bottom(3)
      y_c = (state%global_grid%top(2)+dy-state%global_grid%bottom(2))/2 +state%global_grid%bottom(2)
      
      fluxfac = 2.0_DEFAULT_PRECISION*dt/(dz*dz)
    
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, x_shift, y_shift, z, y_transform, prefac)
      !$OMP DO
      do n=1,nparcels
      
        x_shift = state%parcels%x(n) - x_c
        y_shift = state%parcels%y(n) - y_c
        z = state%parcels%z(n)
      
        if (z .lt. dz) then
      
          y_transform = y_shift - del*sin(x_shift)
        
          prefac = (1.0_DEFAULT_PRECISION + eps*cos(x_shift)) &
		       *(exp(-(y_transform/width)**2.0_DEFAULT_PRECISION)) &
		       *(dz - z)*fluxfac
                
	      state%parcels%b(n) = state%parcels%b(n) + AmpB*prefac
	      state%parcels%h(n) = state%parcels%h(n) + AmpH*prefac
	    
	  
	    endif		

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      

      call timer_stop(handle)
  endif
  
  iteration = iteration + 1
  
  end subroutine




end module
