!A component for generating surface fluxes from an analytic profile
module surface_fluxes_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use state_mod, only: model_state_type
  use parcel_interpolation_mod, only:  x_coords, y_coords, z_coords
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_real
  use MPI
  use timer_mod

  implicit none
  real(kind=DEFAULT_PRECISION) :: pi=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION) 
  real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:) :: bflux, hflux 
  integer :: ierr
  integer :: handle
  integer :: iteration=0
  integer :: nx, ny
  integer :: myrank

contains

  type(component_descriptor_type) function surface_fluxes_get_descriptor()
    surface_fluxes_get_descriptor%name="surface_fluxes"
    surface_fluxes_get_descriptor%version=0.1
    surface_fluxes_get_descriptor%initialisation=>initialisation_callback
    surface_fluxes_get_descriptor%timestep=>timestep_callback
    surface_fluxes_get_descriptor%finalisation=>finalisation_callback
  end function surface_fluxes_get_descriptor
  
  pure real(kind=DEFAULT_PRECISION) function fluxprof(y, x, eps, width, del, y_c, x_c)
    real(kind=DEFAULT_PRECISION), intent(in) :: y, x, eps, width, del, y_c, x_c
    real(kind=DEFAULT_PRECISION) :: y_shift, x_shift, y_transform 

    x_shift = x - x_c
    y_shift = y - y_c
    
    y_transform = y_shift - del*sin(x_shift)
    fluxprof = (1.0_DEFAULT_PRECISION + eps*cos(x_shift)) &
		       *(exp(-(y_transform/width)**2.0_DEFAULT_PRECISION))
  end function fluxprof

  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    
    myrank=state%parallel%my_rank

    if (myrank .eq. 0) print *, "In Surface Flux Initialisation"
    
    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3) 
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    
    allocate(bflux(ny,nx), hflux(ny,nx))
     
    call register_routine_for_timing("surface_fluxes",handle,state)

  end subroutine

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    real(kind=DEFAULT_PRECISION) :: dt, dx, dy, dz, x, y, z, x_c, y_c 
    real(kind=DEFAULT_PRECISION) :: x_shift, y_shift, fluxfac, AmpB, AmpH
    real(kind=DEFAULT_PRECISION) :: eps1, width1, del1, prefac
    real(kind=DEFAULT_PRECISION) :: bfluxsum, bparsum, hfluxsum, hparsum, bfcorr, hfcorr
    integer :: start_x, start_y, end_x, end_y
    integer(kind=PARCEL_INTEGER) :: nparcels, n
    integer :: i, j

    

    !only do this if we're on the final step of the rk integrator
    if (mod(iteration,state%rksteps) == 0) then
      call timer_start(handle)

      nparcels=state%parcels%numparcels_local

      dt = state%dtm
      AmpB = 0.05_DEFAULT_PRECISION
      AmpH = 0.005_DEFAULT_PRECISION
      eps1 = 0.2_DEFAULT_PRECISION
      width1 = 1.0_DEFAULT_PRECISION
      del1 = 0.2_DEFAULT_PRECISION
      
      dx = state%global_grid%resolution(3)
      dy = state%global_grid%resolution(2)
      dz = state%global_grid%resolution(1)
 
      x_c = (state%global_grid%top(3)+dx-state%global_grid%bottom(3))/2 +state%global_grid%bottom(3)
      y_c = (state%global_grid%top(2)+dy-state%global_grid%bottom(2))/2 +state%global_grid%bottom(2)
      
      fluxfac = 2.0_DEFAULT_PRECISION/(dz*dz)
    
      bparsum = 0.0_DEFAULT_PRECISION
      hparsum = 0.0_DEFAULT_PRECISION
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, x, y, z, prefac)
      !$OMP DO REDUCTION(+:bparsum, hparsum)

      do n=1,nparcels
        x = state%parcels%x(n)
        y = state%parcels%y(n)
        z = state%parcels%z(n)
      
        if (z .lt. dz) then
                       
          prefac = fluxprof(y, x, eps1, width1, del1, y_c, x_c)
		  prefac = prefac*(dz - z)*fluxfac
	      
	      bparsum = bparsum + AmpB*prefac*state%parcels%vol(n)
	      hparsum = hparsum + AmpH*prefac*state%parcels%vol(n)
	    endif		

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      call MPI_Allreduce(MPI_IN_PLACE,&
             bparsum,&
             1,&
             PRECISION_TYPE,&
             MPI_SUM,&
             state%parallel%monc_communicator,&
             ierr) 
      
      call MPI_Allreduce(MPI_IN_PLACE,&
             hparsum,&
             1,&
             PRECISION_TYPE,&
             MPI_SUM,&
             state%parallel%monc_communicator,&
             ierr) 
      
      do i=1,nx
        x = x_coords(i)
        do j=1,ny
          y = y_coords(j)
          prefac = fluxprof(y, x, eps1, width1, del1, y_c, x_c) 
          bflux(j,i) = AmpB*prefac
          hflux(j,i) = AmpH*prefac
        enddo
      enddo
      
      start_x=state%local_grid%local_domain_start_index(3)
      end_x=state%local_grid%local_domain_end_index(3)
      
      start_y=state%local_grid%local_domain_start_index(2)
      end_y=state%local_grid%local_domain_end_index(2)
      
      bfluxsum = sum(bflux(start_y:end_y, start_x:end_x))*dx*dy
      hfluxsum = sum(hflux(start_y:end_y, start_x:end_x))*dx*dy
      
      call MPI_Allreduce(MPI_IN_PLACE,&
             bfluxsum,&
             1,&
             PRECISION_TYPE,&
             MPI_SUM,&
             state%parallel%monc_communicator,&
             ierr)  
                     
      call MPI_Allreduce(MPI_IN_PLACE,&
             hfluxsum,&
             1,&
             PRECISION_TYPE,&
             MPI_SUM,&
             state%parallel%monc_communicator,&
             ierr)      
             
      bfcorr = bfluxsum/bparsum
      hfcorr = hfluxsum/hparsum
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, x, y, z, prefac)
      !$OMP DO

      do n=1,nparcels
        x = state%parcels%x(n)
        y = state%parcels%y(n)
        z = state%parcels%z(n)
      
        if (z .lt. dz) then
                       
          prefac = fluxprof(y, x, eps1, width1, del1, y_c, x_c)
		  prefac = prefac*(dz - z)*fluxfac
                
	      state%parcels%b(n) = state%parcels%b(n) + AmpB*prefac*bfcorr*dt
	      state%parcels%h(n) = state%parcels%h(n) + AmpH*prefac*hfcorr*dt
	      
	    endif		

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      write(*,*) "bfcorr = ", bfcorr
      write(*,*) "hfcorr = ", hfcorr

      call timer_stop(handle)
  endif
  
  
  
  iteration = iteration + 1
  
  end subroutine

subroutine finalisation_callback(state)
  type(model_state_type), intent(inout), target :: state
  
  if (myrank .eq. 0) print *, "Deallocating gridded surface flux profiles"
  deallocate(bflux)
  deallocate(hflux)
  
  if (myrank .eq. 0) print *, "Done!"
  
end subroutine finalisation_callback

end module
