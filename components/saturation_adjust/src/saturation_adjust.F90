! Saturation adjustment
! Stripped module

module saturation_adjust_mod
  use datadefn_mod, only : DEFAULT_PRECISION,PARCEL_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use science_constants_mod, only : thref0,q0,l_condense,surf_temp_base,lapse_exp,lapse_ref,p_ref_prefactor,&
  cp_over_rlvap,inv_exn_factor,r_over_cp,rlvap_over_cp
  use saturation_mod, only: qsaturation
  use q_indices_mod, only: get_q_index,standard_q_names
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop

  implicit none

  integer :: ierr
  integer :: handle

contains

  type(component_descriptor_type) function saturation_adjust_get_descriptor()
    saturation_adjust_get_descriptor%name="saturation_adjust"
    saturation_adjust_get_descriptor%version=0.1
    saturation_adjust_get_descriptor%initialisation=>initialisation_callback
    saturation_adjust_get_descriptor%timestep=>timestep_callback
  end function saturation_adjust_get_descriptor

  ! Reference pressure, given surface temperature, surface pressure and a constant lapse rate
  pure real(kind=DEFAULT_PRECISION) function p_ref(z)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in)  :: z
    p_ref=p_ref_prefactor*((surf_temp_base-lapse_ref*z)**lapse_exp)
    return
  end function p_ref

  ! Saturation adjustment, check if 
  real(kind=DEFAULT_PRECISION) FUNCTION error_temp(qv,qc,delta_temp,theta,exn,p)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: qv
    real(kind=DEFAULT_PRECISION), intent(in) :: qc
    real(kind=DEFAULT_PRECISION), intent(in) :: delta_temp
    real(kind=DEFAULT_PRECISION), intent(in) :: theta
    real(kind=DEFAULT_PRECISION), intent(in) :: exn
    real(kind=DEFAULT_PRECISION), intent(in) :: p
    real(kind=DEFAULT_PRECISION) :: temp_guess
    real(kind=DEFAULT_PRECISION) :: qc_guess
    real(kind=DEFAULT_PRECISION) :: qv_guess
    real(kind=DEFAULT_PRECISION) :: delta_qv
    qc_guess=qc+delta_temp*cp_over_rlvap
    if(qc_guess>0.0_DEFAULT_PRECISION) then
      ! predicted regime saturated: prevent qv unequal to qsat
      qv_guess=qv-delta_temp*cp_over_rlvap
      temp_guess=theta*exn+delta_temp
      error_temp=-(qv_guess-qsaturation(temp_guess,0.01_DEFAULT_PRECISION*p))*rlvap_over_cp
    else
      ! predicted regime unsaturated: prevent negative qc
      error_temp= qc_guess*rlvap_over_cp
    end if
  end function error_temp

  subroutine initialisation_callback(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state
    integer:: iqv,iqc
    iqv=get_q_index(standard_q_names%VAPOUR, 'saturation_adjust')
    iqc=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'saturation_adjust')
    
    call register_routine_for_timing("saturation_adjust",handle,current_state)

  end subroutine initialisation_callback
  
  subroutine ridder_temp(qv,qc,delta_temp_a_in,delta_temp_b_in,error_temp_a_in,error_temp_b_in,theta,exn,p,l_fail,delta_temp)
    real(kind=DEFAULT_PRECISION), intent(in) :: qv
    real(kind=DEFAULT_PRECISION), intent(in) :: qc
    real(kind=DEFAULT_PRECISION), intent(in) :: theta
    real(kind=DEFAULT_PRECISION), intent(in) :: exn
    real(kind=DEFAULT_PRECISION), intent(in) :: p
    real(kind=DEFAULT_PRECISION), intent(in) :: delta_temp_a_in
    real(kind=DEFAULT_PRECISION), intent(in) :: delta_temp_b_in
    real(kind=DEFAULT_PRECISION), intent(in) :: error_temp_a_in
    real(kind=DEFAULT_PRECISION), intent(in) :: error_temp_b_in
    real(kind=DEFAULT_PRECISION), intent(out) :: delta_temp
    logical, intent(out) :: l_fail
    real(kind=DEFAULT_PRECISION) :: delta_temp_a
    real(kind=DEFAULT_PRECISION) :: delta_temp_b
    real(kind=DEFAULT_PRECISION) :: delta_temp_c
    real(kind=DEFAULT_PRECISION) :: delta_temp_x
    real(kind=DEFAULT_PRECISION) :: delta_temp_x_old
    real(kind=DEFAULT_PRECISION) :: root
    real(kind=DEFAULT_PRECISION) :: shift
    real(kind=DEFAULT_PRECISION) :: error_temp_a
    real(kind=DEFAULT_PRECISION) :: error_temp_b
    real(kind=DEFAULT_PRECISION) :: error_temp_c
    real(kind=DEFAULT_PRECISION) :: error_temp_x
    real(kind=DEFAULT_PRECISION) :: temp_guess
    real(kind=DEFAULT_PRECISION),parameter :: tolerance=0.0000001_DEFAULT_PRECISION
    integer :: iteration_nr,iqv,iqc
    
    delta_temp_a=delta_temp_a_in
    delta_temp_b=delta_temp_b_in
    error_temp_a=error_temp_a_in
    error_temp_b=error_temp_b_in
    ! Ridder's root finding method    
    do iteration_nr=1,20
      delta_temp_c = 0.5_DEFAULT_PRECISION*(delta_temp_a + delta_temp_b)
      error_temp_c = error_temp(qv,qc,delta_temp_c,theta,exn,p)
      root = sqrt(error_temp_c*error_temp_c - error_temp_a*error_temp_b)
      if(root==0.0_DEFAULT_PRECISION) then
        delta_temp=0.0_DEFAULT_PRECISION
        l_fail=.true.
        return
      end if
      shift = (delta_temp_c - delta_temp_a)*error_temp_c/root
      if(shift==0.0_DEFAULT_PRECISION) then
        delta_temp=0.0_DEFAULT_PRECISION
        l_fail=.true.
        return
      end if
      if ((error_temp_a - error_temp_b) < 0.0_DEFAULT_PRECISION) then
        shift = -shift
        delta_temp_x = delta_temp_c + shift
        error_temp_x = error_temp(qv,qc,delta_temp_x,theta,exn,p)
      end if
      if(iteration_nr > 1) then
        if (abs(delta_temp_x - delta_temp_x_old) < tolerance*max(abs(delta_temp_x),1.0_DEFAULT_PRECISION)) then
          delta_temp=delta_temp_x
          return
        end if
      end if
      delta_temp_x_old = delta_temp_x
      if(error_temp_c*error_temp_x > 0.0_DEFAULT_PRECISION) then
        if(error_temp_a*error_temp_x < 0.0_DEFAULT_PRECISION) then 
          delta_temp_b = delta_temp_x
          error_temp_b = error_temp_x
        else
          delta_temp_a = delta_temp_x
          error_temp_a = error_temp_x
        end if
      else
        delta_temp_a = delta_temp_c 
        delta_temp_b = delta_temp_x 
        error_temp_a = error_temp_c 
        error_temp_b = error_temp_x
      end if
    end do
    
    delta_temp=0.0_DEFAULT_PRECISION
    l_fail=.true.
    return
    
  end subroutine ridder_temp
         
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    ! realistic moist thermodynamics
    ! the goal here is to be very accurate, so that 
    real(kind=DEFAULT_PRECISION) :: qv
    real(kind=DEFAULT_PRECISION) :: qc
    real(kind=DEFAULT_PRECISION) :: b
    real(kind=DEFAULT_PRECISION) :: z
    real(kind=DEFAULT_PRECISION) :: exn
    real(kind=DEFAULT_PRECISION) :: theta
    real(kind=DEFAULT_PRECISION) :: p
    real(kind=DEFAULT_PRECISION) :: delta_q    
    real(kind=DEFAULT_PRECISION) :: delta_temp_a
    real(kind=DEFAULT_PRECISION) :: delta_temp_b
    real(kind=DEFAULT_PRECISION) :: delta_temp
    real(kind=DEFAULT_PRECISION) :: error_temp_a
    real(kind=DEFAULT_PRECISION) :: error_temp_b
    real(kind=DEFAULT_PRECISION) :: temp_guess
    integer(kind=PARCEL_INTEGER):: n
    integer(kind=PARCEL_INTEGER):: n_fail
    real(kind=DEFAULT_PRECISION) :: qv_guess
    real(kind=DEFAULT_PRECISION) :: qc_guess
    real(kind=DEFAULT_PRECISION) :: qv_fail
    real(kind=DEFAULT_PRECISION) :: qc_fail
    real(kind=DEFAULT_PRECISION) :: b_fail
    real(kind=DEFAULT_PRECISION) :: z_fail
    logical :: l_fail,l_any_fail
    integer :: iteration_nr,iqv,iqc
    
    call timer_start(handle)

    iqv=get_q_index(standard_q_names%VAPOUR, 'saturation_adjust')
    iqc=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'saturation_adjust')

    l_any_fail=.false.
    
    do n=1,current_state%parcels%numparcels_local
      l_fail=.false.
      b=current_state%parcels%b(n)
      qv=current_state%parcels%qvalues(iqv,n)
      qc=current_state%parcels%qvalues(iqc,n)
      z=current_state%parcels%z(n)
  
      ! start the saturation adjustment procedure 
      p=p_ref(z)
      exn=inv_exn_factor*p**r_over_cp
      theta=b+thref0
      
      ! First see if we can be lazy
      delta_q=-qc
      qc_guess=qc+delta_q
      qv_guess=qv-delta_q
      temp_guess=theta*exn+delta_q*rlvap_over_cp

      if(qsaturation(temp_guess,0.01_DEFAULT_PRECISION*p)>qv) then
        ! Unsaturated regime when all liquid is evaporated
        qc=qc_guess
        qv=qv_guess
        b=temp_guess/exn-thref0
      else
        ! Not all can be evaporated
        ! Hence we are in a saturated regime
        ! And may just as well perform the whole procedure
        ! Bounds on temperature given by
        ! 1) No change
        ! 2) The first order correction (i.e. No feedbacks)
        
        delta_temp=(qv-qsaturation(theta*exn,0.01*p)-min(qc,0.0_DEFAULT_PRECISION))*rlvap_over_cp
        if(delta_temp>0.0_DEFAULT_PRECISION) then
          delta_temp_a=-0.05_DEFAULT_PRECISION*delta_temp-0.000001_DEFAULT_PRECISION
          delta_temp_b=1.05_DEFAULT_PRECISION*delta_temp+0.000001_DEFAULT_PRECISION        
        else
          delta_temp_a=1.05_DEFAULT_PRECISION*delta_temp-0.000001_DEFAULT_PRECISION  
          delta_temp_b=-0.05_DEFAULT_PRECISION*delta_temp+0.000001_DEFAULT_PRECISION
        end if
        
        error_temp_a=error_temp(qv,qc,delta_temp_a,theta,exn,p)
        error_temp_b=error_temp(qv,qc,delta_temp_b,theta,exn,p)
       
        if(error_temp_a*error_temp_b>0.0_DEFAULT_PRECISION) then
          ! This should not happen: temperature error has same sign before and after first order adjustment      
          l_fail=.true.
        else
          ! Else: use Ridder's algorithm to find root
          call ridder_temp(qv,qc,delta_temp_a,delta_temp_b,error_temp_a,error_temp_b,theta,exn,p,l_fail,delta_temp)
        end if
        
        ! only keep track of most recent failure within time-step
        if(l_fail) then
          l_any_fail=.true.
          n_fail=n
          qv_fail=qv
          qc_fail=qc
          b_fail=b
          delta_temp=0.0_DEFAULT_PRECISION
          z_fail=z
        end if

        qc=qc+delta_temp*cp_over_rlvap
        qv=qv-delta_temp*cp_over_rlvap
        b=theta+delta_temp/exn-thref0

      end if
      
      ! write back to array          
      current_state%parcels%b(n)=b
      current_state%parcels%qvalues(iqv,n)=qv
      current_state%parcels%qvalues(iqc,n)=qc            
    end do
    
    if(l_any_fail) then
      write(*,*) 'saturation adjustment failed for some parcels' 
      write(*,*) 'for example' 
      write(*,*) 'n_fail=',n_fail
      write(*,*) 'qv_fail=',qv_fail
      write(*,*) 'qc_fail=',qc_fail
      write(*,*) 'b_fail=',b_fail
      write(*,*) 'z_fail=',z_fail
      stop
    end if
    
    call timer_stop(handle)

  end subroutine timestep_callback

end module
