!Runge Kutta 4th order integrator component, using low-storage RK method
module rk4_integrator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, dz
  use parcel_haloswap_mod, only: parcel_haloswap
  use MPI
  use timer_mod

  implicit none

  real(kind=DEFAULT_PRECISION) :: originaldt
  integer :: ierr

  integer :: rkstep

  integer :: handle, handle1, handle2, handle3, handle4, handle5
  
  real(kind=DEFAULT_PRECISION), parameter :: cA_1=0.
  real(kind=DEFAULT_PRECISION), parameter :: cA_2=- 567301805773./1357537059087.
  real(kind=DEFAULT_PRECISION), parameter :: cA_3=-2404267990393./2016746695238.
  real(kind=DEFAULT_PRECISION), parameter :: cA_4=-3550918686646./2091501179385.
  real(kind=DEFAULT_PRECISION), parameter :: cA_5=-1275806237668./ 842570457699.

  real(kind=DEFAULT_PRECISION), parameter :: cB_1=1432997174477./ 9575080441755.
  real(kind=DEFAULT_PRECISION), parameter :: cB_2=5161836677717./13612068292357.
  real(kind=DEFAULT_PRECISION), parameter :: cB_3=1720146321549./ 2090206949498.
  real(kind=DEFAULT_PRECISION), parameter :: cB_4=3134564353537./ 4481467310338.
  real(kind=DEFAULT_PRECISION), parameter :: cB_5=2277821191437./14882151754819.

contains

  type(component_descriptor_type) function rk4_integrator_get_descriptor()
    rk4_integrator_get_descriptor%name="rk4_integrator"
    rk4_integrator_get_descriptor%version=0.2
    rk4_integrator_get_descriptor%initialisation=>initialisation_callback
    rk4_integrator_get_descriptor%timestep=>timestep_callback
    rk4_integrator_get_descriptor%finalisation=>finalisation_callback
  end function rk4_integrator_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    if (state%parallel%my_rank .eq. 0) print *, "In RK4 Integrator Initialisation"

    originaldt= options_get_real(state%options_database,"dtm")
    state%dtmax = options_get_real(state%options_database,"dtmax")

    state%dtm = originaldt
    !at present we don't want to update the timestep
    state%update_dtm = .false.

    if (state%parallel%my_rank .eq. 0) print *, "Starting dt=",originaldt

    if (state%rksteps .ne. 0) then
      print *, "Error - another integrator component is present. ABORTING"
      stop
    endif
    state%rksteps = 5

    state%parcels%n_rk=0

    rkstep=1

    call register_routine_for_timing("RK4_all_steps",handle,state)
    call register_routine_for_timing("RK4_step_1",handle1,state)
    call register_routine_for_timing("RK4_step_2",handle2,state)
    call register_routine_for_timing("RK4_step_3",handle3,state)
    call register_routine_for_timing("RK4_step_4",handle4,state)
    call register_routine_for_timing("RK4_step_5",handle5,state)


  end subroutine


! RK4 attempt using low storage (Carpenter and Kennedy) approach

  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    integer(kind=PARCEL_INTEGER) :: nparcels, n
    real(kind=DEFAULT_PRECISION) :: dt
        
    dt=state%dtm

    nparcels=state%parcels%numparcels_local
    !$OMP PARALLEL
    if (rkstep .eq. 1) then
      !$OMP SINGLE

      call timer_start(handle)
      call timer_start(handle1)

      if (state%parallel%my_rank .eq. 0) write(*,"('RK4 integrator: t= ',f10.2,'s -> ',f10.2,'s')") state%time, state%time+dt
      !$OMP END SINGLE

      !$OMP WORKSHARE
      state%parcels%x(1:nparcels) = state%parcels%x(1:nparcels) &
                                  + cB_1*dt*state%parcels%dxdt(1:nparcels)
      state%parcels%y(1:nparcels) = state%parcels%y(1:nparcels) &
                                  + cB_1*dt*state%parcels%dydt(1:nparcels)
      state%parcels%z(1:nparcels) = state%parcels%z(1:nparcels) &
                                  + cB_1*dt*state%parcels%dzdt(1:nparcels)
      state%parcels%p(1:nparcels) = state%parcels%p(1:nparcels) &
                                  + cB_1*dt*state%parcels%dpdt(1:nparcels)
      state%parcels%q(1:nparcels) = state%parcels%q(1:nparcels) &
                                  + cB_1*dt*state%parcels%dqdt(1:nparcels)
      state%parcels%r(1:nparcels) = state%parcels%r(1:nparcels) &
                                  + cB_1*dt*state%parcels%drdt(1:nparcels)
                                  
      state%parcels%dxdt(1:nparcels) = cA_1*state%parcels%dxdt(1:nparcels)
      state%parcels%dydt(1:nparcels) = cA_1*state%parcels%dydt(1:nparcels)
      state%parcels%dzdt(1:nparcels) = cA_1*state%parcels%dzdt(1:nparcels)
      state%parcels%dpdt(1:nparcels) = cA_1*state%parcels%dpdt(1:nparcels)
      state%parcels%dqdt(1:nparcels) = cA_1*state%parcels%dqdt(1:nparcels)
      state%parcels%drdt(1:nparcels) = cA_1*state%parcels%drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      
      !$OMP SINGLE
      call timer_pause(handle)
      call timer_stop(handle1)
      !$OMP END SINGLE


    else if (rkstep .eq. 2) then
      !$OMP SINGLE

      call timer_resume(handle)
      call timer_start(handle2)
      !if (state%parallel%my_rank .eq. 0) print *, "rkstep 2: t=", state%time
      !$OMP END SINGLE
      
      !$OMP WORKSHARE
      state%parcels%x(1:nparcels) = state%parcels%x(1:nparcels) &
                                  + cB_2*dt*state%parcels%dxdt(1:nparcels)
      state%parcels%y(1:nparcels) = state%parcels%y(1:nparcels) &
                                  + cB_2*dt*state%parcels%dydt(1:nparcels)
      state%parcels%z(1:nparcels) = state%parcels%z(1:nparcels) &
                                  + cB_2*dt*state%parcels%dzdt(1:nparcels)
      state%parcels%p(1:nparcels) = state%parcels%p(1:nparcels) &
                                  + cB_2*dt*state%parcels%dpdt(1:nparcels)
      state%parcels%q(1:nparcels) = state%parcels%q(1:nparcels) &
                                  + cB_2*dt*state%parcels%dqdt(1:nparcels)
      state%parcels%r(1:nparcels) = state%parcels%r(1:nparcels) &
                                  + cB_2*dt*state%parcels%drdt(1:nparcels)
                                  
      state%parcels%dxdt(1:nparcels) = cA_2*state%parcels%dxdt(1:nparcels)
      state%parcels%dydt(1:nparcels) = cA_2*state%parcels%dydt(1:nparcels)
      state%parcels%dzdt(1:nparcels) = cA_2*state%parcels%dzdt(1:nparcels)
      state%parcels%dpdt(1:nparcels) = cA_2*state%parcels%dpdt(1:nparcels)
      state%parcels%dqdt(1:nparcels) = cA_2*state%parcels%dqdt(1:nparcels)
      state%parcels%drdt(1:nparcels) = cA_2*state%parcels%drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      
      !$OMP SINGLE
      call timer_pause(handle)
      call timer_stop(handle2)
      !$OMP END SINGLE

    else if (rkstep .eq. 3) then
      !$OMP SINGLE
      call timer_resume(handle)
      call timer_start(handle3)
      !if (state%parallel%my_rank .eq. 0) print *, "rkstep 3: t=", state%time
      !$OMP END SINGLE
      
      !$OMP WORKSHARE
      state%parcels%x(1:nparcels) = state%parcels%x(1:nparcels) &
                                  + cB_3*dt*state%parcels%dxdt(1:nparcels)
      state%parcels%y(1:nparcels) = state%parcels%y(1:nparcels) &
                                  + cB_3*dt*state%parcels%dydt(1:nparcels)
      state%parcels%z(1:nparcels) = state%parcels%z(1:nparcels) &
                                  + cB_3*dt*state%parcels%dzdt(1:nparcels)
      state%parcels%p(1:nparcels) = state%parcels%p(1:nparcels) &
                                  + cB_3*dt*state%parcels%dpdt(1:nparcels)
      state%parcels%q(1:nparcels) = state%parcels%q(1:nparcels) &
                                  + cB_3*dt*state%parcels%dqdt(1:nparcels)
      state%parcels%r(1:nparcels) = state%parcels%r(1:nparcels) &
                                  + cB_3*dt*state%parcels%drdt(1:nparcels)
                                  
      state%parcels%dxdt(1:nparcels) = cA_3*state%parcels%dxdt(1:nparcels)
      state%parcels%dydt(1:nparcels) = cA_3*state%parcels%dydt(1:nparcels)
      state%parcels%dzdt(1:nparcels) = cA_3*state%parcels%dzdt(1:nparcels)
      state%parcels%dpdt(1:nparcels) = cA_3*state%parcels%dpdt(1:nparcels)
      state%parcels%dqdt(1:nparcels) = cA_3*state%parcels%dqdt(1:nparcels)
      state%parcels%drdt(1:nparcels) = cA_3*state%parcels%drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      
      !$OMP SINGLE
      call timer_pause(handle)
      call timer_stop(handle3)
      !$OMP END SINGLE


    else if (rkstep .eq. 4) then
      !$OMP SINGLE
      call timer_resume(handle)
      call timer_start(handle4)
      !if (state%parallel%my_rank .eq. 0) print *, "rkstep 4: t=", state%time
      !$OMP END SINGLE
      
      !$OMP WORKSHARE
      state%parcels%x(1:nparcels) = state%parcels%x(1:nparcels) &
                                  + cB_4*dt*state%parcels%dxdt(1:nparcels)
      state%parcels%y(1:nparcels) = state%parcels%y(1:nparcels) &
                                  + cB_4*dt*state%parcels%dydt(1:nparcels)
      state%parcels%z(1:nparcels) = state%parcels%z(1:nparcels) &
                                  + cB_4*dt*state%parcels%dzdt(1:nparcels)
      state%parcels%p(1:nparcels) = state%parcels%p(1:nparcels) &
                                  + cB_4*dt*state%parcels%dpdt(1:nparcels)
      state%parcels%q(1:nparcels) = state%parcels%q(1:nparcels) &
                                  + cB_4*dt*state%parcels%dqdt(1:nparcels)
      state%parcels%r(1:nparcels) = state%parcels%r(1:nparcels) &
                                  + cB_4*dt*state%parcels%drdt(1:nparcels)
                                  
      state%parcels%dxdt(1:nparcels) = cA_4*state%parcels%dxdt(1:nparcels)
      state%parcels%dydt(1:nparcels) = cA_4*state%parcels%dydt(1:nparcels)
      state%parcels%dzdt(1:nparcels) = cA_4*state%parcels%dzdt(1:nparcels)
      state%parcels%dpdt(1:nparcels) = cA_4*state%parcels%dpdt(1:nparcels)
      state%parcels%dqdt(1:nparcels) = cA_4*state%parcels%dqdt(1:nparcels)
      state%parcels%drdt(1:nparcels) = cA_4*state%parcels%drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      
      !$OMP SINGLE
      call timer_pause(handle)
      call timer_stop(handle4)
      !$OMP END SINGLE
    else if (rkstep .eq. 5) then
      !$OMP SINGLE
      call timer_resume(handle)
      call timer_start(handle5)
      !if (state%parallel%my_rank .eq. 0) print *, "rkstep 4: t=", state%time
      !$OMP END SINGLE
      
      !$OMP WORKSHARE
      state%parcels%x(1:nparcels) = state%parcels%x(1:nparcels) &
                                  + cB_5*dt*state%parcels%dxdt(1:nparcels)
      state%parcels%y(1:nparcels) = state%parcels%y(1:nparcels) &
                                  + cB_5*dt*state%parcels%dydt(1:nparcels)
      state%parcels%z(1:nparcels) = state%parcels%z(1:nparcels) &
                                  + cB_5*dt*state%parcels%dzdt(1:nparcels)
      state%parcels%p(1:nparcels) = state%parcels%p(1:nparcels) &
                                  + cB_5*dt*state%parcels%dpdt(1:nparcels)
      state%parcels%q(1:nparcels) = state%parcels%q(1:nparcels) &
                                  + cB_5*dt*state%parcels%dqdt(1:nparcels)
      state%parcels%r(1:nparcels) = state%parcels%r(1:nparcels) &
                                  + cB_5*dt*state%parcels%drdt(1:nparcels)
                                  
      state%parcels%dxdt(1:nparcels) = cA_5*state%parcels%dxdt(1:nparcels)
      state%parcels%dydt(1:nparcels) = cA_5*state%parcels%dydt(1:nparcels)
      state%parcels%dzdt(1:nparcels) = cA_5*state%parcels%dzdt(1:nparcels)
      state%parcels%dpdt(1:nparcels) = cA_5*state%parcels%dpdt(1:nparcels)
      state%parcels%dqdt(1:nparcels) = cA_5*state%parcels%dqdt(1:nparcels)
      state%parcels%drdt(1:nparcels) = cA_5*state%parcels%drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      
      !$OMP SINGLE
      call timer_stop(handle)
      call timer_stop(handle5)
      !$OMP END SINGLE
    else
      print *, "Error: rkstep beyond maximum value (5)"
      print *, "rkstep=", rkstep
      stop
    endif

    !$OMP END PARALLEL

    if (rkstep .lt. 5) then
      rkstep=rkstep+1
    else
      rkstep=1
    endif

    call parcel_haloswap(state)

  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

  end subroutine



end module
