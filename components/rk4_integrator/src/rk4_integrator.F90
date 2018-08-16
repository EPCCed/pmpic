!Runge Kutta 4th order integrator component
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
  real(kind=DEFAULT_PRECISION), parameter :: cfl=0.1
  integer :: ierr

  integer :: rkstep

  integer :: handle, handle1, handle2, handle3, handle4

contains

  type(component_descriptor_type) function rk4_integrator_get_descriptor()
    rk4_integrator_get_descriptor%name="rk4_integrator"
    rk4_integrator_get_descriptor%version=0.1
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
    state%rksteps = 4

    state%parcels%n_rk=12

    !allocate RK variables
    allocate(state%parcels%xo(state%parcels%maxparcels_local))
    allocate(state%parcels%yo(state%parcels%maxparcels_local))
    allocate(state%parcels%zo(state%parcels%maxparcels_local))
    allocate(state%parcels%xf(state%parcels%maxparcels_local))
    allocate(state%parcels%yf(state%parcels%maxparcels_local))
    allocate(state%parcels%zf(state%parcels%maxparcels_local))
    allocate(state%parcels%po(state%parcels%maxparcels_local))
    allocate(state%parcels%qo(state%parcels%maxparcels_local))
    allocate(state%parcels%ro(state%parcels%maxparcels_local))
    allocate(state%parcels%pf(state%parcels%maxparcels_local))
    allocate(state%parcels%qf(state%parcels%maxparcels_local))
    allocate(state%parcels%rf(state%parcels%maxparcels_local))

    rkstep=1

    call register_routine_for_timing("RK4_all_steps",handle,state)
    call register_routine_for_timing("RK4_step_1",handle1,state)
    call register_routine_for_timing("RK4_step_2",handle2,state)
    call register_routine_for_timing("RK4_step_3",handle3,state)
    call register_routine_for_timing("RK4_step_4",handle4,state)


  end subroutine


!solves dY/dt = F(Y(t,Y)) by calculating:
!
! Y(t+h) = Y(t) + 1/6 (k1 + 2k2 + 2k3 + k4)
!
!where:
! k1 = h f(Y(t    , Y))
! k2 = h f(Y(t+h/2, Y+k1/2))
! k3 = h f(Y(t+h/2, Y+k2/2))
! k4 = h f(Y(t+h  , Y+k3))
!
! the xo, yo, zo etc... variables cache the initial values of x, y, z
! the xf, yf, zf etc... variables store the subtotal for the sum of kn
! x, y, z etc contain Y(t,Y)
! dxdt, dydt, dzdt etc are essentially f(t,Y)
  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state

    integer(kind=PARCEL_INTEGER) :: nparcels, n

    real(kind=DEFAULT_PRECISION) :: dt, dt6, dt3, dt2

    ! define some shorthand
    real(kind=DEFAULT_PRECISION),pointer,dimension(:) :: x,y,z,p,q,r
    real(kind=DEFAULT_PRECISION),pointer,dimension(:) :: xo,yo,zo,po,qo,ro
    real(kind=DEFAULT_PRECISION),pointer,dimension(:) :: xf,yf,zf,pf,qf,rf
    real(kind=DEFAULT_PRECISION),pointer,dimension(:) :: dxdt,dydt,dzdt,dpdt,dqdt,drdt
    x=>state%parcels%x
    y=>state%parcels%y
    z=>state%parcels%z
    p=>state%parcels%p
    q=>state%parcels%q
    r=>state%parcels%r
    xo=>state%parcels%xo
    yo=>state%parcels%yo
    zo=>state%parcels%zo
    po=>state%parcels%po
    qo=>state%parcels%qo
    ro=>state%parcels%ro
    xf=>state%parcels%xf
    yf=>state%parcels%yf
    zf=>state%parcels%zf
    pf=>state%parcels%pf
    qf=>state%parcels%qf
    rf=>state%parcels%rf
    dxdt=>state%parcels%dxdt
    dydt=>state%parcels%dydt
    dzdt=>state%parcels%dzdt
    dpdt=>state%parcels%dpdt
    dqdt=>state%parcels%dqdt
    drdt=>state%parcels%drdt

    dt=state%dtm
    dt6=dt/6
    dt2=dt/2
    dt3=dt/3

    nparcels=state%parcels%numparcels_local
    !$OMP PARALLEL
    if (rkstep .eq. 1) then
      !$OMP SINGLE

      call timer_start(handle)
      call timer_start(handle1)

      if (state%parallel%my_rank .eq. 0) write(*,"('RK4 integrator: t= ',f10.2,'s -> ',f10.2,'s')") state%time, state%time+dt
      !$OMP END SINGLE

      !cache initial parcel positions/vorticities
      !$OMP WORKSHARE
      xo(1:nparcels) = x(1:nparcels)
      yo(1:nparcels) = y(1:nparcels)
      zo(1:nparcels) = z(1:nparcels)

      po(1:nparcels) = p(1:nparcels)
      qo(1:nparcels) = q(1:nparcels)
      ro(1:nparcels) = r(1:nparcels)
      !$OMP END WORKSHARE

      !dxdt * dt = k1

      !add 1/6*k1 to xf
      !k1 is just dt * dxdt so this is already calculated by the previous components called
      !$OMP WORKSHARE
      xf(1:nparcels) = xo(1:nparcels)+dt6*dxdt(1:nparcels)
      yf(1:nparcels) = yo(1:nparcels)+dt6*dydt(1:nparcels)
      zf(1:nparcels) = zo(1:nparcels)+dt6*dzdt(1:nparcels)
      pf(1:nparcels) = po(1:nparcels)+dt6*dpdt(1:nparcels)
      qf(1:nparcels) = qo(1:nparcels)+dt6*dqdt(1:nparcels)
      rf(1:nparcels) = ro(1:nparcels)+dt6*drdt(1:nparcels)

      !update Y to position inside the k2 bracket
      x(1:nparcels) = xo(1:nparcels)+dt2*dxdt(1:nparcels)
      y(1:nparcels) = yo(1:nparcels)+dt2*dydt(1:nparcels)
      z(1:nparcels) = zo(1:nparcels)+dt2*dzdt(1:nparcels)
      p(1:nparcels) = po(1:nparcels)+dt2*dpdt(1:nparcels)
      q(1:nparcels) = qo(1:nparcels)+dt2*dqdt(1:nparcels)
      r(1:nparcels) = ro(1:nparcels)+dt2*drdt(1:nparcels)
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
      !dxdt * dt = k2

      !add 1/3*k2 to xf
      !k1 is just dt * dxdt so this is already calculated by the previous components called
      !$OMP WORKSHARE
      xf(1:nparcels) = xf(1:nparcels)+dt3*dxdt(1:nparcels)
      yf(1:nparcels) = yf(1:nparcels)+dt3*dydt(1:nparcels)
      zf(1:nparcels) = zf(1:nparcels)+dt3*dzdt(1:nparcels)
      pf(1:nparcels) = pf(1:nparcels)+dt3*dpdt(1:nparcels)
      qf(1:nparcels) = qf(1:nparcels)+dt3*dqdt(1:nparcels)
      rf(1:nparcels) = rf(1:nparcels)+dt3*drdt(1:nparcels)

      !update Y to position inside the k3 bracket
      x(1:nparcels) = xo(1:nparcels)+dt2*dxdt(1:nparcels)
      y(1:nparcels) = yo(1:nparcels)+dt2*dydt(1:nparcels)
      z(1:nparcels) = zo(1:nparcels)+dt2*dzdt(1:nparcels)
      p(1:nparcels) = po(1:nparcels)+dt2*dpdt(1:nparcels)
      q(1:nparcels) = qo(1:nparcels)+dt2*dqdt(1:nparcels)
      r(1:nparcels) = ro(1:nparcels)+dt2*drdt(1:nparcels)
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

      !dxdt * dt = k3

      !add 1/3*k3 to xf
      !k1 is just dt * dxdt so this is already calculated by the previous components called
      !$OMP WORKSHARE
      xf(1:nparcels) = xf(1:nparcels)+dt3*dxdt(1:nparcels)
      yf(1:nparcels) = yf(1:nparcels)+dt3*dydt(1:nparcels)
      zf(1:nparcels) = zf(1:nparcels)+dt3*dzdt(1:nparcels)
      pf(1:nparcels) = pf(1:nparcels)+dt3*dpdt(1:nparcels)
      qf(1:nparcels) = qf(1:nparcels)+dt3*dqdt(1:nparcels)
      rf(1:nparcels) = rf(1:nparcels)+dt3*drdt(1:nparcels)

      !update Y to position inside the k4 bracket
      x(1:nparcels) = xo(1:nparcels)+dt*dxdt(1:nparcels)
      y(1:nparcels) = yo(1:nparcels)+dt*dydt(1:nparcels)
      z(1:nparcels) = zo(1:nparcels)+dt*dzdt(1:nparcels)
      p(1:nparcels) = po(1:nparcels)+dt*dpdt(1:nparcels)
      q(1:nparcels) = qo(1:nparcels)+dt*dqdt(1:nparcels)
      r(1:nparcels) = ro(1:nparcels)+dt*drdt(1:nparcels)
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

      !dxdt*dt = k4

      !add 1/6*k4 to xf
      !k1 is just dt * dxdt so this is already calculated by the previous components called
      !$OMP WORKSHARE
      xf(1:nparcels) = xf(1:nparcels)+dt6*dxdt(1:nparcels)
      yf(1:nparcels) = yf(1:nparcels)+dt6*dydt(1:nparcels)
      zf(1:nparcels) = zf(1:nparcels)+dt6*dzdt(1:nparcels)
      pf(1:nparcels) = pf(1:nparcels)+dt6*dpdt(1:nparcels)
      qf(1:nparcels) = qf(1:nparcels)+dt6*dqdt(1:nparcels)
      rf(1:nparcels) = rf(1:nparcels)+dt6*drdt(1:nparcels)
      !$OMP END WORKSHARE
      !$OMP BARRIER
      !$OMP SINGLE
      call timer_stop(handle)
      call timer_stop(handle4)
      !$OMP END SINGLE

    else
      print *, "Error: rkstep beyond maximum value (4)"
      print *, "rkstep=", rkstep
      stop
    endif

    !$OMP END PARALLEL

    if (rkstep .lt. 4) then
      rkstep=rkstep+1
    else
      rkstep=1
    endif

    call parcel_haloswap(state)

  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    deallocate(state%parcels%xo)
    deallocate(state%parcels%yo)
    deallocate(state%parcels%zo)
    deallocate(state%parcels%po)
    deallocate(state%parcels%qo)
    deallocate(state%parcels%ro)
    deallocate(state%parcels%xf)
    deallocate(state%parcels%yf)
    deallocate(state%parcels%zf)
    deallocate(state%parcels%pf)
    deallocate(state%parcels%qf)
    deallocate(state%parcels%rf)

  end subroutine



end module
