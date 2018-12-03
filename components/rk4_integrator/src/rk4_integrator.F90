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
  integer :: ierr

  integer :: rkstep

  integer :: handle, handle1, handle2, handle3, handle4
  integer,parameter :: chunksize=1024

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

    dt=state%dtm
    dt6=dt/6
    dt2=dt/2
    dt3=dt/3

    nparcels=state%parcels%numparcels_local

    if (rkstep .eq. 1) then
      !$OMP SINGLE

      call timer_start(handle)
      call timer_start(handle1)

      if (state%parallel%my_rank .eq. 0) write(*,"('RK4 integrator: t= ',f10.2,'s -> ',f10.2,'s')") state%time, state%time+dt
      !$OMP END SINGLE
      call set_equal_to(state%parcels%xo,state%parcels%x,nparcels)
      call set_equal_to(state%parcels%yo,state%parcels%y,nparcels)
      call set_equal_to(state%parcels%zo,state%parcels%z,nparcels)
      call set_equal_to(state%parcels%po,state%parcels%p,nparcels)
      call set_equal_to(state%parcels%qo,state%parcels%q,nparcels)
      call set_equal_to(state%parcels%ro,state%parcels%r,nparcels)

      !dxdt * dt = k1

      !add 1/6*k1 to xf
      !k1 is just dt * dxdt so this is already calculated by the previous components called
      call add_to_other(state%parcels%xf,state%parcels%xo,state%parcels%dxdt,dt6,nparcels)
      call add_to_other(state%parcels%yf,state%parcels%yo,state%parcels%dydt,dt6,nparcels)
      call add_to_other(state%parcels%zf,state%parcels%zo,state%parcels%dzdt,dt6,nparcels)
      call add_to_other(state%parcels%pf,state%parcels%po,state%parcels%dpdt,dt6,nparcels)
      call add_to_other(state%parcels%qf,state%parcels%qo,state%parcels%dqdt,dt6,nparcels)
      call add_to_other(state%parcels%rf,state%parcels%ro,state%parcels%drdt,dt6,nparcels)

      !update Y to position inside the k2 bracket
      call add_to_other(state%parcels%x,state%parcels%xo,state%parcels%dxdt,dt2,nparcels)
      call add_to_other(state%parcels%y,state%parcels%yo,state%parcels%dydt,dt2,nparcels)
      call add_to_other(state%parcels%z,state%parcels%zo,state%parcels%dzdt,dt2,nparcels)
      call add_to_other(state%parcels%p,state%parcels%po,state%parcels%dpdt,dt2,nparcels)
      call add_to_other(state%parcels%q,state%parcels%qo,state%parcels%dqdt,dt2,nparcels)
      call add_to_other(state%parcels%r,state%parcels%ro,state%parcels%drdt,dt2,nparcels)

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

      call add_to_self(state%parcels%xf,state%parcels%dxdt,dt3,nparcels)
      call add_to_self(state%parcels%yf,state%parcels%dydt,dt3,nparcels)
      call add_to_self(state%parcels%zf,state%parcels%dzdt,dt3,nparcels)
      call add_to_self(state%parcels%pf,state%parcels%dpdt,dt3,nparcels)
      call add_to_self(state%parcels%qf,state%parcels%dqdt,dt3,nparcels)
      call add_to_self(state%parcels%rf,state%parcels%drdt,dt3,nparcels)

      !update Y to position inside the k3 bracket
      call add_to_other(state%parcels%x,state%parcels%xo,state%parcels%dxdt,dt2,nparcels)
      call add_to_other(state%parcels%y,state%parcels%yo,state%parcels%dydt,dt2,nparcels)
      call add_to_other(state%parcels%z,state%parcels%zo,state%parcels%dzdt,dt2,nparcels)
      call add_to_other(state%parcels%p,state%parcels%po,state%parcels%dpdt,dt2,nparcels)
      call add_to_other(state%parcels%q,state%parcels%qo,state%parcels%dqdt,dt2,nparcels)
      call add_to_other(state%parcels%r,state%parcels%ro,state%parcels%drdt,dt2,nparcels)

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

      call add_to_self(state%parcels%xf,state%parcels%dxdt,dt3,nparcels)
      call add_to_self(state%parcels%yf,state%parcels%dydt,dt3,nparcels)
      call add_to_self(state%parcels%zf,state%parcels%dzdt,dt3,nparcels)
      call add_to_self(state%parcels%pf,state%parcels%dpdt,dt3,nparcels)
      call add_to_self(state%parcels%qf,state%parcels%dqdt,dt3,nparcels)
      call add_to_self(state%parcels%rf,state%parcels%drdt,dt3,nparcels)

      !update Y to position inside the k4 bracket
      call add_to_other(state%parcels%x,state%parcels%xo,state%parcels%dxdt,dt,nparcels)
      call add_to_other(state%parcels%y,state%parcels%yo,state%parcels%dydt,dt,nparcels)
      call add_to_other(state%parcels%z,state%parcels%zo,state%parcels%dzdt,dt,nparcels)
      call add_to_other(state%parcels%p,state%parcels%po,state%parcels%dpdt,dt,nparcels)
      call add_to_other(state%parcels%q,state%parcels%qo,state%parcels%dqdt,dt,nparcels)
      call add_to_other(state%parcels%r,state%parcels%ro,state%parcels%drdt,dt,nparcels)
      
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
      call add_to_self(state%parcels%xf,state%parcels%dxdt,dt6,nparcels)
      call add_to_self(state%parcels%yf,state%parcels%dydt,dt6,nparcels)
      call add_to_self(state%parcels%zf,state%parcels%dzdt,dt6,nparcels)
      call add_to_self(state%parcels%pf,state%parcels%dpdt,dt6,nparcels)
      call add_to_self(state%parcels%qf,state%parcels%dqdt,dt6,nparcels)
      call add_to_self(state%parcels%rf,state%parcels%drdt,dt6,nparcels)

      ! finally update xf,yf,zf,pf,qf,rf to x,y,z,p,q,r
      call set_equal_to(state%parcels%x,state%parcels%xf,nparcels)
      call set_equal_to(state%parcels%y,state%parcels%yf,nparcels)
      call set_equal_to(state%parcels%z,state%parcels%zf,nparcels)
      call set_equal_to(state%parcels%p,state%parcels%pf,nparcels)
      call set_equal_to(state%parcels%q,state%parcels%qf,nparcels)
      call set_equal_to(state%parcels%r,state%parcels%rf,nparcels)
      
      !$OMP SINGLE
      call timer_stop(handle)
      call timer_stop(handle4)
      !$OMP END SINGLE

    else
      print *, "Error: rkstep beyond maximum value (4)"
      print *, "rkstep=", rkstep
      stop
    endif


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

  subroutine set_equal_to(varout,varin,nparcels)
      real(kind=DEFAULT_PRECISION), dimension(:), intent(out) :: varout
      real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: varin
      integer(kind=PARCEL_INTEGER), intent(in) :: nparcels
      integer(kind=PARCEL_INTEGER) :: n
  
      !$OMP PARALLEL SHARED(varout,varin) & 
      !$OMP& PRIVATE(n) FIRSTPRIVATE(nparcels) DEFAULT(NONE) 
      !$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
      do n=1,nparcels
        varout(n)=varin(n)
      end do
      !$OMP END DO
      !$OMP END PARALLEL     
  end subroutine

  subroutine add_to_self(varinout,varin2,frac,nparcels)
      real(kind=DEFAULT_PRECISION), dimension(:), intent(inout) :: varinout
      real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: varin2
      real(kind=DEFAULT_PRECISION), intent(in) :: frac
      integer(kind=PARCEL_INTEGER), intent(in) :: nparcels
      integer(kind=PARCEL_INTEGER) :: n
  
      !$OMP PARALLEL SHARED(varinout,varin2) & 
      !$OMP& PRIVATE(n) FIRSTPRIVATE(nparcels,frac) DEFAULT(NONE) 
      !$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
      do n=1,nparcels
        varinout(n)=varinout(n)+frac*varin2(n)
      end do
      !$OMP END DO
      !$OMP END PARALLEL 
  end subroutine
  
  subroutine add_to_other(varout,varin1,varin2,frac,nparcels)
      real(kind=DEFAULT_PRECISION), dimension(:), intent(out) :: varout
      real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: varin1
      real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: varin2
      real(kind=DEFAULT_PRECISION), intent(in) :: frac
      integer(kind=PARCEL_INTEGER), intent(in) :: nparcels
      integer(kind=PARCEL_INTEGER) :: n
  
      !$OMP PARALLEL SHARED(varout,varin1,varin2) & 
      !$OMP& PRIVATE(n) FIRSTPRIVATE(nparcels,frac) DEFAULT(NONE) 
      !$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
      do n=1,nparcels
        varout(n)=varin1(n)+frac*varin2(n)
      end do
      !$OMP END DO
      !$OMP END PARALLEL  
  end subroutine
      
end module
