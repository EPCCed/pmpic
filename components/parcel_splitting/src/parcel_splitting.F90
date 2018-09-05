!Runge Kutta 4th order integrator component
module parcel_splitting_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER, MPI_PARCEL_INT
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

  integer :: iteration=0
  real(kind=DEFAULT_PRECISION), parameter :: stretchmax = 4.
  real(kind=DEFAULT_PRECISION), parameter :: pi=4.*atan(1.d0)

  integer :: handle



contains

  type(component_descriptor_type) function parcel_splitting_get_descriptor()
    parcel_splitting_get_descriptor%name="parcel_splitting"
    parcel_splitting_get_descriptor%version=0.1
    parcel_splitting_get_descriptor%initialisation=>initialisation_callback
    parcel_splitting_get_descriptor%timestep=>timestep_callback
    parcel_splitting_get_descriptor%finalisation=>finalisation_callback
  end function parcel_splitting_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    call register_routine_for_timing("splitting", handle, state)


  end subroutine



  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer(kind=PARCEL_INTEGER) :: i, n, oldtotal
    real(kind=DEFAULT_PRECISION) :: dt, dstretch, volg, r, absw

    volg = state%global_grid%resolution(1)*state%global_grid%resolution(2) &
            *state%global_grid%resolution(3)

    dt=state%dtm

    if (mod(iteration,state%rksteps) == state%rksteps-1) then
      call timer_start(handle)

      n=state%parcels%numparcels_local

      do i=1,state%parcels%numparcels_local

        !update stretch

        dstretch =  state%parcels%p(i)*state%parcels%dpdt(i) &
                  + state%parcels%q(i)*state%parcels%dqdt(i) &
                  + state%parcels%r(i)*state%parcels%drdt(i)
        dstretch = abs(dstretch) ** (1./3.)

        state%parcels%stretch(i) = state%parcels%stretch(i) + dstretch*dt

        !see if the stretch is above the limit. If so, split the parcel
        if ((state%parcels%stretch(i) .gt. stretchmax) .and. (state%parcels%vol(i) .gt. 1./6.**3)) then
          !determine the half separation of the two new parcels
          r = state%parcels%vol(i)*volg/4./pi
          r=r**(1./3.)


          !determine the absolute value of the vorticity
          absw = state%parcels%p(i)**2 + state%parcels%q(i)**2 + state%parcels%r(i)**2
          absw=sqrt(absw)

          !print*, "r=",r," dx,dy,dz=", r*state%parcels%p(i)/absw,r*state%parcels%q(i)/absw,r*state%parcels%r(i)/absw


          !create new parcel at (x,y,z) + r*(p,q,r)/absw
          n = n+1
          state%parcels%x(n) = state%parcels%x(i) + r*state%parcels%p(i)/absw
          state%parcels%y(n) = state%parcels%y(i) + r*state%parcels%q(i)/absw
          state%parcels%z(n) = state%parcels%z(i) + r*state%parcels%r(i)/absw

          state%parcels%p(n) = state%parcels%p(i)
          state%parcels%q(n) = state%parcels%q(i)
          state%parcels%r(n) = state%parcels%r(i)

          state%parcels%dxdt(n) = state%parcels%dxdt(i)
          state%parcels%dydt(n) = state%parcels%dydt(i)
          state%parcels%dzdt(n) = state%parcels%dzdt(i)

          state%parcels%dpdt(n) = state%parcels%dpdt(i)
          state%parcels%dqdt(n) = state%parcels%dqdt(i)
          state%parcels%drdt(n) = state%parcels%drdt(i)

          state%parcels%h(n) = state%parcels%h(i)
          state%parcels%b(n) = state%parcels%b(i)
          state%parcels%vol(n) = state%parcels%vol(i)/2

          state%parcels%stretch(n) = 0.

          state%parcels%tag(n) = state%parcels%tag(i)
          state%parcels%qvalues(:,n) = state%parcels%qvalues(:,i)

          !shift the original parcel by -r and half its volume
          state%parcels%x(i) = state%parcels%x(i) - r*state%parcels%p(i)/absw
          state%parcels%y(i) = state%parcels%y(i) - r*state%parcels%q(i)/absw
          state%parcels%z(i) = state%parcels%z(i) - r*state%parcels%r(i)/absw

          state%parcels%vol(i) = state%parcels%vol(i)/2

          state%parcels%stretch(i) = 0.

          !  print *, i,state%parcels%x(i), state%parcels%y(i), state%parcels%z(i) &
        !           ,  n,state%parcels%x(n), state%parcels%y(n), state%parcels%z(n)

        endif

      enddo


      !update number of locally held parcels to include newly split parcels
      state%parcels%numparcels_local = n

      oldtotal = state%parcels%numparcels_global
      !update total number of parcels
      call MPI_Allreduce(state%parcels%numparcels_local,&
                         state%parcels%numparcels_global,&
                         1,&
                         MPI_PARCEL_INT,&
                         MPI_SUM,&
                         state%parallel%monc_communicator,&
                         ierr)

      if (state%parallel%my_rank .eq. 0) then
        if (state%parcels%numparcels_global .gt. oldtotal) then
          write(*,"('split ', i6, ' parcels')") state%parcels%numparcels_global - oldtotal
        else
          write(*,*) "No split parcels"
        endif
      endif

    call timer_stop(handle)

    endif


    iteration = iteration + 1


  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state


  end subroutine



end module
