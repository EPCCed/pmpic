!component that determines the vorticity tendency
! Solves:
! dp/dt = vort . grad(u) + db/dy
! dq=dt = vort . grad(v) - db/dx
! dr/dt = vort . grad(w)
! on a grid (derivatives are calculated spectrally) then interpolates to the parcels
module vorticity_tendency_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use pencil_fft_mod, only : initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  use MPI
  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords, grid2par, par2grid
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use fftops_mod, only: fftops_init, diffx, diffy, diffz, laplinv, spectral_filter
  use prognostics_mod, only: prognostic_field_type
  implicit none

#ifndef TEST_MODE
  private
#endif


  real(kind=DEFAULT_PRECISION) :: PI
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: as, bs !spectral variables
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), ALLOCATABLE :: ap, dudx, dwdx, dvdy, dwdy
  type(prognostic_field_type) :: dp, dq, dr
  integer :: fourier_space_sizes(3)
  integer :: ierr
  integer :: handle
  integer :: nx, ny, nz
  integer :: iteration



  public vorticity_tendency_get_descriptor
contains


  type(component_descriptor_type) function vorticity_tendency_get_descriptor()
    vorticity_tendency_get_descriptor%name="vorticity_tendency"
    vorticity_tendency_get_descriptor%version=0.1
    vorticity_tendency_get_descriptor%initialisation=>initialisation_callback
    vorticity_tendency_get_descriptor%timestep=>timestep_callback
    vorticity_tendency_get_descriptor%finalisation=>finalisation_callback
  end function vorticity_tendency_get_descriptor

  !> This initialisation callback sets up the pencil fft module, allocates data for the fourier space variables
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state


    integer :: my_y_start, my_x_start
    integer :: nxp, nyp, nzp

    iteration=0


    PI=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

    fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)

    !initialise spectral derivatives module
    call fftops_init(current_state,my_x_start,my_y_start,fourier_space_sizes)


    nx=fourier_space_sizes(X_INDEX)
    ny=fourier_space_sizes(Y_INDEX)
    nz=fourier_space_sizes(Z_INDEX)



    allocate(as(nz,ny,nx))
    allocate(bs(nz,ny,nx))


    if (.not. allocated(current_state%u_s%data)) then
      allocate(current_state%u_s%data(nz,ny,nx))
      allocate(current_state%v_s%data(nz,ny,nx))
      allocate(current_state%w_s%data(nz,ny,nx))
    endif

    nzp = current_state%local_grid%size(Z_INDEX)
    nyp = current_state%local_grid%size(Y_INDEX) + 2*current_state%local_grid%halo_size(Y_INDEX)
    nxp = current_state%local_grid%size(X_INDEX) + 2*current_state%local_grid%halo_size(X_INDEX)

    allocate(ap(nzp,nyp,nxp),dudx(nzp,nyp,nxp), dwdx(nzp,nyp,nxp),dvdy(nzp,nyp,nxp),dwdy(nzp,nyp,nxp))
    allocate(dp%data(nzp,nyp,nxp), dq%data(nzp,nyp,nxp),dr%data(nzp,nyp,nxp))





    call register_routine_for_timing("vort_tend",handle,current_state)

  end subroutine initialisation_callback







  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i, xi, xf, zi, zf, yi, yf, j, k
    integer(kind=PARCEL_INTEGER):: n

    real(kind=DEFAULT_PRECISION) :: omax, omaxglobal, dtmax

    !print *, "Entering vorticity_tendency"
    call timer_start(handle)



    !determine total parcel buoyancy
    !$OMP PARALLEL DO
    do n=1,current_state%parcels%numparcels_local
      current_state%parcels%btot(n) = current_state%parcels%b(n) + 12.5 &
      * max(0.0, current_state%parcels%h(n)-exp(-current_state%parcels%z(n)))
    enddo
    !$OMP END PARALLEL DO


    call par2grid(current_state,current_state%parcels%btot,current_state%b)

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP SINGLE
    do i=1,3
      start_loc(i)=current_state%local_grid%local_domain_start_index(i)
      end_loc(i)=current_state%local_grid%local_domain_end_index(i)
    end do
    zi=start_loc(Z_INDEX)
    zf=end_loc(Z_INDEX)
    yi=start_loc(Y_INDEX)
    yf=end_loc(Y_INDEX)
    xi=start_loc(X_INDEX)
    xf=end_loc(X_INDEX)

  !$OMP END SINGLE

    ! get fft of b and put it in bs
    call perform_forward_3dfft(current_state, current_state%b%data(zi:zf, &
         yi:yf, xi:xf), bs)


    ! x component:
    ! dp/dt = vort . grad(u) + db/dy = p*du/dx + q*du/dy + r*du/dz +db/dy
    ! ** but ** we express dw/dz as q + dw/dx (we can do this since q = du/dz - dw/dx)

    call diffy(bs,as) !db/dy (spectral)
    ! take db/dy and convert to positional space, put into the dp array
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dp%data(zi:zf,yi:yf,xi:xf))
  !  !$OMP END SINGLE

    !calculate du/dx (spectral)

    call diffx(current_state%u_s%data,as)
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dudx(zi:zf,yi:yf,xi:xf))
  !  !$OMP END SINGLE
    !add p*du/dx to dp
    !$OMP WORKSHARE
    dp%data(zi:zf,yi:yf,xi:xf) = dp%data(zi:zf,yi:yf,xi:xf) &
     + current_state%p%data(zi:zf,yi:yf,xi:xf)*dudx(zi:zf,yi:yf,xi:xf)
     !$OMP END WORKSHARE

    !calculate du/dy (spectral)
    call diffy(current_state%u_s%data,as)
!    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
  !  !$OMP END SINGLE
    !add q*du/dy to dp
    !$OMP WORKSHARE
    dp%data(zi:zf,yi:yf,xi:xf) = dp%data(zi:zf,yi:yf,xi:xf) &
      + current_state%q%data(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dw/dx (spectral)
    call diffx(current_state%w_s%data,as)
    !add r*(q + dw/dx) to dp
    !!$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dwdx(zi:zf,yi:yf,xi:xf))
    !!$OMP END SINGLE
    !$OMP WORKSHARE
    dp%data(zi:zf,yi:yf,xi:xf) = dp%data(zi:zf,yi:yf,xi:xf) &
      + current_state%r%data(zi:zf,yi:yf,xi:xf) &
      * ( current_state%q%data(zi:zf,yi:yf,xi:xf) + dwdx(zi:zf,yi:yf,xi:xf))
    !$OMP END WORKSHARE

    !y component:
    !dq/dt = vort . grad(v) - db/dx = p*dv/dx + q*dv/dy + v*dv/dz -db/dx
    !! **but** we express dv/dz as dw/dy-p  (since p = dw/dy - dv/dz)

    call diffx(bs,as) !db/dy (spectral)
    ! take db/dx and convert to positional space, put into the dq array
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dq%data(zi:zf,yi:yf,xi:xf))
    !!$OMP END SINGLE

    !calculate dv/dx (spectral)
    call diffx(current_state%v_s%data,as)
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
  !  !$OMP END SINGLE
    !add p*dv/dx to dq
    !$OMP WORKSHARE
    dq%data(zi:zf,yi:yf,xi:xf) = - dq%data(zi:zf,yi:yf,xi:xf) &
     + current_state%p%data(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dv/dy (spectral)
    call diffy(current_state%v_s%data,as)
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dvdy(zi:zf,yi:yf,xi:xf))
  !  !$OMP END SINGLE
    !add p*dv/dy to dq
    !$OMP WORKSHARE
    dq%data(zi:zf,yi:yf,xi:xf) = dq%data(zi:zf,yi:yf,xi:xf) &
      + current_state%q%data(zi:zf,yi:yf,xi:xf)*dvdy(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dw/dy (spectral)
    call diffy(current_state%w_s%data,as)
    !add r*(dw/dy - p) to dq
  !  !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dwdy(zi:zf,yi:yf,xi:xf))
!    !$OMP END SINGLE
    !$OMP WORKSHARE
    dq%data(zi:zf,yi:yf,xi:xf) = dq%data(zi:zf,yi:yf,xi:xf) &
      + current_state%r%data(zi:zf,yi:yf,xi:xf) &
      * ( -current_state%p%data(zi:zf,yi:yf,xi:xf) + dwdy(zi:zf,yi:yf,xi:xf))
    !$OMP END WORKSHARE



    ! dr/dt = vort . grad(w) = p*dw/dx + q*dw/dy + r*dw/dz
    ! **but** using dw/dz= -du/dx - dv/dy (since div(velocity) = 0)

    !since we precomputed all these terms previously, we don't need to carry out
    ! any differentiation or inverse FFT

    !$OMP WORKSHARE
    !add p*dw/dx to dr
    dr%data(zi:zf,yi:yf,xi:xf) = current_state%p%data(zi:zf,yi:yf,xi:xf)*dwdx(zi:zf,yi:yf,xi:xf)

    !add q*dw/dy to dr
    dr%data(zi:zf,yi:yf,xi:xf) = dr%data(zi:zf,yi:yf,xi:xf) &
      + current_state%q%data(zi:zf,yi:yf,xi:xf)*dwdy(zi:zf,yi:yf,xi:xf)

    !add r*(-du/dx) to dr
    dr%data(zi:zf,yi:yf,xi:xf) = dr%data(zi:zf,yi:yf,xi:xf) &
      + current_state%r%data(zi:zf,yi:yf,xi:xf) &
      * (-dudx(zi:zf,yi:yf,xi:xf))

    !add r*(-dv/dy) to dr
    dr%data(zi:zf,yi:yf,xi:xf) = dr%data(zi:zf,yi:yf,xi:xf) &
      + current_state%r%data(zi:zf,yi:yf,xi:xf) &
      * (-dvdy(zi:zf,yi:yf,xi:xf))
    !$OMP END WORKSHARE

    !$OMP END PARALLEL


    !interpolate back to parcels


    call grid2par(current_state,dp,current_state%parcels%dpdt)
    call grid2par(current_state,dq,current_state%parcels%dqdt)
    call grid2par(current_state,dr,current_state%parcels%drdt)

    if ( mod(iteration,current_state%rksteps) == 0) then
      !We now want to determine the maximum vorticity
      omax = maxval(current_state%p%data(zi:zf,yi:yf,xi:xf)**2 &
                  + current_state%q%data(zi:zf,yi:yf,xi:xf)**2 &
                  + current_state%r%data(zi:zf,yi:yf,xi:xf)**2)
      omax=sqrt(omax)

      !This is the local maximum. We want the global maximum so we do a MPI reduction operation
      call MPI_Allreduce(omax,&
                         omaxglobal,&
                         1,&
                         PRECISION_TYPE,&
                         MPI_MAX,&
                         current_state%parallel%monc_communicator,&
                         ierr)

      !this is the maximum timestep in vorticity
      if (omaxglobal .gt. 0.) then
        dtmax = 0.5/omaxglobal
      else
        dtmax=current_state%dtmax
      endif

      if (current_state%dtm .gt. dtmax) then
        current_state%dtm = dtmax
        if (current_state%parallel%my_rank .eq. 0) then
           print*, "Timestep limited by vorticity", current_state%dtm
        endif
      else if (current_state%dtm .lt. current_state%dtmax) then
        if (current_state%parallel%my_rank .eq. 0) then
           print*, "Timestep limited by velocity", current_state%dtm
        endif
      else
        if (current_state%parallel%my_rank .eq. 0) then
           print*, "Timestep limited by dtmax", current_state%dtm
        endif
      endif



      ! print *, "vorticity tendency"
      ! print *, "dtmax=",dtmax
      !print *, maxval(dp%data), maxval(dq%data), maxval(dr%data)
      !print *, minval(dp%data), minval(dq%data), minval(dr%data)
    endif

    iteration=iteration+1

    call timer_stop(handle)

    !print *, "exiting vorticity_tendency"

    !call MPI_Finalize(ierr)
    !stop



  end subroutine timestep_callback






  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    !call finalise_pencil_fft(current_state%parallel%monc_communicator)
    !deallocate(a, b, c, d, e, f)
    !deallocate(current_state%u_s%data, current_state%v_s%data, current_state%w_s%data)

  end subroutine finalisation_callback



end module
