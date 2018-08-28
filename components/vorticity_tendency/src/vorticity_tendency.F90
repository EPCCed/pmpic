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
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), ALLOCATABLE :: ap
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

    allocate(ap(nzp,nyp,nxp))
    allocate(current_state%dp%data(nzp,nyp,nxp), current_state%dq%data(nzp,nyp,nxp),current_state%dr%data(nzp,nyp,nxp))





    call register_routine_for_timing("vort_tend",handle,current_state)

  end subroutine initialisation_callback







  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i, xi, xf, zi, zf, yi, yf, j, k
    integer(kind=PARCEL_INTEGER):: n

    real(kind=DEFAULT_PRECISION) :: omax, omaxglobal, dtmax

    ! define some shorthand
    real(kind=DEFAULT_PRECISION),pointer,dimension(:,:,:) :: p,q,r,dp,dq,dr,hgliq,hg,u_s,v_s,w_s
    real(kind=DEFAULT_PRECISION),pointer,dimension(:) :: btot,b,h,z
    btot=>current_state%parcels%btot
    b=>current_state%parcels%b
    h=>current_state%parcels%h
    z=>current_state%parcels%z
    
    p=>current_state%p%data
    q=>current_state%q%data
    r=>current_state%r%data
    hg=>current_state%hg%data
    hgliq=>current_state%hgliq%data
    u_s=>current_state%u_s%data
    v_s=>current_state%v_s%data
    w_s=>current_state%w_s%data
    dp=>current_state%dp%data
    dq=>current_state%dq%data
    dr=>current_state%dr%data
    
    !print *, "Entering vorticity_tendency"
    call timer_start(handle)

    !$OMP PARALLEL DEFAULT(SHARED)

    !determine total parcel buoyancy
    !$OMP DO
    do n=1,current_state%parcels%numparcels_local
      btot(n) = b(n) + 12.5*max(0.0,h(n)-exp(-z(n)))
    enddo
    !$OMP END DO

    !$OMP SINGLE
    call par2grid(current_state,btot,current_state%b)
    call par2grid(current_state,h,current_state%hg)
    !$OMP END SINGLE

    !$OMP DO
    call par2grid(current_state,btot,current_state%b)
    call par2grid(current_state,h,current_state%hg)
    do i=1,size(hg,3)
      do j=1,size(hg,2)
        do k=1,size(hg,1)
          hgliq(k,j,i) = max(0.,hg(k,j,i) -exp(-z_coords(k)))
        enddo
      enddo
    enddo
    !$OMP END DO
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



    ! get fft of b and put it in bs
    call perform_forward_3dfft(current_state, current_state%b%data(zi:zf, &
         yi:yf, xi:xf), bs)
    !$OMP END SINGLE

    ! x component:
    ! dp/dt = vort . grad(u) + db/dy = p*du/dx + q*du/dy + r*du/dz +db/dy
    ! ** but ** we express dw/dz as q + dw/dz (we can do this since q = du/dz - dw/dx)

    call diffy(bs,as) !db/dy (spectral)
    ! take db/dy and convert to positional space, put into the dp array

    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dp(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE

    !calculate du/dx (spectral)

    call diffx(u_s,as)
    !$OMP SINGLE

    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*du/dx to dp

    !$OMP WORKSHARE
    dp(zi:zf,yi:yf,xi:xf) = dp(zi:zf,yi:yf,xi:xf) &
     + p(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
     !$OMP END WORKSHARE
     
    !calculate du/dy (spectral)
    call diffy(u_s,as)
    !$OMP SINGLE

    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*du/dx to dp

    !$OMP WORKSHARE
    dp(zi:zf,yi:yf,xi:xf) = dp(zi:zf,yi:yf,xi:xf) &
      + q(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dw/dx (spectral)
    call diffx(w_s,as)
    !add r*(q + dw/dx) to dp
    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))

    !$OMP END SINGLE
    !$OMP WORKSHARE
    dp(zi:zf,yi:yf,xi:xf) = dp(zi:zf,yi:yf,xi:xf) &
      + r(zi:zf,yi:yf,xi:xf)*( q(zi:zf,yi:yf,xi:xf) + ap(zi:zf,yi:yf,xi:xf))
    !$OMP END WORKSHARE

    !y component:
    !dq/dt = vort . grad(v) - db/dx = p*dv/dx + q*dv/dy + v*dv/dz -db/dx
    !! **but** we express dv/dz as dw/dy-p  (since p = dw/dy - dv/dz)

    call diffx(bs,as) !db/dy (spectral)
    ! take db/dx and convert to positional space, put into the dq array

    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, dq(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE

    !calculate dv/dx (spectral)
    call diffx(v_s,as)
    !$OMP SINGLE

    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*dv/dx to dq

    !$OMP WORKSHARE
    dq(zi:zf,yi:yf,xi:xf) = - dq(zi:zf,yi:yf,xi:xf) &
     + p(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dv/dy (spectral)
    call diffy(v_s,as)
    !$OMP SINGLE

    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*dv/dx to dq

    !$OMP WORKSHARE
    dq(zi:zf,yi:yf,xi:xf) = dq(zi:zf,yi:yf,xi:xf) &
      + q(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dw/dy (spectral)
    call diffy(w_s,as)
    !add r*(dw/dy - p) to dq
    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))

    dq(zi:zf,yi:yf,xi:xf) = dq(zi:zf,yi:yf,xi:xf) &
      + r(zi:zf,yi:yf,xi:xf)*( -p(zi:zf,yi:yf,xi:xf) + ap(zi:zf,yi:yf,xi:xf))

    !$OMP END SINGLE
    !$OMP WORKSHARE
    dq(zi:zf,yi:yf,xi:xf) = dq(zi:zf,yi:yf,xi:xf) &
      + r(zi:zf,yi:yf,xi:xf)*( -p(zi:zf,yi:yf,xi:xf) + ap(zi:zf,yi:yf,xi:xf))
    !$OMP END WORKSHARE


    ! dr/dt = vort . grad(w) = p*dw/dx + q*dw/dy + r*dw/dz
    ! **but** using dw/dz= -du/dx - dv/dy (since div(velocity) = 0)

    !calculate dw/dx (spectral)

    call diffx(w_s,as)
    !$OMP SINGLE

    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*dw/dx to dr

    !$OMP WORKSHARE
    dr(zi:zf,yi:yf,xi:xf) = p(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dw/dy (spectral)
    call diffy(w_s,as)
    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))
    !$OMP END SINGLE
    !add p*dw/dx to dr

    !$OMP WORKSHARE
    dr(zi:zf,yi:yf,xi:xf) = dr(zi:zf,yi:yf,xi:xf) &
      + q(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate du/dx (spectral)
    call diffx(u_s,as)
    !add r*(-du/dx) to dr
    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))

    !$OMP END SINGLE
    !$OMP WORKSHARE
    dr(zi:zf,yi:yf,xi:xf) = dr(zi:zf,yi:yf,xi:xf) &
      - r(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !calculate dv/dy (spectral)
    call diffy(v_s,as)
    !add r*(-du/dx) to dr
    !$OMP SINGLE
    call perform_backwards_3dfft(current_state, as, ap(zi:zf,yi:yf,xi:xf))

    !$OMP END SINGLE

    !$OMP WORKSHARE
    dr(zi:zf,yi:yf,xi:xf) = dr(zi:zf,yi:yf,xi:xf) &
      - r(zi:zf,yi:yf,xi:xf)*ap(zi:zf,yi:yf,xi:xf)
    !$OMP END WORKSHARE

    !$OMP END PARALLEL


    !interpolate back to parcels


    call grid2par(current_state,current_state%dp,current_state%parcels%dpdt)
    call grid2par(current_state,current_state%dq,current_state%parcels%dqdt)
    call grid2par(current_state,current_state%dr,current_state%parcels%drdt)

    if ( mod(iteration,current_state%rksteps) == 0) then
      !We now want to determine the maximum vorticity
      omax = maxval(dp**2 + dq**2 + dr**2)
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

      if (current_state%dtmax .gt. dtmax) then
        current_state%dtm = dtmax
      endif

      !print *, "vorticity tendency"
      !print *, "dtmax=",dtmax
      !print *, maxval(dp), maxval(dq), maxval(dr)
      !print *, minval(dp), minval(dq), minval(dr)
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
