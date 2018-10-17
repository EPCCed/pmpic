!component that obtains the velocity from the voticity via
! laplacian(A) = vort, where velocity = -curl(A)
module vort2vel_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use pencil_fft_mod, only : initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  use MPI
  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords, par2grid, cache_parcel_interp_weights, grid2par
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  use fftops_mod, only: fftops_init, diffx, diffy, diffz, laplinv, spectral_filter
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION) :: PI
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: a, b, c, d, e, f
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: kx, ky, kz
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: k2
  integer :: fourier_space_sizes(3)
  integer :: ierr
  integer :: handle
  integer :: nx, ny, nz
  real(kind=DEFAULT_PRECISION) :: dz, hdzi
  logical :: k2eq0
  !type(halo_communication_type), save :: halo_swap_state
  integer :: iteration

  public vort2vel_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The fft solver component descriptor
  type(component_descriptor_type) function vort2vel_get_descriptor()
    vort2vel_get_descriptor%name="vort2vel"
    vort2vel_get_descriptor%version=0.1
    vort2vel_get_descriptor%initialisation=>initialisation_callback
    vort2vel_get_descriptor%timestep=>timestep_callback
    vort2vel_get_descriptor%finalisation=>finalisation_callback
  end function vort2vel_get_descriptor

  !> This initialisation callback sets up the pencil fft module, allocates data for the fourier space variables
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), ALLOCATABLE :: kx_global, ky_global, kz_global
    real(kind=DEFAULT_PRECISION) :: kzmax, kymax, kxmax

    integer :: i, my_y_start, my_x_start, rank, j, k

    iteration=0

    PI=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

    fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)

    if (my_x_start .eq. 1 .and. my_y_start .eq. 1) then
      k2eq0 = .true.
    else
      k2eq0 = .false.
    endif

    !initialise spectral derivatives module
    call fftops_init(current_state,my_x_start,my_y_start,fourier_space_sizes)

    nx=fourier_space_sizes(X_INDEX)
    ny=fourier_space_sizes(Y_INDEX)
    nz=fourier_space_sizes(Z_INDEX)

    dz=current_state%global_grid%resolution(Z_INDEX)
    hdzi = 1./2./dz

    allocate(a(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(b(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(c(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(d(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(e(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(f(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

    allocate(current_state%u_s%data(nz,ny,nx))
    allocate(current_state%v_s%data(nz,ny,nx))
    allocate(current_state%w_s%data(nz,ny,nx))



    call register_routine_for_timing("vort2vel",handle,current_state)

  end subroutine initialisation_callback







  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i,j,k, xi, xf, yi, yf, zi, zf
    real(kind=DEFAULT_PRECISION) :: L, pi
    real(kind=DEFAULT_PRECISION), allocatable, dimension(:,:) :: atop, abot, btop,bbot
    real(kind=DEFAULT_PRECISION), ALLOCATABLE, dimension(:) :: ubar, vbar
    real(kind=DEFAULT_PRECISION) :: uavg, vavg
    real(kind=DEFAULT_PRECISION) :: umax, dtmax, dtmaxglobal, vmax,wmax, Lx, Ly, Lz

    allocate(atop(ny,nx), abot(ny,nx), btop(ny,nx), bbot(ny,nx))
    if (k2eq0) then
      allocate(ubar(nz), vbar(nz))
    endif

    !print *, "Entering vort2vel"
    call timer_start(handle)


    call cache_parcel_interp_weights(current_state)

    call par2grid(current_state,current_state%parcels%p,current_state%p)
    call par2grid(current_state,current_state%parcels%q,current_state%q)
    call par2grid(current_state,current_state%parcels%r,current_state%r)

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

    !Take fft of vorticities to get them into semi-spectral space



    ! get fft of p and put it in a
    call perform_forward_3dfft(current_state, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), a)

    ! get fft of q and put it in b
    call perform_forward_3dfft(current_state, current_state%q%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
        start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), b)

    ! get fft of r and put it in c
    call perform_forward_3dfft(current_state, current_state%r%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
        start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), c)



    ! we now want to correct vorticity so div(vort) = 0
    ! To do this we introduce a scalar lambda where vort = vort_p - grad(lambda)
    ! (this doesn't affect the velocity as curl(gradient) = 0)
    ! Thus div(vort) = div(vort_p) - laplacian(lambda) = 0, hence laplacian(lambda) = div(vort)

    !$OMP PARALLEL default(shared)

    ! First we calculate div(vort_p)
    call diffx(a,d) !d = da/dx
    call diffy(b,e) !e = db/dy

    !f = dc/dz
    !$OMP WORKSHARE
    f(1,:,:) = hdzi * ( 4.*c(2,:,:) - c(3,:,:) - 3.*c(1,:,:) )
    f(2:nz-1,:,:) = (c(3:nz,:,:) - c(1:nz-2,:,:))*hdzi
    f(nz,:,:) = hdzi * ( c(nz-2,:,:) + 3.*c(nz,:,:) - 4.*c(nz-1,:,:) )


    !f -> d + e + f    (f -> div(vort))
    f(:,:,:) = d(:,:,:) + e(:,:,:) + f(:,:,:)
    !$OMP END WORKSHARE

    !$OMP SINGLE
    if (k2eq0) f(:,1, 1) = 0.0 !set constant part of f to zero
    !$OMP END SINGLE

    ! invert laplacian (f = div(vort) -> f = lambda)
    call laplinv(f,df_zero_on_boundary=.true.)

    !spectrally filter lambda
    call spectral_filter(f)

    ! now we want to correct vorticity: vort -> vort - grad(lambda)

    ! d = df/dx
    call diffx(f,d)
    !$OMP WORKSHARE
    a(:,:,:) = a(:,:,:) - d(:,:,:)
    !$OMP END WORKSHARE

    !d = df/dy
    call diffy(f,d)
    !$OMP WORKSHARE
    b(:,:,:) = b(:,:,:) - d(:,:,:)
    !$OMP END WORKSHARE

    !d = df/dz (df/fz=0 on boundaries)
    call diffz(f,d)
    !$OMP WORKSHARE
    c(:,:,:) = c(:,:,:) - d(:,:,:)
    !$OMP END WORKSHARE

    !ensure vertical average of vorticity is zero
    !$OMP WORKSHARE
    c(1,:,:) = 0.
    !$OMP END WORKSHARE

    !spectrally filter a, b and c
    call spectral_filter(a,out=d)
    call spectral_filter(b,out=e)
    call spectral_filter(c,out=f)

    !cache top and bottom of a and b for use with z derivatives of A and B potentials later
    !$OMP WORKSHARE
    atop(:,:) = a(nz,:,:)
    abot(:,:) = a(1,:,:)
    btop(:,:) = b(nz,:,:)
    bbot(:,:) = b(1,:,:)
    !$OMP END WORKSHARE

    !$OMP END PARALLEL

    !inverse fft corrected and filtered vorticities to physical space
    call perform_backwards_3dfft(current_state, d, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))
    call perform_backwards_3dfft(current_state, e, current_state%q%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))
    call perform_backwards_3dfft(current_state, f, current_state%r%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))



    !if we are the process at the bottom left hand corner then calculate the mean velocity
    if (k2eq0) then
      ubar(1)=0
      vbar(1)=0
      do k=1,nz-1
        ubar(k+1) = ubar(k) + dz/2. * (e(k,1,1)+e(k+1,1,1))
        vbar(k+1) = vbar(k) - dz/2. * (d(k,1,1)+d(k+1,1,1))
      enddo
      uavg=1./nz * sum(ubar(1:nz-1)+0.5*ubar(nz))
      vavg=1./nz * sum(vbar(1:nz-1)+0.5*vbar(nz))
      ubar(:) = ubar(:) - uavg
      vbar(:) = vbar(:) - vavg
    endif


    !$OMP PARALLEL default(shared)

    !invert vorticity (not filtered!!) to get potentials a, b and c
    call laplinv(a, f_zero_on_boundary=.true.)
    call laplinv(b, f_zero_on_boundary=.true.)
    call laplinv(c, df_zero_on_boundary=.true.)

    ! Now we want to compute the velocity

    ! u = db/dz - dc/dy

    ! e=db/dz, d=dc/dy
    call diffz(b,e,bbot,btop)
    call diffy(c,d)

    !$OMP WORKSHARE
    f(:,:,:) = e(:,:,:) - d(:,:,:)
    !$OMP END WORKSHARE

    !$OMP SINGLE
    if (k2eq0) then
      f(:,1,1) = ubar(:)
      !these steps should be unneccessary but let's put them in to be safe
    !  f(:,1,2) = 0.
  !    f(:,2,1) = 0.
  !    f(:,2,2) = 0.
    endif
    !$OMP END SINGLE

    call spectral_filter(f, out=current_state%u_s%data)

    !$OMP END PARALLEL


    call perform_backwards_3dfft(current_state, f, current_state%u%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))



    !$OMP PARALLEL default(shared)

    ! v = dc/dx - da/dz

    ! e=db/dz, d=dc/dy
    call diffz(a,e,abot,atop)
    call diffx(c,d)
    !$OMP WORKSHARE
    f(:,:,:) =  d(:,:,:) - e(:,:,:)
    !$OMP END WORKSHARE

    !$OMP SINGLE
    if (k2eq0) then
      f(:,1,1) = vbar(:)
      !these steps should be unneccessary but let's put them in to be safe
    !  f(:,1,2) = 0.
  !    f(:,2,1) = 0.
!      f(:,2,2) = 0.
    endif
    !$OMP END SINGLE

    call spectral_filter(f, out=current_state%v_s%data)

    !$OMP END PARALLEL


    call perform_backwards_3dfft(current_state, f, current_state%v%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))


    !$OMP parallel default(shared)
    ! w = da/dy - db/dx

    ! d=db/dx, e=da/dy
    call diffx(b,d)
    call diffy(a,e)

    !$OMP WORKSHARE
    f(:,:,:) = e(:,:,:) - d(:,:,:)
    !$OMP END WORKSHARE

    call spectral_filter(f, out=current_state%w_s%data)

    !$OMP END PARALLEL

    call perform_backwards_3dfft(current_state, f, current_state%w%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))


    !interpolate velocity to parcels
    call grid2par(current_state, current_state%u, current_state%parcels%dxdt)
    call grid2par(current_state, current_state%v, current_state%parcels%dydt)
    call grid2par(current_state, current_state%w, current_state%parcels%dzdt)

    if (mod(iteration,current_state%rksteps) ==0 ) then
      !determine the maximum velocity in each direction
      umax = maxval(current_state%u%data(zi:zf,yi:yf,xi:xf)**2)
      umax = sqrt(umax)

      vmax = maxval(current_state%v%data(zi:zf,yi:yf,xi:xf)**2)
      vmax = sqrt(vmax)

      wmax = maxval(current_state%w%data(zi:zf,yi:yf,xi:xf)**2)
      wmax = sqrt(wmax)

      !if the velocoties are zero, make them small (but non-zero) to avoid div zero errors
      if (umax .eq. 0) umax = 1.E-15
      if (vmax .eq. 0) vmax = 1.E-15
      if (wmax .eq. 0) wmax = 1.E-15

      !get the length of our local computational domain
      Lx = current_state%local_grid%size(3)*current_state%global_grid%resolution(3)
      Ly = current_state%local_grid%size(2)*current_state%global_grid%resolution(2)
      Lz = current_state%local_grid%size(1)*current_state%global_grid%resolution(1)

      !We don't want the parcels to move more than one computational domain per timestep
      !so we determine the maximum allowed timestep (i.e. the crossing time for a domain)
      ! This is because when doing a parcel haloswap we can only swap to adjacent processes,
      ! so we don't want a parcel to be able to move further than a single process-length per timestep
      dtmax = minval ((/ Lx/umax, Ly/vmax, Lz/wmax /))


      !This is the local maximum. We want the global maximum so we do a MPI reduction operation
      ! (We want to find the smallert maximum, so we want the minimum of all the maximums)
      call MPI_Allreduce(dtmax,&
                         dtmaxglobal,&
                         1,&
                         PRECISION_TYPE,&
                         MPI_MIN,&
                         current_state%parallel%monc_communicator,&
                         ierr)


      if (current_state%dtmax .gt. dtmaxglobal) then
        current_state%dtm = dtmaxglobal
      else
        current_state%dtm = current_state%dtmax
      endif

      ! print *, "Velocities"
      ! print *, "dtmax=",dtmax
      !print *, maxval(current_state%u%data), maxval(current_state%v%data), maxval(current_state%w%data)
      !print *, minval(current_state%u%data), minval(current_state%v%data), minval(current_state%w%data)

      !print *, "Leaving vort2vel"

    endif


    call timer_stop(handle)

    deallocate(atop, abot, btop, bbot)
    if (k2eq0) then
      deallocate(ubar, vbar)
    endif




    !call MPI_Finalize(ierr)
    !stop

    iteration=iteration+1


  end subroutine timestep_callback




  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_pencil_fft(current_state%parallel%monc_communicator)
    deallocate(a, b, c, d, e, f)
    deallocate(current_state%u_s%data, current_state%v_s%data, current_state%w_s%data)

  end subroutine finalisation_callback



end module
