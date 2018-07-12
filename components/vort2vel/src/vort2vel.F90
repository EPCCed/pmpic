!component that obtains the velocity from the voticity via
! laplacian(A) = vort, where velocity = -curl(A)
module vort2vel_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use pencil_fft_mod, only : initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  !use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  !use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
  !     init_halo_communication, finalise_halo_communication, blocking_halo_swap, get_single_field_per_halo_cell
  !use registry_mod, only : is_component_enabled
  !use logging_mod, only : LOG_ERROR, log_master_log
  !use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE, mpi_wtime
  use MPI
  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords
  use timer_mod, only: register_routine_for_timing, timer_start, timer_stop
  implicit none

#ifndef TEST_MODE
  private
#endif

  !real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: cos_x, cos_y
  real(kind=DEFAULT_PRECISION) :: PI
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: p_s, q_s, r_s, lambda_s
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: kx, ky, kz
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: k2
  integer :: fourier_space_sizes(3)
  integer :: ierr
  integer :: handle
  !type(halo_communication_type), save :: halo_swap_state

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
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION), dimension(:), ALLOCATABLE :: kx_global, ky_global, kz_global
    real(kind=DEFAULT_PRECISION) :: kzmax, kymax, kxmax

    integer :: i, my_y_start, my_x_start, rank, j, k

    PI=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

    fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)

    !call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)

    allocate(p_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(q_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(r_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(lambda_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

    rank=current_state%parallel%my_rank
    !call MPI_Barrier(current_state%parallel%monc_communicator,ierr)
    !print *, rank, "Fourier sizes=", fourier_space_sizes
    !call MPI_Barrier(current_state%parallel%monc_communicator,ierr)
    !print *, rank, "real sizes   =", current_state%local_grid%size
    !call MPI_Barrier(current_state%parallel%monc_communicator,ierr)

    current_state%p%data(:,:,:) = 1.

    !construct k vectors. First we construct them globally, then the local ones

    nx=current_state%global_grid%size(X_INDEX)
    ny=current_state%global_grid%size(Y_INDEX)
    nz=current_state%global_grid%size(Z_INDEX)

    allocate(kx_global( (nx/2+1)*2 ))
    allocate(ky_global( (ny/2+1)*2 ))
    allocate(kz_global(nz))

    do i=0,nx/2
      kx_global(2*i+1)=i/dble(nx)/current_state%global_grid%resolution(X_INDEX)
      kx_global(2*i+2)=i/dble(nx)/current_state%global_grid%resolution(X_INDEX)
    enddo

    do i=0,ny/2
      ky_global(2*i+1)=i/dble(ny)/current_state%global_grid%resolution(Y_INDEX)
      ky_global(2*i+2)=i/dble(ny)/current_state%global_grid%resolution(Y_INDEX)
    enddo

    kz_global(:) = 0.

    allocate(kx(fourier_space_sizes(X_INDEX)))
    allocate(ky(fourier_space_sizes(Y_INDEX)))
    allocate(kz(fourier_space_sizes(Z_INDEX)))

    kx(:) = kx_global(my_x_start:my_x_start+fourier_space_sizes(X_INDEX)-1)
    ky(:) = ky_global(my_y_start:my_y_start+fourier_space_sizes(Y_INDEX)-1)
    kz(:) = kz_global(:)

    kxmax=maxval(kx_global)
    kymax=maxval(ky_global)
    kzmax=maxval(kz_global)

    allocate(k2(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

    deallocate(kx_global, ky_global, kz_global)

    !get k^2 array
    do i=1,fourier_space_sizes(X_INDEX)
      do j=1,fourier_space_sizes(Y_INDEX)
        do k=1,fourier_space_sizes(Z_INDEX)
          k2(k,j,i) = kz(k)*kz(k) + ky(j)*ky(j) + kx(i)*kx(i)
        enddo
      enddo
    enddo

    call register_routine_for_timing("vort2vel",handle,current_state)

  end subroutine initialisation_callback

  !> Timestep call back, which will transform to Fourier space, do a tridiagonal solve and then back into time space
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i
    real(kind=DEFAULT_PRECISION) :: L, pi

    integer :: xi, xf, yi, yf, zi, zf


    call timer_start(handle)


    do i=1,3
      start_loc(i)=current_state%local_grid%local_domain_start_index(i)
      end_loc(i)=current_state%local_grid%local_domain_end_index(i)
    end do

    xi=start_loc(X_INDEX)
    xf=  end_loc(X_INDEX)
    yi=start_loc(y_INDEX)
    yf=  end_loc(y_INDEX)
    zi=start_loc(z_INDEX)
    zf=  end_loc(z_INDEX)

    pi=3.141592

    L = current_state%global_grid%top(X_INDEX)-current_state%global_grid%bottom(X_INDEX)

    !the current aim just takes p, differentiates it in spectral space via wavenumber multiplication,
    ! then returns this (in real space) as q

    ! p = input array, set to sin(2*pi*x/L) where L is the length of the domain in the x direction
    ! q = output from dp/dx
    ! r reference values for dp/dz = 2*pi/L * cos(2*pi*x/L)

    do i=1,size(current_state%p%data,X_INDEX)
      current_state%p%data(:,:,i) = sin(2*pi*x_coords(i)/L)
      current_state%r%data(:,:,i) = 2*pi/L*cos(2*pi*x_coords(i)/L)
    enddo

  !  print *, "Gradient should be (analytically derived)", current_state%r%data(1,3,3:10),current_state%r%data(1,3,xf-5:xf)


    call perform_forward_3dfft(current_state, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), p_s)

    !calculate q=dp/dx in spectral space ( q_s = 2*pi*i*kx * p_s)
    call diffx(p_s,q_s,kx)

    call perform_backwards_3dfft(current_state, q_s, current_state%q%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))

    !print *, "value of gradient (numerically derived)=",current_state%q%data(1,3,3:10), current_state%q%data(1,3,xf-5:xf)

    !output some data for graphing
    open(unit=10,file="fft.dat")
    do i=xi,xf
      !print *, x_coords(i), current_state%r%data(5,5,i), current_state%q%data(5,5,i),&
       !abs(current_state%r%data(5,5,i)-current_state%q%data(5,5,i))
       write(10,*) x_coords(i), current_state%r%data(5,5,i), current_state%q%data(5,5,i),&
        abs(current_state%r%data(5,5,i)-current_state%q%data(5,5,i))
    enddo
    close(10)

   call timer_stop(handle)


    ! call MPI_Finalize(ierr)
    ! stop

    ! call blocking_halo_swap(current_state, halo_swap_state, copy_p_to_halo_buffer, &
    !      perform_local_data_copy_for_p, copy_halo_buffer_to_p)

  end subroutine timestep_callback





  !> Called at MONC finalisation, will call to the pencil fft module to clean itself up and free the pressure term
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_pencil_fft(current_state%parallel%monc_communicator)
    deallocate(p_s, q_s, r_s, lambda_s)

  end subroutine finalisation_callback


  ! calculates i*k_x * in
  subroutine diffx(in, out, kx)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: in
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: kx
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: out
    integer :: oddeven
    integer :: i, j, k, nx, ny, nz
    real(kind=DEFAULT_PRECISION), parameter :: twopi=2.*3.1415926535

    !we have 3 cases...
    ! 1) size(k) is even -No haloswapping needed
    ! 2) size(k) is odd - we need haloswapping
    !     a) first element is the imaginary point - haloswap down
    !     b) first element is the real part - haloswap up

    oddeven=mod(size(kx),2)

    nx=size(in,X_INDEX)
    ny=size(in,Y_INDEX)
    nz=size(in,Z_INDEX)

    select case(oddeven)
    case(0) !we're even

      !print *, kx

      do i=1,nx,2
        do j=1,ny
          do k=1,nz
            out(k,j,i+1) = in(k,j,i)*kx(i)*twopi
            out(k,j,i) = -1.*in(k,j,i+1)*kx(i+1)*twopi
          enddo
        enddo
      enddo




    case(1) !we're odd
      error stop "Not implemented yet"

    end select


  end subroutine

end module
