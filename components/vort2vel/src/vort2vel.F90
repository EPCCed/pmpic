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
  use fftops_mod, only: fftops_init, diffx, diffy
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

    !initialise spectral derivatives module
    call fftops_init(current_state,my_x_start,my_y_start,fourier_space_sizes)


    !call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)

    allocate(p_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(q_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(r_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(lambda_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

    rank=current_state%parallel%my_rank
    !specify p and r p=exp(-((x-x0)/r0)^2), r=dp/dx = -2(x-xo)/r0^2*exp(-((x-x0)/r0)^2)
    do i=1,size(current_state%p%data,3)
      current_state%p%data(:,:,i) = exp(-((x_coords(i)-3000.)/1000.)**2)
      current_state%r%data(:,:,i) = exp(-((x_coords(i)-3000.)/1000.)**2)*(-2.)*(x_coords(i)-3000.)/1000./1000.
    enddo

    do j=1,size(current_state%p%data,2)
      current_state%p%data(:,j,:)= current_state%p%data(:,j,:) + exp(-((y_coords(j)-3000.)/1000.)**2)
      current_state%r%data(:,j,:) = exp(-((y_coords(j)-3000.)/1000.)**2)*(-2.)*(y_coords(j)-3000.)/1000./1000.
    enddo

    call register_routine_for_timing("vort2vel",handle,current_state)

  end subroutine initialisation_callback






  !> Timestep call back, which will transform to Fourier space, do a tridiagonal solve and then back into time space
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i,j
    real(kind=DEFAULT_PRECISION) :: L, pi




    call timer_start(handle)

    do i=1,3
      start_loc(i)=current_state%local_grid%local_domain_start_index(i)
      end_loc(i)=current_state%local_grid%local_domain_end_index(i)
    end do





  ! get fft of p

    call perform_forward_3dfft(current_state, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), p_s)

    !calculate q=dp/dx in spectral space ( q_s = 2*pi*i*kx * p_s)
    call diffy(current_state,p_s,q_s)

    ! undo fft of q and put it into q
    call perform_backwards_3dfft(current_state, q_s, current_state%q%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))
    !undo fft of p_s and put it in p
    call perform_backwards_3dfft(current_state, p_s, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
        start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))



    !output some data for graphing/verification
    do j=0,current_state%parallel%processes-1
      call MPI_Barrier(current_state%parallel%monc_communicator,ierr)
      if (current_state%parallel%my_rank .eq. j) then
        if (current_state%parallel%my_rank .eq. 0) then
           open(unit=10,file="fft.dat")
        else
           open(unit=10,file="fft.dat",access="append")
        endif
        do i=3,current_state%local_grid%size(2)+2
          write(10,*) y_coords(i), current_state%p%data(5,i,5), current_state%q%data(5,i,5),&
           current_state%r%data(5,i,5)
        enddo
        close(10)
      endif
    enddo

   call timer_stop(handle)


     call MPI_Finalize(ierr)
     stop

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



end module
