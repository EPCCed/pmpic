!component that tests the tridiagonal laplacian inversion method
module laplinv_test_mod
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
  use fftops_mod, only: fftops_init, diffx, diffy, diffz, laplinv
  implicit none

#ifndef TEST_MODE
  private
#endif

  !real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: cos_x, cos_y
  real(kind=DEFAULT_PRECISION) :: PI
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: p_s, q_s, r_s, a, b, c, aref, bref, cref
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: kx, ky
  real(kind=DEFAULT_PRECISION), dimension(:,:), ALLOCATABLE :: k2
  integer :: fourier_space_sizes(3)
  integer :: ierr
  integer :: handle
  real(kind=DEFAULT_PRECISION) :: p1, p2, p3, q1, q2, q3, alpha, beta, gamma, lz
  !type(halo_communication_type), save :: halo_swap_state

  public laplinv_test_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The fft solver component descriptor
  type(component_descriptor_type) function laplinv_test_get_descriptor()
    laplinv_test_get_descriptor%name="laplinv_test"
    laplinv_test_get_descriptor%version=0.1
    laplinv_test_get_descriptor%initialisation=>initialisation_callback
    laplinv_test_get_descriptor%timestep=>timestep_callback
    laplinv_test_get_descriptor%finalisation=>finalisation_callback
  end function laplinv_test_get_descriptor

  !> This initialisation callback sets up the pencil fft module, allocates data for the fourier space variables
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION), dimension(:), ALLOCATABLE :: kx_global, ky_global, kz_global
    real(kind=DEFAULT_PRECISION) :: kzmax, kymax, kxmax

    integer :: i, my_y_start, my_x_start, rank, j, k
    integer :: istart, jstart

    print *, "laplinv test init"

    PI=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

    lz=2.*pi

    ! We're going to modify/override the default domain sizes to fit the analytical form of z=[0:2*pi]
    current_state%global_grid%resolution(1) = lz/(current_state%global_grid%size(1)-1)
    current_state%global_grid%resolution(2) = lz/current_state%global_grid%size(2)
    current_state%global_grid%resolution(3) = lz/current_state%global_grid%size(3)

    fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)

    !initialise spectral derivatives module
    call fftops_init(current_state,my_x_start,my_y_start,fourier_space_sizes)

    !allocate spectral arrays
    allocate(p_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(q_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(r_s(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    ! allocate(a(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    ! allocate(b(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    ! allocate(c(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(aref(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(bref(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    allocate(cref(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

    !allocate wavenumbers
    allocate(kx(fourier_space_sizes(X_INDEX)))
    allocate(ky(fourier_space_sizes(Y_INDEX)))
    allocate(k2(fourier_space_sizes(Y_INDEX),fourier_space_sizes(X_INDEX)))

    !set parameters


    alpha=-1.d0
    beta=2.d0
    p1=0.2*lz
    p2=0.7*lz
    p3=lz
    q1=0.3*lz
    q2=0.8*lz
    q3=lz

    nx=fourier_space_sizes(X_INDEX)
    ny=fourier_space_sizes(Y_INDEX)
    nz=fourier_space_sizes(Z_INDEX)

    z_coords=z_coords/z_coords(nz)
    z_coords=z_coords*lz

    !setup k arrays
    do i=1,nx
      kx(i) = ((my_x_start+i-2)/2) /(lz)
      do j=1,ny
        if (i .eq. 1) ky(j) = ((my_y_start+j-2)/2) / (lz)
        k2(j,i) = kx(i)*kx(i) + ky(j)*ky(j)
      enddo
    enddo

    ! print *, "laplinv ks"
    ! do i=1,10
    !   print *, i, kx(i), ky(i)
    ! enddo

    p_s=0.
    q_s=0.
    r_s=0.
    Aref=0.
    Bref=0.
    Cref=0.



    ! setup fourier space arrays

    do i=1,nx
      do j=1,ny
        if ( (mod(my_x_start+i-1,2) .ne. 0) .and. (mod(my_y_start+i-1,2) .ne. 0)) then
          !real parts
          Aref(:,j,i) = alpha*z_coords(:)*(z_coords(:)-p1)*(z_coords(:)-p2)*(z_coords(:)-p3)
          Bref(:,j,i) = beta*z_coords(:)*(z_coords(:)-q1)*(z_coords(:)-q2)*(z_coords(:)-q3)

          p_s(:,j,i) = 2*alpha*(6*z_coords(:)*z_coords(:) - 3*(p1+p2+p3)*z_coords(:) &
                       + p1*p2 + p1*p3 + p2*p3) - k2(j,i)*Aref(:,j,i)
          q_s(:,j,i) = 2*beta*(6*z_coords(:)*z_coords(:) - 3*(q1+q2+q3)*z_coords(:) &
                       + q1*q2 + q1*q3 + q2*q3) - k2(j,i)*Bref(:,j,i)
        else if ( (mod(my_x_start+i-1,2) .eq. 0) .and. (mod(my_y_start+i-1,2) .eq. 0)) then
          !imaginary parts
          Cref(:,j,i) = -kx(i)*alpha*z_coords(:)**2 &
                        * ( 0.2*z_coords(:)**3 - 0.25*(p1+p2+p3)*z_coords(:)**2 &
                        +  1./3.*(p1*p2 + p2*p3 + p1*p3)*z_coords(:) - 0.5*p1*p2*p3) &
                        -ky(j)*beta*z_coords(:)**2 &
                        * ( 0.2*z_coords(:)**3 - 0.25*(q1+q2+q3)*z_coords(:)**2 &
                        +  1./3.*(q1*q2 + q2*q3 + q1*q3)*z_coords(:) - 0.5*q1*q2*q3) + 20.
          r_s(:,j,i) = -kx(i)*alpha &
                        * ( 4.*z_coords(:)**3 - 3.*(p1+p2+p3)*z_coords(:)**2 &
                        +  2.*(p1*p2 + p2*p3 + p1*p3)*z_coords(:) - p1*p2*p3) &
                        -ky(j)*beta &
                        * ( 4.*z_coords(:)**3 - 3.*(q1+q2+q3)*z_coords(:)**2 &
                        +  2.*(q1*q2 + q2*q3 + q1*q3)*z_coords(:) - q1*q2*q3) &
                        - k2(j,i)*Cref(:,j,i)
        endif



      enddo
    enddo

    !write vorticities to file as they will be overwritten during solving the laplacian
    open(unit=10,file="vort.txt")
    do i=1,nz
      write(10,*) z_coords(i), p_s(i,5,3), q_s(i,5,3), r_s(i,6,4)
    enddo
    close(10)

  end subroutine initialisation_callback







  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i,j,k
    real(kind=DEFAULT_PRECISION) :: L, pi

    print *, "laplinv_test"


    call laplinv(p_s, f_zero_on_boundary=.true.)
    call laplinv(q_s, f_zero_on_boundary=.true.)
    call laplinv(r_s, df_zero_on_boundary=.true.)

    !write computed A and B (stored in p and q) as well as the reference (analytical) values
    open(unit=10,file="pot.dat")
    do i=1,size(p_s,1)
      write(10,*) z_coords(i), p_s(i,5,3), Aref(i,5,3), q_s(i,5,3), Bref(i,5,3), r_s(i,6,4), Cref(i,6,4)
    enddo
    close(10)



     call MPI_Finalize(ierr)
     stop



  end subroutine timestep_callback






  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_pencil_fft(current_state%parallel%monc_communicator)


  end subroutine finalisation_callback



end module
