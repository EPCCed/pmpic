!module that provides wrapper functions for FFTE calls
!
! NOTE: This implementation is thread safe. Many real-to-complex and complex-to-real
!       FFTs can be called in tandem, one per thread. See pencilfft.F90 for an
!       example of parallel usage 
!
! NOTE: FFTE only has functionality for complex to complex in-place FFTs.
!       This module allows for real to complex FFTs by providing wrappers for
!       the complex to complex calls
!
! example usage: 
!       Takes the r2c of a real array "input", producing complex array "output"
!       Then takes the c2r of "output", returning the real array "input"
!  ________________________________________________________________________________
! |integer :: n                                                                    |
! |double precision :: input(n)                                                    |
! |complex*16 :: output(n/2+1)                                                     |
! |                                                                                |
! |call ffte_init(n) !set up ffte for this problem size                            | 
! |call ffte_r2c(input, output, n) ! real to complex (forward) FFT                 |
! |call ffte_c2r(output,input, n) ! complex to real (reverse) FFT                  |
! |call ffte_finalise() !clean up work arrays                                      |
!  --------------------------------------------------------------------------------
!
module ffte_mod
    use omp_lib 

    implicit none
    
    complex*16, allocatable, dimension(:,:) :: wk ! work array for the FFT calculations
    complex*16, allocatable, dimension(:,:) :: data !the array that will have the in-place FFT applied to it

    contains
    
    !initialises the FFTE routines for the real problem size (n)
    ! NOTE: CALL ME IN AN OPENMP SINGLE REGION!!!!!!!!!!!!!!!!!
    subroutine ffte_init(n)
        integer, intent(in) :: n
        integer :: nthreads
        integer :: i
        
        !Get the maximum number of OpenMP threads that can be used
        nthreads = omp_get_max_threads()
        
        !check that n is a valid size for FFTE to handle
        if (.not. ffte_check_factors(n)) then
            stop "FFTE: FFT size, n, can only contain prime factors of 2, 3 and 5"
        endif
        
        !allocate work and data arrays 
        allocate(wk(2*n,nthreads))
        allocate(data(n,nthreads))
        
        !initialise FFTE (for each thread's work array)
        do i =1,nthreads
            call ZFFT1D(data(:,i),n,0,wk(:,i))
        enddo
    end subroutine
    
    !computes a real-to-complex (e.g. forward) FFT
    ! in  : double precision real array of size n (input)
    ! out : double precision complex array of size n/2+1 (output)
    ! n   : integer - size of in (input)
    !
    ! NOTE: I can be called from within an OpenMP Parallel region
    subroutine ffte_r2c(in, out, n)
        integer, intent(in) :: n
        double precision, intent(in) :: in(n)
        complex*16, intent(out) :: out(n/2+1)

        integer :: i
        integer :: mythread
        
        mythread = omp_get_thread_num()+1

        !copy real input into complex "data" array
        do i=1,n
            data(i,mythread) = dcmplx(in(i), 0.d0)
        enddo
        
        !compute forward FFT (in-place on data)
        call ZFFT1D(data(1,mythread),n,-1,wk(1,mythread))

        !extract the first n/2+1 terms from the FFT and return them as output
        !(The last n/2-1 terms of r2c are complex conjugates of the previous 
        ! ones and so are redundant)
        out(:) = data(1:n/2+1,mythread)

    end subroutine
    
    !computes a complex-to-real (e.g. inverse) FFT
    ! in  : double precision complex array of size n/2+1 (input)
    ! out : double precision real array of size n (output)
    ! n   : integer - size of out (input)
    !
    ! NOTE: I can be called from within an OpenMP Parallel region
    subroutine ffte_c2r(in,out,n)
        integer, intent(in) :: n
        complex*16, intent(in) :: in(n/2+1)
        double precision, intent(out) :: out(n)

        integer :: i, mythread
        
        mythread = omp_get_thread_num()+1

        !construct the array to be inverse FFT'd

        !copy the first n/2+1 terms
        data(1:n/2+1,mythread) = in(:)
        !set the remaining entries to the complex conjugates of the previous ones
        do i=n/2+2,n
            data(i,mythread) = dconjg(in(n-i+2))
        enddo
        
        
        !do the inverse fft (in place on data)
        call ZFFT1D(data(:,mythread),n,1,wk(:,mythread))
        
        !extract the real part of data and place it in out
        do i=1,n
            out(i) = real(data(i,mythread))
        enddo

    end subroutine
    
    !Finalises the FFTE routines
    !NOTE: CALL ME IN AN OMP SINGLE REGION!!!!!
    subroutine ffte_finalise()
        
        !deallocate the work arrays
        deallocate(wk,data)

    end subroutine
    
    !Checks to see if the input, n, only has prime factors of 2, 3 and 5
    ! If it does, return true. Otherwise, return false
    logical function ffte_check_factors(n)
        integer, intent(in) :: n
        integer :: m

        m=n

        do while (mod(m,5) .eq. 0)
            m = m/5
        enddo

        do while (mod(m,3) .eq. 0)
            m = m/3
        enddo

        do while (mod(m,2) .eq. 0)
            m = m/2
        enddo

        if (m .eq. 1) then
            ffte_check_factors= .true.
        else
            ffte_check_factors= .false.
        endif

    end function

end module ffte_mod
