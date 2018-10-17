! Allows for parcel mixing - both the absorbtion and removal of small parcels and adding in new parcels to fill holes
module parcel_mixing_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE, PARCEL_INTEGER, MPI_PARCEL_INT
  use state_mod, only: model_state_type
  use monc_component_mod, only: component_descriptor_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_integer_array, options_get_real_array
  use parcel_interpolation_mod, only: nx, ny, nz, dx, dy, dz, cache_parcel_interp_weights, is, js, ks, &
                                      delxs, delys, delzs, perform_halo_swap, x_coords, y_coords, z_coords
  use parcel_haloswap_mod, only: parcel_haloswap
  use MPI
  use timer_mod
  use prognostics_mod, only : prognostic_field_type

  implicit none


  integer :: ierr

  integer :: iteration=0
  real(kind=DEFAULT_PRECISION), parameter :: pi=4.*atan(1.d0)

  integer :: handle

  real(kind=DEFAULT_PRECISION) :: vmin

  integer, allocatable, dimension(:) :: keep
  integer, ALLOCATABLE, dimension(:,:,:) :: npercell
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: vres, bres, hres, pres, qres, rres



contains

  type(component_descriptor_type) function parcel_mixing_get_descriptor()
    parcel_mixing_get_descriptor%name="parcel_mixing"
    parcel_mixing_get_descriptor%version=0.1
    parcel_mixing_get_descriptor%initialisation=>initialisation_callback
    parcel_mixing_get_descriptor%timestep=>timestep_callback
    parcel_mixing_get_descriptor%finalisation=>finalisation_callback
  end function parcel_mixing_get_descriptor


  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state

    call register_routine_for_timing("mixing", handle, state)

    allocate(npercell(nz-1,ny-1,nx-1))
    allocate(vres(nz,ny,nx),bres(nz,ny,nx),hres(nz,ny,nx),pres(nz,ny,nx),qres(nz,ny,nx),rres(nz,ny,nx))

    vmin = dx*dy*dz/6./6./6.


  end subroutine



  subroutine timestep_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer(kind=PARCEL_INTEGER) :: n, norig, last, nadd, nremove
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz
    real(kind=DEFAULT_PRECISION) :: w000, w001, w010, w011, w100, w101, w110, w111
    real(kind=DEFAULT_PRECISION) :: v, portion, vold, qold, vnew, qnew

    !only do this if we're on the last step of an rk timestep
    if (mod(iteration,state%rksteps) == state%rksteps-1) then
      call timer_start(handle)

      call cache_parcel_interp_weights(state)

      allocate(keep(state%parcels%numparcels_local))


      norig = state%parcels%numparcels_local

      !! TEST CODE - DO NOT UNCOMMENT IN PRODUCTION RUN
      !!randomly reduce the volume of some parcels (to test removal)
       ! do i=1,10000
       !   do
       !     call RANDOM_NUMBER(v)
       !     n=v*norig
       !     exit
       !     !if (state%parcels%z(n) .gt. 1. .or. state%parcels%z(n) .lt. 4.) exit
       !   enddo
       !   state%parcels%vol(n) = 0.5*vmin
       !   print *, "Reduced volume of parcel", n
       ! enddo
       !vold=sum(state%parcels%vol(1:norig))
       !qold=sum(state%parcels%h(1:norig)*state%parcels%vol(1:norig))
       !print *, "vold=", vold
       !print *, "qold=", qold
      !! END TEST CODE

      !!$OMP PARALLEL default(shared) private(i,j,k,delx, dely, delz,w000, w001, w010, w011, w100, w101, w110, w111, v) &
      !!$OMP& firstprivate(norig)

      !$OMP parallel WORKSHARE
      npercell(:,:,:) = 0.
      state%vol%data(:,:,:) = 0.
      state%b%data(:,:,:) =0.
      state%hg%data(:,:,:) =0.
      state%p%data(:,:,:)=0.
      state%q%data(:,:,:)=0.
      state%r%data(:,:,:)=0.
      vres(:,:,:) = 0.
      bres(:,:,:) = 0.
      hres(:,:,:) = 0.
      pres(:,:,:) = 0.
      qres(:,:,:) = 0.
      rres(:,:,:) = 0.
      !$OMP END parallel WORKSHARE


      nremove=0






      !loop over all parcels, determine those to be removed, and construct the residual and keep grid
      !called in a subroutine so we can have a !$OMP DO REDUCTION() on state%var%data variables
      call construct_grids(state,state%vol%data, state%b%data, state%hg%data, state%p%data, state%q%data, state%r%data, nremove)




      !$OMP parallel WORKSHARE
      state%vol%data(1,:,:) = state%vol%data(1,:,:)*2.
      state%vol%data(nz,:,:) = state%vol%data(nz,:,:)*2
      state%b%data(1,:,:) = state%b%data(1,:,:)*2.
      state%b%data(nz,:,:) = state%b%data(nz,:,:)*2
      state%hg%data(1,:,:) = state%hg%data(1,:,:)*2.
      state%hg%data(nz,:,:) = state%hg%data(nz,:,:)*2
      state%p%data(1,:,:) = state%p%data(1,:,:)*2.
      state%p%data(nz,:,:) = state%p%data(nz,:,:)*2
      state%q%data(1,:,:) = state%q%data(1,:,:)*2.
      state%q%data(nz,:,:) = state%q%data(nz,:,:)*2
      state%r%data(1,:,:) = state%r%data(1,:,:)*2.
      state%r%data(nz,:,:) = state%r%data(nz,:,:)*2

      vres(1,:,:) = vres(1,:,:)*2
      vres(nz,:,:) = vres(nz,:,:)*2
      bres(1,:,:) = bres(1,:,:)*2
      bres(nz,:,:) = bres(nz,:,:)*2
      hres(1,:,:) = hres(1,:,:)*2
      hres(nz,:,:) = hres(nz,:,:)*2
      pres(1,:,:) = pres(1,:,:)*2
      pres(nz,:,:) = pres(nz,:,:)*2
      qres(1,:,:) = qres(1,:,:)*2
      qres(nz,:,:) = qres(nz,:,:)*2
      rres(1,:,:) = rres(1,:,:)*2
      rres(nz,:,:) = rres(nz,:,:)*2
      !$OMP END parallel WORKSHARE







      !haloswap original grids
      call perform_halo_swap(state,state%vol%data,perform_sum=.true.)
      call perform_halo_swap(state,state%b%data,perform_sum=.true.)
      call perform_halo_swap(state,state%hg%data,perform_sum=.true.)
      call perform_halo_swap(state,state%p%data,perform_sum=.true.)
      call perform_halo_swap(state,state%q%data,perform_sum=.true.)
      call perform_halo_swap(state,state%r%data,perform_sum=.true.)



      !print *, "Parcel mixing: number to be removed =      ", nremove


      nadd=0
      n = 1
      !now loop over internal grid cells and see if they have enough parcels. If not, create a new parcel (backfilling)
      do i=3,nx-2
        do j=3,ny-2
          do k=1,nz-1
            if (npercell(k,j,i) < 3) then
              nadd=nadd+1
              !create new parcel
              if (n .le. state%parcels%numparcels_local) then
                do while(keep(n) .ne. 0)
                  !iterate n forward until we find an empty spot to put it in
                  n=n+1
                  if (n .gt. state%parcels%numparcels_local) exit
                enddo
              endif
              !total volume of surrounding grid points
              v = sum(state%vol%data(k:k+1,j:j+1,i:i+1))

              portion = 2*vmin/v !proportion to remove from each neighbour cell

              !set properties of the new parcel
              state%parcels%x(n) = 0.5*(x_coords(i) + x_coords(i+1))
              state%parcels%y(n) = 0.5*(y_coords(j) + y_coords(j+1))
              state%parcels%z(n) = 0.5*(z_coords(k) + z_coords(k+1))

              state%parcels%vol(n) = 2*vmin
              state%parcels%stretch(n) = 0
              ! set the parcel properties to the volume-averaged value for that cell
              state%parcels%b(n) = sum(state%b%data(k:k+1,j:j+1,i:i+1))/v
              state%parcels%h(n) = sum(state%hg%data(k:k+1,j:j+1,i:i+1))/v
              state%parcels%p(n) = sum(state%p%data(k:k+1,j:j+1,i:i+1))/v
              state%parcels%q(n) = sum(state%q%data(k:k+1,j:j+1,i:i+1))/v
              state%parcels%r(n) = sum(state%r%data(k:k+1,j:j+1,i:i+1))/v
              if (n .lt. state%parcels%numparcels_local) then
                keep(n) = 2 !mark this cell as being one to keep, but not as an un-modified cell
              endif
              n=n+1

              if (n .gt. state%parcels%maxparcels_local) then
                print *, "Error: reached maximum number of parcels in rank", state%parallel%my_rank
                error stop "Aborting: maxparcels reached"
              endif

              !adjust residues to take new parcel into account
              vres(k:k+1,j:j+1,i:i+1) = vres(k:k+1,j:j+1,i:i+1) - portion*state%vol%data(k:k+1,j:j+1,i:i+1)
              bres(k:k+1,j:j+1,i:i+1) = bres(k:k+1,j:j+1,i:i+1) - portion*state%b%data(k:k+1,j:j+1,i:i+1)
              hres(k:k+1,j:j+1,i:i+1) = hres(k:k+1,j:j+1,i:i+1) - portion*state%hg%data(k:k+1,j:j+1,i:i+1)
              pres(k:k+1,j:j+1,i:i+1) = pres(k:k+1,j:j+1,i:i+1) - portion*state%p%data(k:k+1,j:j+1,i:i+1)
              qres(k:k+1,j:j+1,i:i+1) = qres(k:k+1,j:j+1,i:i+1) - portion*state%q%data(k:k+1,j:j+1,i:i+1)
              rres(k:k+1,j:j+1,i:i+1) = rres(k:k+1,j:j+1,i:i+1) - portion*state%r%data(k:k+1,j:j+1,i:i+1)

            endif
          enddo
        enddo
      enddo

      !print *, "Parcel mixing: number to be added =        ", nadd


      !we only end out with more than numparcels_local if we added more parcels than we removed
      !if not, we need to backfill later
      state%parcels%numparcels_local = maxval((/ n-1, state%parcels%numparcels_local/))

      !halo swap residual fields
      call perform_halo_swap(state,vres,perform_sum=.true.)
      call perform_halo_swap(state,bres,perform_sum=.true.)
      call perform_halo_swap(state,hres,perform_sum=.true.)
      call perform_halo_swap(state,pres,perform_sum=.true.)
      call perform_halo_swap(state,qres,perform_sum=.true.)
      call perform_halo_swap(state,rres,perform_sum=.true.)


      !now loop over each original parcel and add the residuals

      if (nadd .ne. 0 .or. nremove .ne. 0) then
        !$OMP PARALLEL default(shared) private(i,j,k,delx,dely,delz,w000,w001,w010,w011,w100,w101,w110,w111,v)
        !$OMP DO
        do n=1,norig
          if (keep(n) .eq. 1) then
            !calculate tridiagonl weights
            i=is(n)
            j=js(n)
            k=ks(n)
            delx=delxs(n)
            dely=delys(n)
            delz=delzs(n)
            v = state%parcels%vol(n)
            w000 = (1-delz)*(1-dely)*(1-delx)*v
            w001 = (1-delz)*(1-dely)*(delx)*v
            w010 = (1-delz)*(dely)*(1-delx)*v
            w011 = (1-delz)*(dely)*(delx)*v
            w100 = (delz)*(1-dely)*(1-delx)*v
            w101 = (delz)*(1-dely)*(delx)*v
            w110 = (delz)*(dely)*(1-delx)*v
            w111 = (delz)*(dely)*(delx)*v

            !update parcel values with residual values

            state%parcels%vol(n) = state%parcels%vol(n) &
                    + w000*vres(k,j,i)/state%vol%data(k,j,i) &
                    + w001*vres(k,j,i+1)/state%vol%data(k,j,i+1) &
                    + w010*vres(k,j+1,i)/state%vol%data(k,j+1,i) &
                    + w011*vres(k,j+1,i+1)/state%vol%data(k,j+1,i+1) &
                    + w100*vres(k+1,j,i)/state%vol%data(k+1,j,i) &
                    + w101*vres(k+1,j,i+1)/state%vol%data(k+1,j,i+1) &
                    + w110*vres(k+1,j+1,i)/state%vol%data(k+1,j+1,i) &
                    + w111*vres(k+1,j+1,i+1)/state%vol%data(k+1,j+1,i+1)


            state%parcels%b(n) = state%parcels%b(n)*v &
                    + w000/state%vol%data(k,j,i)*bres(k,j,i) &
                    + w001/state%vol%data(k,j,i+1)*bres(k,j,i+1) &
                    + w010/state%vol%data(k,j+1,i)*bres(k,j+1,i) &
                    + w011/state%vol%data(k,j+1,i+1)*bres(k,j+1,i+1) &
                    + w100/state%vol%data(k+1,j,i)*bres(k+1,j,i) &
                    + w101/state%vol%data(k+1,j,i+1)*bres(k+1,j,i+1) &
                    + w110/state%vol%data(k+1,j+1,i)*bres(k+1,j+1,i) &
                    + w111/state%vol%data(k+1,j+1,i+1)*bres(k+1,j+1,i+1)
            state%parcels%b(n) = state%parcels%b(n)/state%parcels%vol(n)


            state%parcels%h(n) = state%parcels%h(n)*v &
                    + w000/state%vol%data(k,j,i)*hres(k,j,i) &
                    + w001/state%vol%data(k,j,i+1)*hres(k,j,i+1) &
                    + w010/state%vol%data(k,j+1,i)*hres(k,j+1,i) &
                    + w011/state%vol%data(k,j+1,i+1)*hres(k,j+1,i+1) &
                    + w100/state%vol%data(k+1,j,i)*hres(k+1,j,i) &
                    + w101/state%vol%data(k+1,j,i+1)*hres(k+1,j,i+1) &
                    + w110/state%vol%data(k+1,j+1,i)*hres(k+1,j+1,i) &
                    + w111/state%vol%data(k+1,j+1,i+1)*hres(k+1,j+1,i+1)
            state%parcels%h(n) = state%parcels%h(n)/state%parcels%vol(n)


            state%parcels%p(n) = state%parcels%p(n)*v &
                    + w000/state%vol%data(k,j,i)*pres(k,j,i) &
                    + w001/state%vol%data(k,j,i+1)*pres(k,j,i+1) &
                    + w010/state%vol%data(k,j+1,i)*pres(k,j+1,i) &
                    + w011/state%vol%data(k,j+1,i+1)*pres(k,j+1,i+1) &
                    + w100/state%vol%data(k+1,j,i)*pres(k+1,j,i) &
                    + w101/state%vol%data(k+1,j,i+1)*pres(k+1,j,i+1) &
                    + w110/state%vol%data(k+1,j+1,i)*pres(k+1,j+1,i) &
                    + w111/state%vol%data(k+1,j+1,i+1)*pres(k+1,j+1,i+1)
            state%parcels%p(n) = state%parcels%p(n)/state%parcels%vol(n)


            state%parcels%q(n) = state%parcels%q(n)*v &
                    + w000/state%vol%data(k,j,i)*qres(k,j,i) &
                    + w001/state%vol%data(k,j,i+1)*qres(k,j,i+1) &
                    + w010/state%vol%data(k,j+1,i)*qres(k,j+1,i) &
                    + w011/state%vol%data(k,j+1,i+1)*qres(k,j+1,i+1) &
                    + w100/state%vol%data(k+1,j,i)*qres(k+1,j,i) &
                    + w101/state%vol%data(k+1,j,i+1)*qres(k+1,j,i+1) &
                    + w110/state%vol%data(k+1,j+1,i)*qres(k+1,j+1,i) &
                    + w111/state%vol%data(k+1,j+1,i+1)*qres(k+1,j+1,i+1)
            state%parcels%q(n) = state%parcels%q(n)/state%parcels%vol(n)


            state%parcels%r(n) = state%parcels%r(n)*v &
                    + w000/state%vol%data(k,j,i)*rres(k,j,i) &
                    + w001/state%vol%data(k,j,i+1)*rres(k,j,i+1) &
                    + w010/state%vol%data(k,j+1,i)*rres(k,j+1,i) &
                    + w011/state%vol%data(k,j+1,i+1)*rres(k,j+1,i+1) &
                    + w100/state%vol%data(k+1,j,i)*rres(k+1,j,i) &
                    + w101/state%vol%data(k+1,j,i+1)*rres(k+1,j,i+1) &
                    + w110/state%vol%data(k+1,j+1,i)*rres(k+1,j+1,i) &
                    + w111/state%vol%data(k+1,j+1,i+1)*rres(k+1,j+1,i+1)
            state%parcels%r(n) = state%parcels%r(n)/state%parcels%vol(n)
          endif


        enddo
        !$OMP END DO
        !$OMP END PARALLEL
      endif


      !now backfill removed parcels (if necessary)
      if (state%parcels%numparcels_local .eq. norig) then
        !print *, "Backfilling"
        last=norig
        do n=1,norig

          !if (n.ge.last) exit

          !if we wish to remove this parcel...
          if (keep(n) .eq. 0) then
            !go back from the last parcel until we find one that can be moved
            do
              if (keep(last) .eq. 1) exit
              last=last-1
            enddo
            !print *, "backfilling", last, "to", n

            !copy parcel properties

            state%parcels%x(n) = state%parcels%x(last)
            state%parcels%y(n) = state%parcels%y(last)
            state%parcels%z(n) = state%parcels%z(last)

            state%parcels%p(n) = state%parcels%p(last)
            state%parcels%q(n) = state%parcels%q(last)
            state%parcels%r(n) = state%parcels%r(last)

            state%parcels%b(n) = state%parcels%b(last)
            state%parcels%h(n) = state%parcels%h(last)
            state%parcels%vol(n) = state%parcels%vol(last)

            state%parcels%tag(n) = state%parcels%tag(last)
            state%parcels%stretch(n) = state%parcels%stretch(last)

            last=last-1


            !decrease our parcel count
            state%parcels%numparcels_local = state%parcels%numparcels_local-1

          endif

        enddo
      endif

      deallocate(keep)

      n = state%parcels%numparcels_global


      !update global parcel count
      call MPI_Allreduce(state%parcels%numparcels_local,&
                         state%parcels%numparcels_global,&
                         1,&
                         MPI_PARCEL_INT,&
                         MPI_SUM,&
                         state%parallel%monc_communicator,&
                         ierr)




      if ( state%parallel%my_rank .eq. 0 ) then
        print *, "Parcel mixing: net change in parcel number=", state%parcels%numparcels_global-n
      endif

      !!TEST CODE (keep commented for production runs)
      ! n=state%parcels%numparcels_local
      ! vnew=sum(state%parcels%vol(1:n))
      ! qnew=sum(state%parcels%h(1:n)*state%parcels%vol(1:n))
      ! print *, "Final volume", vnew
      ! print *, "Final volume-integrated value", qnew
      !
      ! print *, "vdiff=", vnew-vold, "frac=", abs(vnew-vold)/vold
      ! print *, "qdiff=", qnew-qold, "frac=", abs(qnew-qold)/qold
      ! stop
      !!END TEST CODE

      call timer_stop(handle)

    endif

    iteration = iteration +1


  end subroutine

  subroutine finalisation_callback(state)
    type(model_state_type), intent(inout), target :: state


  end subroutine


  subroutine construct_grids(state,vol, b, h, p, q, r, nremove)
    implicit none

    type(model_state_type), intent(inout), target :: state
    integer(kind=PARCEL_INTEGER), intent(inout) :: nremove
    real(kind=DEFAULT_PRECISION), dimension(:,:,:) :: vol,b,h,p,q,r
    integer(kind=PARCEL_INTEGER) :: n
    integer :: i, j, k
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz
    real(kind=DEFAULT_PRECISION) :: w000, w001, w010, w011, w100, w101, w110, w111
    real(kind=DEFAULT_PRECISION) :: v






      !loop over all parcels, determine those to be removed, and construct the residual and keep grids

      ! where for a quantitiy, x:
      !  xres(i,j,k) =  sum_n x(n)*tridiag_weight*vol(n) (n is range of parcels to be removed)
      !  xkeep(i,j,k) = sum_n x(n)*tridiag_weight*vol(n) (n is range of parcels to keep)

    !$OMP PARALLEL default(shared) private(i,j,k,delx,dely,delz,w000,w001,w010,w011,w100,w101,w110,w111,v)
    !$OMP DO REDUCTION(+:vol, b, h, p, q, r, vres, bres, hres, pres, qres, rres, nremove, npercell)
    do n=1,state%parcels%numparcels_local
      !calculate tridiagonal weights
      i=is(n)
      j=js(n)
      k=ks(n)
      delx=delxs(n)
      dely=delys(n)
      delz=delzs(n)
      v = state%parcels%vol(n)
      w000 = (1-delz)*(1-dely)*(1-delx)*v
      w001 = (1-delz)*(1-dely)*(delx)*v
      w010 = (1-delz)*(dely)*(1-delx)*v
      w011 = (1-delz)*(dely)*(delx)*v
      w100 = (delz)*(1-dely)*(1-delx)*v
      w101 = (delz)*(1-dely)*(delx)*v
      w110 = (delz)*(dely)*(1-delx)*v
      w111 = (delz)*(dely)*(delx)*v

      !determine if we need to remove the parcel
      if (state%parcels%vol(n) .lt. vmin) then
        keep(n) = 0
        nremove=nremove+1

        !interpolate parcel onto residual grids
        !volume
        vres(k,j,i) = vres(k,j,i) + w000
        vres(k,j,i+1) = vres(k,j,i+1) + w001
        vres(k,j+1,i) = vres(k,j+1,i) + w010
        vres(k,j+1,i+1) = vres(k,j+1,i+1) + w011
        vres(k+1,j,i) = vres(k+1,j,i) + w100
        vres(k+1,j,i+1) = vres(k+1,j,i+1) + w101
        vres(k+1,j+1,i) = vres(k+1,j+1,i) + w110
        vres(k+1,j+1,i+1) = vres(k+1,j+1,i+1) + w111

        !buoyancy
        bres(k,j,i) = bres(k,j,i) + w000*state%parcels%b(n)
        bres(k,j,i+1) = bres(k,j,i+1) + w001*state%parcels%b(n)
        bres(k,j+1,i) = bres(k,j+1,i) + w010*state%parcels%b(n)
        bres(k,j+1,i+1) = bres(k,j+1,i+1) + w011*state%parcels%b(n)
        bres(k+1,j,i) = bres(k+1,j,i) + w100*state%parcels%b(n)
        bres(k+1,j,i+1) = bres(k+1,j,i+1) + w101*state%parcels%b(n)
        bres(k+1,j+1,i) = bres(k+1,j+1,i) + w110*state%parcels%b(n)
        bres(k+1,j+1,i+1) = bres(k+1,j+1,i+1) + w111*state%parcels%b(n)

        !humidity
        hres(k,j,i) = hres(k,j,i) + w000*state%parcels%h(n)
        hres(k,j,i+1) = hres(k,j,i+1) + w001*state%parcels%h(n)
        hres(k,j+1,i) = hres(k,j+1,i) + w010*state%parcels%h(n)
        hres(k,j+1,i+1) = hres(k,j+1,i+1) + w011*state%parcels%h(n)
        hres(k+1,j,i) = hres(k+1,j,i) + w100*state%parcels%h(n)
        hres(k+1,j,i+1) = hres(k+1,j,i+1) + w101*state%parcels%h(n)
        hres(k+1,j+1,i) = hres(k+1,j+1,i) + w110*state%parcels%h(n)
        hres(k+1,j+1,i+1) = hres(k+1,j+1,i+1) + w111*state%parcels%h(n)

        !p
        pres(k,j,i) = pres(k,j,i) + w000*state%parcels%p(n)
        pres(k,j,i+1) = pres(k,j,i+1) + w001*state%parcels%p(n)
        pres(k,j+1,i) = pres(k,j+1,i) + w010*state%parcels%p(n)
        pres(k,j+1,i+1) = pres(k,j+1,i+1) + w011*state%parcels%p(n)
        pres(k+1,j,i) = pres(k+1,j,i) + w100*state%parcels%p(n)
        pres(k+1,j,i+1) = pres(k+1,j,i+1) + w101*state%parcels%p(n)
        pres(k+1,j+1,i) = pres(k+1,j+1,i) + w110*state%parcels%p(n)
        pres(k+1,j+1,i+1) = pres(k+1,j+1,i+1) + w111*state%parcels%p(n)

        !q
        qres(k,j,i) = qres(k,j,i) + w000*state%parcels%q(n)
        qres(k,j,i+1) = qres(k,j,i+1) + w001*state%parcels%q(n)
        qres(k,j+1,i) = qres(k,j+1,i) + w010*state%parcels%q(n)
        qres(k,j+1,i+1) = qres(k,j+1,i+1) + w011*state%parcels%q(n)
        qres(k+1,j,i) = qres(k+1,j,i) + w100*state%parcels%q(n)
        qres(k+1,j,i+1) = qres(k+1,j,i+1) + w101*state%parcels%q(n)
        qres(k+1,j+1,i) = qres(k+1,j+1,i) + w110*state%parcels%q(n)
        qres(k+1,j+1,i+1) = qres(k+1,j+1,i+1) + w111*state%parcels%q(n)

        !r
        rres(k,j,i) = rres(k,j,i) + w000*state%parcels%r(n)
        rres(k,j,i+1) = rres(k,j,i+1) + w001*state%parcels%r(n)
        rres(k,j+1,i) = rres(k,j+1,i) + w010*state%parcels%r(n)
        rres(k,j+1,i+1) = rres(k,j+1,i+1) + w011*state%parcels%r(n)
        rres(k+1,j,i) = rres(k+1,j,i) + w100*state%parcels%r(n)
        rres(k+1,j,i+1) = rres(k+1,j,i+1) + w101*state%parcels%r(n)
        rres(k+1,j+1,i) = rres(k+1,j+1,i) + w110*state%parcels%r(n)
        rres(k+1,j+1,i+1) = rres(k+1,j+1,i+1) + w111*state%parcels%r(n)


      else
        keep(n) = 1
        !add this parcel to its cell's count
        npercell(k,j,i) = npercell(k,j,i) + 1

        !interpolate parcel onto grid
        !volume
        vol(k,j,i) = vol(k,j,i) + w000
        vol(k,j,i+1) = vol(k,j,i+1) + w001
        vol(k,j+1,i) = vol(k,j+1,i) + w010
        vol(k,j+1,i+1) = vol(k,j+1,i+1) + w011
        vol(k+1,j,i) = vol(k+1,j,i) + w100
        vol(k+1,j,i+1) = vol(k+1,j,i+1) + w101
        vol(k+1,j+1,i) = vol(k+1,j+1,i) + w110
        vol(k+1,j+1,i+1) = vol(k+1,j+1,i+1) + w111

        !buoyancy
        b(k,j,i) = b(k,j,i) + w000*state%parcels%b(n)
        b(k,j,i+1) = b(k,j,i+1) + w001*state%parcels%b(n)
        b(k,j+1,i) = b(k,j+1,i) + w010*state%parcels%b(n)
        b(k,j+1,i+1) = b(k,j+1,i+1) + w011*state%parcels%b(n)
        b(k+1,j,i) = b(k+1,j,i) + w100*state%parcels%b(n)
        b(k+1,j,i+1) = b(k+1,j,i+1) + w101*state%parcels%b(n)
        b(k+1,j+1,i) = b(k+1,j+1,i) + w110*state%parcels%b(n)
        b(k+1,j+1,i+1) = b(k+1,j+1,i+1) + w111*state%parcels%b(n)

        !humidity
        h(k,j,i) = h(k,j,i) + w000*state%parcels%h(n)
        h(k,j,i+1) = h(k,j,i+1) + w001*state%parcels%h(n)
        h(k,j+1,i) = h(k,j+1,i) + w010*state%parcels%h(n)
        h(k,j+1,i+1) = h(k,j+1,i+1) + w011*state%parcels%h(n)
        h(k+1,j,i) = h(k+1,j,i) + w100*state%parcels%h(n)
        h(k+1,j,i+1) = h(k+1,j,i+1) + w101*state%parcels%h(n)
        h(k+1,j+1,i) = h(k+1,j+1,i) + w110*state%parcels%h(n)
        h(k+1,j+1,i+1) = h(k+1,j+1,i+1) + w111*state%parcels%h(n)

        !p
        p(k,j,i) = p(k,j,i) + w000*state%parcels%p(n)
        p(k,j,i+1) = p(k,j,i+1) + w001*state%parcels%p(n)
        p(k,j+1,i) = p(k,j+1,i) + w010*state%parcels%p(n)
        p(k,j+1,i+1) = p(k,j+1,i+1) + w011*state%parcels%p(n)
        p(k+1,j,i) = p(k+1,j,i) + w100*state%parcels%p(n)
        p(k+1,j,i+1) = p(k+1,j,i+1) + w101*state%parcels%p(n)
        p(k+1,j+1,i) = p(k+1,j+1,i) + w110*state%parcels%p(n)
        p(k+1,j+1,i+1) = p(k+1,j+1,i+1) + w111*state%parcels%p(n)

        !q
        q(k,j,i) = q(k,j,i) + w000*state%parcels%q(n)
        q(k,j,i+1) = q(k,j,i+1) + w001*state%parcels%q(n)
        q(k,j+1,i) = q(k,j+1,i) + w010*state%parcels%q(n)
        q(k,j+1,i+1) = q(k,j+1,i+1) + w011*state%parcels%q(n)
        q(k+1,j,i) = q(k+1,j,i) + w100*state%parcels%q(n)
        q(k+1,j,i+1) = q(k+1,j,i+1) + w101*state%parcels%q(n)
        q(k+1,j+1,i) = q(k+1,j+1,i) + w110*state%parcels%q(n)
        q(k+1,j+1,i+1) = q(k+1,j+1,i+1) + w111*state%parcels%q(n)

        !r
        r(k,j,i) = r(k,j,i) + w000*state%parcels%r(n)
        r(k,j,i+1) = r(k,j,i+1) + w001*state%parcels%r(n)
        r(k,j+1,i) = r(k,j+1,i) + w010*state%parcels%r(n)
        r(k,j+1,i+1) = r(k,j+1,i+1) + w011*state%parcels%r(n)
        r(k+1,j,i) = r(k+1,j,i) + w100*state%parcels%r(n)
        r(k+1,j,i+1) = r(k+1,j,i+1) + w101*state%parcels%r(n)
        r(k+1,j+1,i) = r(k+1,j+1,i) + w110*state%parcels%r(n)
        r(k+1,j+1,i+1) = r(k+1,j+1,i+1) + w111*state%parcels%r(n)

      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine



end module
