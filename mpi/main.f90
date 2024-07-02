!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! main.f90 --- 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program hydro_main
  use hydro_commons
  use hydro_parameters
  use hydro_IO
  use hydro_mpi
  use hydro_principal
  use Ptim
  implicit none

  real(kind=prec_real)   :: dt, tps_start_iter, tps_stop_iter

  !--MPI--!
#if defined(key_gpu)
  ! Initialization of OpenACC environment
  call initialisation_openacc()
#endif
  ! Initialization of MPI environment
  call MPI_INIT(code)
  ! MPI timing
  call PTIM_start(LABEL="Hydro MPI2D")
  !--MPI--!

  ! Read run parameters
  call read_params

  ! Initialize MPI domains and communicator
  call init_mpi
  call print_io_versions

  ! Initialize hydro grid
  call init_hydro

  ! Main time loop
  do while (t < tend .and. nstep < nstepmax)

     ! Start timing iteration
     if (rang==0) tps_start_iter = MPI_WTIME()

     ! Output results
     if( on_output .and. nstep>nstepstart .and. MOD(nstep,noutput)==0)then
        !$acc update host(uold(imin+2:imax-2,jmin+2:jmax-2,1:nvar))
        call output
     end if

     ! Compute new time-step
     if(MOD(nstep,2)==0)then
        call cmpdt(dt)
        if(nstep==0)dt=dt/2._prec_real
     endif

     ! Directional splitting
     if(MOD(nstep,2)==0)then
        call godunov(1,dt)
        call godunov(2,dt)
     else
        call godunov(2,dt)
        call godunov(1,dt)
     end if

     nstep=nstep+1
     t=t+dt
     !--MPI--!
      if (rang==0) then
        tps_stop_iter = MPI_WTIME()
        write(*,'("step=",I6," t=",1pe10.3," dt=",1pe10.3," (",0pf8.3," s.)")') &
             nstep,t,dt,tps_stop_iter-tps_start_iter
     endif
      !--MPI--!

  end do

  ! Final output
  if (on_output) then
   !$acc update host(uold(imin+2:imax-2,jmin+2:jmax-2,1:nvar))
   call output
  endif

  !$acc exit data delete(uold(:,:,:))
  !--MPI--!
  ! Free MPI structures
  call clean_mpi
#ifdef HDF5
  call clean_hdf5
#endif
#ifdef ADIOS
  call clean_adios_write
#endif
  ! MPI timing result
  call PTIM_stop(LABEL="Hydro MPI2D")
  ! Shutting down MPI environment
  call MPI_FINALIZE(code)
  !--MPI--!
end program hydro_main
