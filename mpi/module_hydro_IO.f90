!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module_hydro_IO.f90 --- 
!!!!
!! subroutine read_params
!! subroutine output 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hydro_IO
  use hydro_commons
  use hydro_utils
  use hydro_adios
  use hydro_hdf5
  use hydro_mpiio
  use hydro_netcdf
  use hydro_pnetcdf
  use hydro_sionlib
  implicit none

  logical :: restart = .false.

contains

subroutine read_params
  use hydro_parameters
  implicit none

  ! Local variables
  integer(kind=prec_int) :: narg,iargc
  character(LEN=256) :: infile

  ! Namelists
  namelist/run/nstepmax,tend,noutput,on_output, &
               fileformat_in,fileformat_out,restart,restart_file
  namelist/mesh/nx,ny,dx,boundary_left,boundary_right,boundary_down,boundary_up, &
                idimbloc,jdimbloc
  namelist/hydro/gamma,courant_factor,smallr,smallc,niter_riemann, &
                 iorder,scheme,slope_type

  narg = iargc()
  IF(narg .NE. 1)THEN
     write(*,*)'You should type: a.out input.nml'
     write(*,*)'File input.nml should contain a parameter namelist'
     STOP
  END IF
  CALL getarg(1,infile)
  open(1,file=infile)
  read(1,NML=run)
  read(1,NML=mesh)
  read(1,NML=hydro)
  close(1)

end subroutine read_params


subroutine output
  use hydro_mpi,only:rang
  implicit none

  integer :: ierr

  select case(fileformat_out)
    case('std')
      call output_std
#ifdef ADIOS
    case('adios')
      call output_adios
#endif
#ifdef HDF5
    case('hdf5')
      call output_hdf5
#endif
#ifdef MPIIO
    case('mpiio')
      call output_mpiio
#endif
#ifdef NETCDF4
    case('netcdf')
      call output_netcdf
#endif
#ifdef PNETCDF
    case('pnetcdf')
      call output_pnetcdf
#endif
#ifdef SIONLIB
    case('sionlib')
      call output_sionlib
#endif
    case default
      if (rang==0) print *,'ERROR: Invalid fileformat_out in namelist: ',fileformat_out
      call MPI_Finalize(ierr)
      stop
  end select
end subroutine output


subroutine output_std
  use hydro_mpi
  use hydro_commons
  use hydro_parameters
  implicit none

  ! Local variables
  character(LEN=80) :: filename
  character(LEN=5)  :: char,charpe
  integer(kind=prec_int) :: nout,MYPE=0

  if (rang==0) print *,'Std: Outputting array of size=',nx,ny,nvar

  call timer_start()

  nout=nstep/noutput
  call title(nout,char)
  call title(rang,charpe)
  filename='output_'//TRIM(char)//'.'//TRIM(charpe)
  open(10,file=filename,form='unformatted')
  rewind(10)
  write(10)real(t,kind=prec_output),real(gamma,kind=prec_output)
  write(10)nx,ny,nvar,nstep
  write(10)real(uold(imin+2:imax-2,jmin+2:jmax-2,1:nvar),kind=prec_output)
  close(10)

  call timer_stop("Std",.true.,int(ntx,kind=8)*nty*nvar*8)

end subroutine output_std


subroutine input
  use hydro_commons,only:nstep,nstepstart
  use hydro_mpi,only:rang
  implicit none

  integer :: ierr

  select case(fileformat_in)
    case('std')
      call input_std
#ifdef ADIOS
    case('adios')
      call input_adios
#endif
#ifdef HDF5
    case('hdf5')
      call input_hdf5
#endif
#ifdef MPIIO
    case('mpiio')
      call input_mpiio
#endif
#ifdef NETCDF4
    case('netcdf')
      call input_netcdf
#endif
#ifdef PNETCDF
    case('pnetcdf')
      call input_pnetcdf
#endif
#ifdef SIONLIB
    case('sionlib')
      call input_sionlib
#endif
    case default
      if (rang==0) print *,'ERROR: Invalid fileformat_in in namelist: ',fileformat_in
      call MPI_Finalize(ierr)
      stop
  end select

  nstepstart=nstep
end subroutine input


subroutine input_std
  use hydro_mpi
  use hydro_commons
  use hydro_parameters
  use hydro_utils,only:title
  implicit none

  ! Local variables
  character(LEN=80) :: filename
  character(LEN=5)  :: charpe
  integer :: nx_prov,ny_prov,nvar_prov
  real(kind=prec_real) :: gamma_prov

  if (rang==0) print *,'Std: Reading file ',trim(restart_file)

  call timer_start()

  call title(rang,charpe)
  filename=trim(restart_file)//'.'//TRIM(charpe)
  open(10,file=filename,form='unformatted')
  read(10) t,gamma_prov
  read(10) nx_prov,ny_prov,nvar_prov,nstep

  if(gamma_prov/=gamma .or. nx_prov/=nx .or. ny_prov/=ny .or. nvar_prov/=nvar)then
    close(10)
    print *,'ERROR: Invalid parameters in restart_file:'
    print *,rang,'In nml:          gamma=',gamma,' nx=',nx,' ny=',ny,'nvar=',nvar
    print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',nx_prov,' ny=',ny_prov,'nvar=',nvar_prov
    call MPI_Finalize(code)
    stop
  endif

  read(10) uold(imin+2:imax-2,jmin+2:jmax-2,1:nvar)
  close(10)

  call timer_stop("Std",.false.,int(ntx,kind=8)*nty*nvar*8)
end subroutine input_std


subroutine print_io_versions
  use hydro_mpi,only:rang
  implicit none

  integer::ierr
  integer::majnum, minnum, relnum
  integer::main_version, sub_version, patch_level, fileformat_version
  character(len=2)::main_versionc, sub_versionc, patch_levelc, fileformat_versionc
  integer::mpiversion,mpisubversion
  character(len=4)::majnumc,minnumc,relnumc
  character(len=2)::mpiversionc,mpisubversionc

  if (rang==0) then
#ifdef ADIOS
     print *,'  Compiled     with ADIOS support'
#else
     print *,'  Not compiled with ADIOS support'
#endif
#ifdef HDF5
     call H5get_libversion_f(majnum,minnum,relnum,ierr)
     write(majnumc,'(I3)') majnum;write(minnumc,'(I3)') minnum;write(relnumc,'(I3)') relnum
     print *,'  Compiled     with parallel HDF5 support (version ', &
             trim(adjustl(majnumc)),'.',trim(adjustl(minnumc)),'.',trim(adjustl(relnumc)),')'
#else
     print *,'  Not compiled with parallel HDF5 support'
#endif
     call MPI_GET_VERSION(mpiversion,mpisubversion,ierr)
     write(mpiversionc,'(I2)') mpiversion;write(mpisubversionc,'(I2)') mpisubversion
#ifdef MPIIO
     print *,'  Compiled     with MPI-I/O support (MPI-',&
             trim(adjustl(mpiversionc)),'.',trim(adjustl(mpisubversionc)),')'
#else
     print *,'  Not compiled with MPI-I/O support (MPI-',&
            trim(adjustl(mpiversionc)),'.',trim(adjustl(mpisubversionc)),')'
#endif
#ifdef NETCDF4
     print *,'  Compiled     with NetCDF-4 support (',trim(nf90_inq_libvers()),')'
#else
     print *,'  Not compiled with NetCDF-4 support'
#endif
#ifdef PNETCDF
     print *,'  Compiled     with Parallel-NetCDF support (',trim(nfmpi_inq_libvers()),')'
#else
     print *,'  Not compiled with Parallel-NetCDF support'
#endif
#ifdef SIONLIB
     !Not (yet) implemented in the library
     !call fsion_get_version(main_version,sub_version,patch_level,fileformat_version)
     print *,'  Compiled     with SIONlib support'
#else
     print *,'  Not compiled with SIONlib support'
#endif
     print *,''
  endif

end subroutine print_io_versions

end module hydro_IO
