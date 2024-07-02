module hydro_adios
#ifdef ADIOS
  use hydro_commons
  use hydro_mpi
  use hydro_parameters
  use hydro_utils
  use adios_write_mod
  implicit none

  integer,private :: ierr
  logical,private :: is_initialized

contains

subroutine output_adios
  implicit none

  integer :: nout
  character(LEN=5)  :: charoutnum
  character(LEN=80) :: filename

  integer :: sz_double, sz_int
  integer(kind=8) :: fh, group_size, total_size

  if (rang==0) print *,'ADIOS: Outputting array of size=',ntx,nty,nvar

  call check_adios_write_initialization

  call timer_start()

  !Determine filename
  nout=nstep/noutput
  call title(nout,charoutnum)
  filename='output_'//trim(charoutnum)//'.ad'

  !Create file
  call adios_open(fh,"output",filename,"w",comm2d,ierr)

  !Provide group_size to ADIOS
  call MPI_Type_size(MPI_INTEGER,sz_int,ierr)
  call MPI_Type_size(MPI_DOUBLE_PRECISION,sz_double,ierr)
  group_size = (2+nx*ny*nvar)*int(sz_double,kind=8)+8*int(sz_int,kind=8)
  call adios_group_size(fh,group_size,total_size,ierr)

  !Write run parameters
  call adios_write(fh,"t",    t,ierr)
  call adios_write(fh,"gamma",gamma,ierr)
  call adios_write(fh,"ntx",  ntx,ierr)
  call adios_write(fh,"nty",  nty,ierr)
  call adios_write(fh,"nvar", nvar,ierr)
  call adios_write(fh,"nstep",nstep,ierr)

  !Write data structure (necessary to provide to ADIOS all necessary information to write u)
  call adios_write(fh,"nx",nx,ierr)
  call adios_write(fh,"ny",ny,ierr)
  call adios_write(fh,"startx",startx,ierr)
  call adios_write(fh,"starty",starty,ierr)

  !Write data
  call adios_write(fh,"u",uold(3:nx+2,3:ny+2,1:nvar),ierr)

  !Close file
  call adios_close(fh,ierr)

  call timer_stop("ADIOS",.true.,int(ntx,kind=8)*nty*nvar*8)

end subroutine output_adios


subroutine input_adios
  !This implementation works only if restarting with the same number of processes
  implicit none

  integer(kind=8) :: fh, sz_double, sz_int
  integer :: ntx_prov,nty_prov,nvar_prov,sz_prov
  real(kind=prec_real) :: gamma_prov

  call check_adios_write_initialization

  if (rang==0) print *,'ADIOS: Reading file ',trim(restart_file)

  call timer_start()

  !Open file
  call adios_open(fh,"output",trim(restart_file),"r",comm2d,ierr)

  call MPI_Type_size(MPI_DOUBLE_PRECISION,sz_prov,ierr)
  sz_double = sz_prov
  call MPI_Type_size(MPI_INTEGER,sz_prov,ierr)
  sz_int = sz_prov

  !Read run parameters
  call adios_read(fh,"t",    t,sz_double,ierr)
  call adios_read(fh,"gamma",gamma_prov,sz_double,ierr)
  call adios_read(fh,"ntx",  ntx_prov,sz_int,ierr)
  call adios_read(fh,"nty",  nty_prov,sz_int,ierr)
  call adios_read(fh,"nvar", nvar_prov,sz_int,ierr)
  call adios_read(fh,"nstep",nstep,sz_int,ierr)

  !Read data (get data into a contiguous region)
  call adios_read(fh,"u",uold(3:nx+2,3:ny+2,1:nvar),sz_double*nx*ny*nvar,ierr)

  !Close file
  !Warning: data is only available after this
  call adios_close(fh,ierr)

  if(gamma_prov/=gamma .or. ntx_prov/=ntx .or. nty_prov/=nty .or. nvar_prov/=nvar)then
    print *,'ERROR: Invalid parameters in restart_file:'
    print *,rang,'In nml:          gamma=',gamma,' nx=',ntx,' ny=',nty,'nvar=',nvar
    print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',ntx_prov,' ny=',nty_prov,'nvar=',nvar_prov
    call MPI_Abort(comm2d,1,code)
  endif

  call timer_stop("ADIOS",.false.,int(ntx,kind=8)*nty*nvar*8)
end subroutine input_adios


subroutine check_adios_write_initialization
  implicit none

  if(is_initialized) return

  call adios_init("config_adios.xml",comm2d,ierr)
  is_initialized = .true.

end subroutine check_adios_write_initialization



subroutine clean_adios_write
  implicit none

   if(.not.is_initialized) return

   call adios_finalize(rang,ierr)
   is_initialized = .false.

end subroutine clean_adios_write

#endif
end module hydro_adios
