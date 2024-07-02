module hydro_pnetcdf
#ifdef PNETCDF
  use hydro_commons
  use hydro_mpi
  use hydro_parameters
  use hydro_utils
  use pnetcdf
  implicit none

  integer(kind=mpi_offset_kind),dimension(1),parameter,private::one1d=1

  integer,private :: fh, mode, status

contains

subroutine output_pnetcdf
  implicit none

  integer :: nout
  character(LEN=5)  :: charoutnum
  character(LEN=80) :: filename

  integer :: ntx_dim_id,nty_dim_id,nvar_dim_id
  integer :: gamma_id,nstep_id,ntx_id,nty_id,nvar_id,t_id,u_id
  integer(kind=MPI_OFFSET_KIND),dimension(3) :: starts

  if (rang==0) print *,'Parallel netCDF: Outputting array of size=',ntx,nty,nvar

  call timer_start()

  !Determine filename
  nout=nstep/noutput
  call title(nout,charoutnum)
  filename='output_'//trim(charoutnum)//'.nc'

  !Create file
  mode = NF_64BIT_DATA
  status = nf90mpi_create(comm2d,filename,mode,MPI_INFO_NULL,fh);call CHKERR(status)

  !Define dimensions
  status = nf90mpi_def_dim(fh,'ntx', int(ntx, kind=mpi_offset_kind),ntx_dim_id)
  status = nf90mpi_def_dim(fh,'nty', int(nty, kind=mpi_offset_kind),nty_dim_id)
  status = nf90mpi_def_dim(fh,'nvar',int(nvar,kind=mpi_offset_kind),nvar_dim_id)

  !Define variables
  status = nf90mpi_def_var(fh,'t',    NF90_DOUBLE,t_id)
  status = nf90mpi_def_var(fh,'gamma',NF90_DOUBLE,gamma_id)
  status = nf90mpi_def_var(fh,'ntx',  NF90_INT,   ntx_id)
  status = nf90mpi_def_var(fh,'nty',  NF90_INT,   nty_id)
  status = nf90mpi_def_var(fh,'nvar', NF90_INT,   nvar_id)
  status = nf90mpi_def_var(fh,'nstep',NF90_INT,   nstep_id)

  status = nf90mpi_def_var(fh,'u',NF90_DOUBLE,(/ntx_dim_id,nty_dim_id,nvar_dim_id/),u_id)

  !Definition of dimensions and variables finished
  status = nf90mpi_enddef(fh);call CHKERR(status)

  !Write run parameters
  status = nf90mpi_begin_indep_data(fh);call CHKERR(status)
  if (rang==0) then
     status = nf90mpi_put_var(fh,t_id,    t);    call CHKERR(status)
     status = nf90mpi_put_var(fh,gamma_id,gamma);call CHKERR(status)
     status = nf90mpi_put_var(fh,ntx_id,  ntx);  call CHKERR(status)
     status = nf90mpi_put_var(fh,nty_id,  nty);  call CHKERR(status)
     status = nf90mpi_put_var(fh,nvar_id, nvar); call CHKERR(status)
     status = nf90mpi_put_var(fh,nstep_id,nstep);call CHKERR(status)
  end if

  status = nf90mpi_end_indep_data(fh);call CHKERR(status)

  starts(1) = startx+1 ; starts(2) = starty+1 ; starts(3) = 1
  status = nf90mpi_put_var_all(fh,u_id,uold(3:nx+2,3:ny+2,1:nvar),start=starts);call CHKERR(status)

  status = nf90mpi_close(fh);call CHKERR(status)

  call timer_stop("Parallel netCDF",.true.,int(ntx,kind=8)*nty*nvar*8)

end subroutine output_pnetcdf


subroutine input_pnetcdf
  implicit none

  integer :: ntx_prov,nty_prov,nvar_prov
  real(kind=prec_real) :: gamma_prov
  real(kind=8),dimension(1) :: td
  integer,dimension(1) :: ti
  integer :: var_id
  integer(kind=MPI_OFFSET_KIND),dimension(3) :: starts

  if (rang==0) print *,'Parallel netCDF: Reading file ',trim(restart_file)

  call timer_start()

  mode = NF_NOWRITE
  status = nf90mpi_open(comm2d,trim(restart_file),mode,MPI_INFO_NULL,fh)

  !Read run parameters
  status = nf90mpi_inq_varid(fh,'t',    var_id);status = nf90mpi_get_var_all(fh,var_id,td);t=td(1)
  status = nf90mpi_inq_varid(fh,'gamma',var_id);status = nf90mpi_get_var_all(fh,var_id,td);gamma_prov=td(1)
  status = nf90mpi_inq_varid(fh,'ntx',  var_id);status = nf90mpi_get_var_all(fh,var_id,ti);ntx_prov=ti(1)
  status = nf90mpi_inq_varid(fh,'nty',  var_id);status = nf90mpi_get_var_all(fh,var_id,ti);nty_prov=ti(1)
  status = nf90mpi_inq_varid(fh,'nvar', var_id);status = nf90mpi_get_var_all(fh,var_id,ti);nvar_prov=ti(1)
  status = nf90mpi_inq_varid(fh,'nstep',var_id);status = nf90mpi_get_var_all(fh,var_id,ti);nstep=ti(1)
  if(gamma_prov/=gamma .or. ntx_prov/=ntx .or. nty_prov/=nty .or. nvar_prov/=nvar)then
    print *,'ERROR: Invalid parameters in restart_file:'
    print *,rang,'In nml:          gamma=',gamma,' nx=',ntx,' ny=',nty,'nvar=',nvar
    print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',ntx_prov,' ny=',nty_prov,'nvar=',nvar_prov
    call MPI_Abort(comm2d,1,code)
  endif

  starts(1) = startx+1 ; starts(2) = starty+1 ; starts(3) = 1
  status = nf90mpi_inq_varid(fh,'u',var_id)
  status = nf90mpi_get_var_all(fh,var_id,uold(3:nx+2,3:ny+2,1:nvar),start=starts);call CHKERR(status)

  status = nf90mpi_close(fh);call CHKERR(status)

  call timer_stop("Parallel netCDF",.false.,int(ntx,kind=8)*nty*nvar*8)

end subroutine input_pnetcdf


subroutine CHKERR(status)
  use hydro_mpi
  implicit none

  integer::status

  if ((status)/=NF_NOERR) then
    write(*,*) nf90mpi_strerror(status)
    call MPI_Abort(comm2d,1,code)
  end if
end subroutine CHKERR

#endif
end module hydro_pnetcdf
