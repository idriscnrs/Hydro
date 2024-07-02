module hydro_netcdf
#ifdef NETCDF4
  use hydro_commons
  use hydro_mpi
  use hydro_parameters
  use hydro_utils
  use netcdf
  implicit none

  integer,private :: fh, status
contains

subroutine output_netcdf
  implicit none

  integer :: nout
  character(LEN=5)  :: charoutnum
  character(LEN=80) :: filename

  integer :: old_mode
  integer :: ntx_dim_id,nty_dim_id,nvar_dim_id
  integer :: gamma_id,nstep_id,ntx_id,nty_id,nvar_id,t_id,u_id
  integer,dimension(3) :: dims,starts

  if (rang==0) print *,'NetCDF: Outputting array of size=',ntx,nty,nvar

  call timer_start()

  !Determine filename
  nout=nstep/noutput
  call title(nout,charoutnum)
  filename='output_'//trim(charoutnum)//'.nc4'

  !Create file
  status = nf90_create(trim(filename),IOR(NF90_NETCDF4,NF90_MPIIO),fh,comm=comm2d,info=MPI_INFO_NULL);call CHKERR(status)

  !Define dimensions
  status = nf90_def_dim(fh,"ntx", ntx, ntx_dim_id)
  status = nf90_def_dim(fh,"nty", nty, nty_dim_id)
  status = nf90_def_dim(fh,"nvar",nvar,nvar_dim_id)

  !Define variables
  status = nf90_def_var(fh,"t",    NF90_DOUBLE,t_id)    ;call CHKERR(status)
  status = nf90_def_var(fh,"gamma",NF90_DOUBLE,gamma_id);call CHKERR(status)
  status = nf90_def_var(fh,"ntx",  NF90_INT,   ntx_id)  ;call CHKERR(status)
  status = nf90_def_var(fh,"nty",  NF90_INT,   nty_id)  ;call CHKERR(status)
  status = nf90_def_var(fh,"nvar", NF90_INT,   nvar_id) ;call CHKERR(status)
  status = nf90_def_var(fh,"nstep",NF90_INT,   nstep_id);call CHKERR(status)

  dims(1)=ntx_dim_id
  dims(2)=nty_dim_id
  dims(3)=nvar_dim_id
  status = nf90_def_var(fh,"u",NF90_DOUBLE,dims,u_id);call CHKERR(status)

  !Deactivate prefilling of variables (performance impact)
  status = nf90_set_fill(fh,NF90_NOFILL,old_mode);call CHKERR(status)

  !Definition of dimensions and variables finished
  status = nf90_enddef(fh);call CHKERR(status)

  !Write run parameters
  !By default all accesses are independent
  if (rang==0) then
     status = nf90_put_var(fh,t_id,    t);    call CHKERR(status)
     status = nf90_put_var(fh,gamma_id,gamma);call CHKERR(status)
     status = nf90_put_var(fh,ntx_id,  ntx);  call CHKERR(status)
     status = nf90_put_var(fh,nty_id,  nty);  call CHKERR(status)
     status = nf90_put_var(fh,nvar_id, nvar); call CHKERR(status)
     status = nf90_put_var(fh,nstep_id,nstep);call CHKERR(status)
  end if

  !Set access to collective for this variable
  status = nf90_var_par_access(fh,u_id, NF90_COLLECTIVE);call CHKERR(status)

  !Write data
  starts(1) = startx+1 ; starts(2) = starty+1 ; starts(3) = 1

  status = nf90_put_var(fh,u_id,uold(3:nx+2,3:ny+2,1:nvar),start=starts);call CHKERR(status)

  !Close file
  status = nf90_close(fh);call CHKERR(status)

  call timer_stop("NetCDF",.true.,int(ntx,kind=8)*nty*nvar*8)
end subroutine output_netcdf


subroutine input_netcdf
  implicit none

  integer :: ntx_prov,nty_prov,nvar_prov
  real(kind=prec_real) :: gamma_prov
  real(kind=8),dimension(1) :: td
  integer,dimension(1) :: ti
  integer :: var_id
  integer,dimension(3) :: starts

  if (rang==0) print *,'NetCDF: Reading file ',trim(restart_file)

  call timer_start()

  status =  nf90_open(trim(restart_file),IOR(NF90_NOWRITE,NF90_MPIIO),fh,comm=comm2d,info=MPI_INFO_NULL);call CHKERR(status)

  !Read run parameters
  status = nf90_inq_varid(fh,'t',    var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,t);         call CHKERR(status)
  status = nf90_inq_varid(fh,'gamma',var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,gamma_prov);call CHKERR(status)
  status = nf90_inq_varid(fh,'ntx',  var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,ntx_prov);  call CHKERR(status)
  status = nf90_inq_varid(fh,'nty',  var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,nty_prov);  call CHKERR(status)
  status = nf90_inq_varid(fh,'nvar', var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,nvar_prov); call CHKERR(status)
  status = nf90_inq_varid(fh,'nstep',var_id);call CHKERR(status);status = nf90_get_var(fh,var_id,nstep);     call CHKERR(status)
  if(gamma_prov/=gamma .or. ntx_prov/=ntx .or. nty_prov/=nty .or. nvar_prov/=nvar)then
    print *,'ERROR: Invalid parameters in restart_file:'
    print *,rang,'In nml:          gamma=',gamma,' nx=',ntx,' ny=',nty,'nvar=',nvar
    print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',ntx_prov,' ny=',nty_prov,'nvar=',nvar_prov
    call MPI_Abort(comm2d,1,code)
  endif

  !Get variable information
  status = nf90_inq_varid(fh,'u',var_id);call CHKERR(status)

  !Set access to collective for this variable
  status = nf90_var_par_access(fh,var_id, NF90_COLLECTIVE);call CHKERR(status)

  !Read data
  starts(1) = startx+1 ; starts(2) = starty+1 ; starts(3) = 1
  status = nf90_get_var(fh,var_id,uold(3:nx+2,3:ny+2,1:nvar),start=starts);call CHKERR(status)

  status = nf90_close(fh);call CHKERR(status)

  call timer_stop("NetCDF",.false.,int(ntx,kind=8)*nty*nvar*8)
end subroutine input_netcdf


subroutine CHKERR(status)
  use hydro_mpi
  implicit none

  integer::status

  if ((status)/=NF90_NOERR) then
    write(*,*) nf90_strerror(status)
    call MPI_Abort(comm2d,1,code)
  end if
end subroutine CHKERR

#endif
end module hydro_netcdf
