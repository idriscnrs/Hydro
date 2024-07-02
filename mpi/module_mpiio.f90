module hydro_mpiio
#ifdef MPIIO
  use hydro_utils
  implicit none

contains

subroutine output_mpiio
  use hydro_mpi
  use hydro_commons
  use hydro_parameters
  implicit none

  integer :: fh, nout, sz
  integer,parameter :: ndim=3
  integer,dimension(ndim) :: glob_pos, glob_uold_shape, pos_noghost, uold_shape, uold_noghost_shape
  integer :: type_in, type_out
  integer(kind=MPI_OFFSET_KIND)::disp
  character(LEN=5)  :: charoutnum
  character(LEN=80) :: filename

  if (rang==0) print *,'MPI-IO: Outputting array of size=',ntx,nty,nvar

  call timer_start()

  !Determine filename
  nout=nstep/noutput
  call title(nout,charoutnum)
  filename='output_'//trim(charoutnum)//'.mp'

  !Open file
  call MPI_File_open(comm2d,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,code)

  !Write run parameters
  if (rang==0) then
    call MPI_File_write(fh,t    ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,code)
    call MPI_File_write(fh,gamma,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,code)
    call MPI_File_write(fh,ntx  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_write(fh,nty  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_write(fh,nvar ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_write(fh,nstep,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_get_position(fh,disp,code)
  end if

  !Broadcast offset
  call MPI_Sizeof(disp,sz,code)
  call MPI_Bcast(disp,sz,MPI_BYTE,0,comm2d,code)

  !Create subarray to represent data to write inside each process without the ghost cells
  uold_shape=shape(uold)
  uold_noghost_shape=(/nx,ny,nvar/)
  pos_noghost=(/2,2,0/)
  call MPI_Type_create_subarray(ndim,uold_shape,uold_noghost_shape,pos_noghost,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_in,code)
  call MPI_Type_commit(type_in,code)

  !Create subarray to represent data to write inside the global array
  glob_uold_shape=(/ntx,nty,nvar/)
  glob_pos=(/startx,starty,0/)
  call MPI_Type_create_subarray(ndim,glob_uold_shape,uold_noghost_shape,glob_pos,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_out,code)
  call MPI_Type_commit(type_out,code)

  !Set file view (each process sees only its own part of uold)
  call MPI_File_set_view(fh,disp,MPI_DOUBLE_PRECISION,type_out,"native",MPI_INFO_NULL,code)

  !Write data
  call MPI_File_write_all(fh,uold,1,type_in,MPI_STATUS_IGNORE,code)

  !Close file
  call MPI_File_close(fh,code)

  call timer_stop("MPI-IO",.true.,int(ntx,kind=8)*nty*nvar*8)
end subroutine output_mpiio


subroutine input_mpiio
  use hydro_mpi
  use hydro_commons
  use hydro_parameters
  implicit none

  integer :: fh, sz
  integer,parameter :: ndim=3
  integer,dimension(ndim) :: glob_pos, glob_uold_shape, pos_noghost, uold_shape, uold_noghost_shape
  integer :: type_in, type_out
  integer(kind=MPI_OFFSET_KIND)::disp
  integer :: ntx_prov,nty_prov,nvar_prov
  real(kind=prec_real) :: gamma_prov

  if (rang==0) print *,'MPI-IO: Reading file ',trim(restart_file)

  call timer_start()

  !Open file
  call MPI_File_open(MPI_COMM_WORLD,trim(restart_file),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,code)

  !Read run parameters
!  if (rang==0) then
    call MPI_File_read_all(fh,t         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,code)
    call MPI_File_read_all(fh,gamma_prov,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,code)
    call MPI_File_read_all(fh,ntx_prov  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_read_all(fh,nty_prov  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_read_all(fh,nvar_prov ,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)
    call MPI_File_read_all(fh,nstep,1,MPI_INTEGER,MPI_STATUS_IGNORE,code)

    call MPI_File_get_position(fh,disp,code)

    if(gamma_prov/=gamma .or. ntx_prov/=ntx .or. nty_prov/=nty .or. nvar_prov/=nvar)then
      print *,'ERROR: Invalid parameters in restart_file:'
      print *,rang,'In nml:          gamma=',gamma,' nx=',ntx,' ny=',nty,'nvar=',nvar
      print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',ntx_prov,' ny=',nty_prov,'nvar=',nvar_prov
      call MPI_Abort(comm2d,1,code)
    endif
!  end if

!  !Broadcast offset
!  call MPI_Sizeof(disp,sz,code)
!  call MPI_Bcast(disp,sz,MPI_BYTE,0,comm2d,code)

  !Create subarray to represent data to read inside each process without the ghost cells
  uold_shape=shape(uold)
  uold_noghost_shape=(/nx,ny,nvar/)
  pos_noghost=(/2,2,0/)
  call MPI_Type_create_subarray(ndim,uold_shape,uold_noghost_shape,pos_noghost,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_in,code)
  call MPI_Type_commit(type_in,code)

  !Create subarray to represent data to read inside the global array
  glob_uold_shape=(/ntx,nty,nvar/)
  glob_pos=(/startx,starty,0/)
  call MPI_Type_create_subarray(ndim,glob_uold_shape,uold_noghost_shape,glob_pos,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_out,code)
  call MPI_Type_commit(type_out,code)

  !Set file view (each process sees only its own part of uold)
  call MPI_File_set_view(fh,disp,MPI_DOUBLE_PRECISION,type_out,"native",MPI_INFO_NULL,code)

  !Read data
  call MPI_File_read_all(fh,uold,1,type_in,MPI_STATUS_IGNORE,code)

  !Close file
  call MPI_File_close(fh,code)

  call timer_stop("MPI-IO",.false.,int(ntx,kind=8)*nty*nvar*8)
end subroutine input_mpiio

#endif
end module hydro_mpiio
