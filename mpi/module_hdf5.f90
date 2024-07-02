module hydro_hdf5
#ifdef HDF5
  use h5lt
  use hdf5
  use hydro_commons
  use hydro_mpi
  use hydro_parameters
  use hydro_utils,only:title
  implicit none

  integer,private :: status
  logical,private :: is_initialized
  integer(hid_t),private :: mode_id, mode_indep_id
  integer(hid_t),private :: fh,prp_id
  integer(hid_t),private :: da_id,ds_id,fs_id

contains

subroutine output_hdf5
  use hydro_commons
  use hydro_mpi,only:rang
  use hydro_utils
  implicit none

  integer :: nout
  character(LEN=5)  :: charoutnum
  character(LEN=80) :: filename
  integer(size_t) :: sz
  integer(hsize_t),dimension(3) :: counts,dims,starts

  call check_hdf5_initialization

  if (rang==0) print *,'HDF5: Outputting array of size=',ntx,nty,nvar

  call timer_start()

  nout=nstep/noutput
  call title(nout,charoutnum)
  filename='output_'//trim(charoutnum)//'.h5'

  call H5Pcreate_f(H5P_FILE_ACCESS_F,prp_id,status)
  call H5Pset_fapl_mpio_f(prp_id,comm2d,MPI_INFO_NULL,status)
  call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,fh,status,access_prp=prp_id)
  call H5Pclose_f(prp_id,status)

  !Create property list for access mode (used when switching between
  !independent and collective MPI-I/O modes
  call H5Pcreate_f(H5P_DATASET_XFER_F,mode_id,status)
!  call H5Pcreate_f(H5P_DATASET_XFER_F,mode_indep_id,status)

  !Put mode_id in collective mode
  call H5Pset_dxpl_mpio_f(mode_id,H5FD_MPIO_COLLECTIVE_F,status)
!  !Put mode_indep_id in independent mode
!  call H5Pset_dxpl_mpio_f(mode_indep_id,H5FD_MPIO_INDEPENDENT_F,status)

  !Write run parameters
  sz=1
  call H5LTset_attribute_double_f(fh,'/','t',    (/t/),sz,    status)
  call H5LTset_attribute_double_f(fh,'/','gamma',(/gamma/),sz,status)
  call H5LTset_attribute_int_f(fh,'/','ntx',  (/ntx/),  sz, status)
  call H5LTset_attribute_int_f(fh,'/','nty',  (/nty/),  sz, status)
  call H5LTset_attribute_int_f(fh,'/','nvar', (/nvar/), sz, status)
  call H5LTset_attribute_int_f(fh,'/','nstep',(/nstep/),sz, status)

  dims(1)   = ntx ;    dims(2)   = nty ;    dims(3)   = nvar
  counts(1) = nx ;     counts(2) = ny ;     counts(3) = nvar
  starts(1) = startx ; starts(2) = starty ; starts(3) = 0

  !Description of data in memory
  call H5Screate_simple_f(3,counts,ds_id,status)

  !Description of data in file
  call H5Screate_simple_f(3,dims,fs_id,status)
  call H5Dcreate_f(fh,'u',H5T_NATIVE_DOUBLE,fs_id,da_id,status)
  call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)

  call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,uold(3:nx+2,3:ny+2,1:nvar),counts,status,&
                  file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)

  call H5Dclose_f(da_id,status)
  call H5Sclose_f(ds_id,status)
  call H5Sclose_f(fs_id,status)

  call H5Pclose_f(mode_id,status)
!  call H5Pclose_f(mode_indep_id,status)
  call H5Fclose_f(fh,status)

  call timer_stop("HDF5",.true.,int(ntx,kind=8)*nty*nvar*8)

end subroutine output_hdf5



subroutine input_hdf5
  use hydro_commons
  use hydro_mpi,only:rang
  use hydro_utils
  implicit none

  integer :: ntx_prov,nty_prov,nvar_prov
  real(kind=prec_real) :: gamma_prov
  integer,dimension(1) :: ti
  real(kind=8),dimension(1) :: td
  integer(hsize_t),dimension(3) :: counts,starts

  call check_hdf5_initialization

  if (rang==0) print *,'HDF5: Reading file ',trim(restart_file)

  call timer_start()

  call H5Pcreate_f(H5P_FILE_ACCESS_F,prp_id,status)
  call H5Pset_fapl_mpio_f(prp_id,comm2d,MPI_INFO_NULL,status)
  call H5Fopen_f(restart_file,H5F_ACC_RDONLY_F,fh,status,access_prp=prp_id)
  call H5Pclose_f(prp_id,status)

  !Create property list for access mode (used when switching between
  !independent and collective MPI-I/O modes
  call H5Pcreate_f(H5P_DATASET_XFER_F,mode_id,status)

  !Put mode_id in collective mode
  call H5Pset_dxpl_mpio_f(mode_id,H5FD_MPIO_COLLECTIVE_F,status)

  call H5LTget_attribute_double_f(fh,'/','t',    td,status);t=td(1)
  call H5LTget_attribute_double_f(fh,'/','gamma',td,status);gamma_prov=td(1)
  call H5LTget_attribute_int_f(fh,'/','ntx',  ti,status);ntx_prov=ti(1)
  call H5LTget_attribute_int_f(fh,'/','nty',  ti,status);nty_prov=ti(1)
  call H5LTget_attribute_int_f(fh,'/','nvar', ti,status);nvar_prov=ti(1)
  call H5LTget_attribute_int_f(fh,'/','nstep',ti,status);nstep=ti(1)
  if(gamma_prov/=gamma .or. ntx_prov/=ntx .or. nty_prov/=nty .or. nvar_prov/=nvar)then
    print *,'ERROR: Invalid parameters in restart_file:'
    print *,rang,'In nml:          gamma=',gamma,' nx=',ntx,' ny=',nty,'nvar=',nvar
    print *,rang,'In restart_file: gamma=',gamma_prov,' nx=',ntx_prov,' ny=',nty_prov,'nvar=',nvar_prov
    call MPI_Abort(comm2d,1,code)
  endif

  counts(1) = nx ;     counts(2) = ny ;     counts(3) = nvar
  starts(1) = startx ; starts(2) = starty ; starts(3) = 0

  !Description of data in memory
  call H5Screate_simple_f(3,counts,ds_id,status)

  !Description of data in file
  call H5Dopen_f(fh,'u',da_id,status)
  call H5Dget_space_f(da_id,fs_id,status)
  call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
  call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,uold(3:nx+2,3:ny+2,1:nvar),counts,status,&
                 file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
  call H5Dclose_f(da_id,status)
  call H5Sclose_f(fs_id,status)
  call H5Sclose_f(ds_id,status)


  call H5Pclose_f(mode_id,status)
  call H5Fclose_f(fh,status)

  call timer_stop("HDF5",.false.,int(ntx,kind=8)*nty*nvar*8)
end subroutine input_hdf5



subroutine check_hdf5_initialization
  implicit none

  if(is_initialized) return

  call h5open_f(status)
  is_initialized = .true.

end subroutine check_hdf5_initialization



subroutine clean_hdf5
  implicit none

   if(.not.is_initialized) return

   call h5close_f(status)
   is_initialized = .false.

end subroutine clean_hdf5

#endif
end module hydro_hdf5
