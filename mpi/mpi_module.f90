module hydro_mpi
  use mpi
  use hydro_precision

  implicit none

  integer, parameter :: NORTH=1
  integer, parameter :: SOUTH=2
  integer, parameter :: WEST=3
  integer, parameter :: EAST=4

  integer, parameter :: NDIMS=2

  integer, parameter                  :: etiquette=100

  integer :: ntx !Global nx
  integer :: nty !Global ny
  integer :: startx !Global x-position
  integer :: starty !Global y-position

  integer                             :: nb_procs, rang, code
  integer                             :: comm2D    ! 2D cartesian topology communicator
  integer                             :: bloc_dim1 ! Derived datatypes to send boundary blocks on dimension 1 (X)
  integer                             :: bloc_dim2 ! Derived datatypes to send boundary blocks on dimension 2 (Y)
  logical                             :: reorganize
  integer, dimension(NDIMS)           :: coords, dims
  logical, dimension(NDIMS)           :: periods
  integer, dimension(2*NDIMS)         :: voisins


! Ajout RD communication RDMA contigus
  !real(kind=prec_real),allocatable,dimension(:,:,:) :: tabx1,tabx2,tabx3,tabx4
  !real(kind=prec_real),allocatable,dimension(:,:,:) :: taby1,taby2,taby3,taby4
  real(kind=prec_real),allocatable,dimension(:,:,:) :: tabx1,tabx2,tabx3,tabx4
  real(kind=prec_real),allocatable,dimension(:,:,:) :: taby1,taby2,taby3,taby4

contains

subroutine init_mpi
  use hydro_parameters

  integer :: type_temp, size_double
  integer(kind=MPI_ADDRESS_KIND) :: stride
  integer :: alloc_err

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

  ! Create 2D cartesian topology
  dims(:)=0
  call MPI_DIMS_CREATE(nb_procs,ndims,dims,code)
  periods(:)=.false.
  if(boundary_left==3 .AND. boundary_right==3) periods(1)=.true.
  if(boundary_up  ==3 .AND. boundary_down ==3) periods(2)=.true.
  reorganize=.true.
  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods,reorganize,comm2d,code)
  call MPI_COMM_RANK(comm2d,rang,code)
  call MPI_CART_COORDS(comm2d,rang,ndims,coords,code)

  if (rang==0) then
     print *
     print *,' MPI Execution with ',nb_procs,' processes (',dims(1),'x',dims(2),')'
     print *,' Starting time integration, nx = ',nx,' ny = ',ny
     print *
  endif

  ! Calcul de la taille des sous-domaines locaux
  ntx = nx
  nty = ny

  startx=coords(1)*(ntx/dims(1))+min(coords(1),mod(ntx,dims(1)))
  if (coords(1)<mod(nx,dims(1))) then
     nx=nx/dims(1)+1
  else
     nx=nx/dims(1)
  endif

  starty=coords(2)*(nty/dims(2))+min(coords(2),mod(nty,dims(2)))
  if (coords(2)<mod(ny,dims(2))) then
     ny=ny/dims(2)+1
  else
     ny=ny/dims(2)
  endif

  ! Modify boundary conditions and determine neighbors list
  ! 0 stands for MPI boundary conditions
  voisins(:)=MPI_PROC_NULL
  call MPI_CART_SHIFT(comm2d,0,1,voisins(WEST), voisins(EAST), code)
  call MPI_CART_SHIFT(comm2d,1,1,voisins(SOUTH),voisins(NORTH),code)
  if(voisins(NORTH)/=MPI_PROC_NULL) boundary_up   =0
  if(voisins(SOUTH)/=MPI_PROC_NULL) boundary_down =0
  if(voisins(WEST) /=MPI_PROC_NULL) boundary_left =0
  if(voisins(EAST) /=MPI_PROC_NULL) boundary_right=0

  ! Create derived datatypes
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_double,code)
  stride = (nx+4)*(ny+4)*size_double

  call MPI_TYPE_VECTOR(ny,2,nx+4,MPI_DOUBLE_PRECISION,type_temp,code)
  call MPI_TYPE_CREATE_HVECTOR(nvar,1,stride,type_temp,bloc_dim1,code)
  call MPI_TYPE_COMMIT(bloc_dim1,code)

  call MPI_TYPE_VECTOR(2,nx,nx+4,MPI_DOUBLE_PRECISION,type_temp,code)
  call MPI_TYPE_CREATE_HVECTOR(nvar,1,stride,type_temp,bloc_dim2,code)
  call MPI_TYPE_COMMIT(bloc_dim2,code)

! ajout RD communication contigus RDMA
!  allocate(tabx1(2*nvar*ny),tabx2(2*nvar*ny),tabx3(2*nvar*ny),tabx4(2*nvar*ny))
!  allocate(taby1(2*nvar*nx),taby2(2*nvar*nx),taby3(2*nvar*nx),taby4(2*nvar*nx))

  allocate(tabx1(1:2,3:ny+2,1:nvar),tabx2(1:2,3:ny+2,1:nvar),stat=alloc_err)
  allocate(tabx3(1:2,3:ny+2,1:nvar),tabx4(1:2,3:ny+2,1:nvar),stat=alloc_err)
  allocate(taby1(3:nx+2,1:2,1:nvar),taby2(3:nx+2,1:2,1:nvar),stat=alloc_err)
  allocate(taby3(3:nx+2,1:2,1:nvar),taby4(3:nx+2,1:2,1:nvar),stat=alloc_err)
  if (alloc_err .ne. 0) then
   print *, "Allocation error in tabxy subroutine init_mpi"
   STOP
  endif
  !$acc enter data create(tabx1(1:2,3:ny+2,1:nvar),tabx2(1:2,3:ny+2,1:nvar), &
  !$acc                   tabx3(1:2,3:ny+2,1:nvar),tabx4(1:2,3:ny+2,1:nvar), &
  !$acc                   taby1(3:nx+2,1:2,1:nvar),taby2(3:nx+2,1:2,1:nvar), &
  !$acc                   taby3(3:nx+2,1:2,1:nvar),taby4(3:nx+2,1:2,1:nvar)  )
end subroutine init_mpi

subroutine clean_mpi

! ajout RD deallocate temporary varaibles for sendrecv with compiler standard
! action
  !$acc exit data delete(tabx1,tabx2,tabx3,tabx4,taby1,taby2,taby3,taby4)
  if (allocated(tabx1)) deallocate(tabx1) 
  if (allocated(tabx2)) deallocate(tabx2)
  if (allocated(tabx3)) deallocate(tabx3)
  if (allocated(tabx4)) deallocate(tabx4)
  if (allocated(taby1)) deallocate(taby1)
  if (allocated(taby2)) deallocate(taby2)
  if (allocated(taby3)) deallocate(taby3)
  if (allocated(taby4)) deallocate(taby4)

  call MPI_TYPE_FREE(bloc_dim1,code)
  call MPI_TYPE_FREE(bloc_dim2,code)
  call MPI_COMM_FREE(comm2d,code)

end subroutine clean_mpi

#ifdef key_gpu
subroutine initialisation_openacc
  use openacc
  character(len=6) :: local_rank_env
  integer          :: local_rank_env_status, local_rank
 
  ! Initialisation d'OpenACC
  !$acc init
 
  ! Récupération du rang local du processus via la variable d'environnement
  ! positionnée par Slurm, l'utilisation de MPI_Comm_rank n'étant pas encore
  ! possible puisque cette routine est utilisée AVANT l'initialisation de MPI
  call get_environment_variable(name="SLURM_LOCALID", value=local_rank_env, &
                                status=local_rank_env_status)
 
  if (local_rank_env_status == 0) then
      read(local_rank_env, *) local_rank
      ! Définition du GPU à utiliser via OpenACC
      call acc_set_device_num(local_rank, acc_get_device_type())
  else
      print *, "Erreur : impossible de déterminer le rang local du processus"
      stop 1
  end if
end subroutine initialisation_openacc
#endif

end module hydro_mpi
