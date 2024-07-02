!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module_hydro_principal.f90 --- 
!!!!
!! subroutine init_hydro
!! subroutine cmpdt(dt)
!! subroutine godunov(idim,dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hydro_principal

contains

subroutine init_hydro
  use hydro_io,only:input,restart
  use hydro_mpi,only:rang
  use hydro_commons
  use hydro_const
  use hydro_parameters
  implicit none

  ! Local variables
  integer(kind=prec_int) :: i,j

  imin=1
  imax=nx+4
  jmin=1
  jmax=ny+4

  allocate(uold(imin:imax,jmin:jmax,1:nvar))
  !$acc enter data create(uold(:,:,:))
  if (restart) then
    call input
    !$acc update device(uold(:,:,:))
    return
  end if

  ! Initial conditions in grid interior
  ! Warning: conservative variables U = (rho, rhou, rhov, E)

!!$  ! Jet
!!$  do j=jmin+2,jmax-2
!!$     do i=imin+2,imax-2
!!$        uold(i,j,ID)=1.0
!!$        uold(i,j,IU)=0.0
!!$        uold(i,j,IV)=0.0
!!$        uold(i,j,IP)=1.0/(gamma-1.0)
!!$     end do
!!$  end do

  ! Wind tunnel with point explosion
  !$acc parallel present(uold(:,:,:))
  !$acc loop gang
  do j=jmin+2,jmax-2
     !$acc loop vector
     do i=imin+2,imax-2
        uold(i,j,ID)=1.0_prec_real
        uold(i,j,IU)=0.0_prec_real
        uold(i,j,IV)=0.0_prec_real
        uold(i,j,IP)=1.d-5
     end do
  end do
  !--MPI--!
  if (rang==0) uold(imin+2,jmin+2,IP)=1./dx/dx
  !--MPI--!
  !$acc end parallel

!!$  ! 1D Sod test
!!$  do j=jmin+2,jmax-2
!!$     do i=imin+2,imax/2
!!$        uold(i,j,ID)=1.0
!!$        uold(i,j,IU)=0.0
!!$        uold(i,j,IV)=0.0
!!$        uold(i,j,IP)=1.0/(gamma-1.0)
!!$     end do
!!$     do i=imax/2+1,imax-2
!!$        uold(i,j,ID)=0.125
!!$        uold(i,j,IU)=0.0
!!$        uold(i,j,IV)=0.0
!!$        uold(i,j,IP)=0.1/(gamma-1.0)
!!$     end do
!!$  end do

end subroutine init_hydro


subroutine cmpdt(dt)
  use hydro_mpi
  use hydro_commons
  use hydro_const
  use hydro_parameters
  use hydro_utils
  implicit none

  ! Dummy arguments
  real(kind=prec_real), intent(out) :: dt  
  ! Local variables
  real(kind=prec_real)   :: dt_local
  integer(kind=prec_int) :: i,j
  real(kind=prec_real)   :: cournox,cournoy,eken
  real(kind=prec_real),  dimension(:,:), allocatable   :: q
  real(kind=prec_real),  dimension(:)  , allocatable   :: e,c

  ! compute time step on grid interior
  cournox = zero
  cournoy = zero

  allocate(q(1:nx,1:IP),e(1:nx),c(1:nx))
  !$acc parallel loop gang num_gangs(gang_number) present(uold(:,:,:)) &
  !$acc private(q(1:nx,1:IP),e(1:nx),c(1:nx)) reduction(max:cournox,cournoy) 
  do j=jmin+2,jmax-2
     !$acc loop vector
     do i=1,nx
        q(i,ID) = max(uold(i+2,j,ID),smallr)
        q(i,IU) = uold(i+2,j,IU)/q(i,ID)
        q(i,IV) = uold(i+2,j,IV)/q(i,ID)
        eken = half*(q(i,IU)**2+q(i,IV)**2)
        q(i,IP) = uold(i+2,j,IP)/q(i,ID) - eken
        e(i)=q(i,IP)
     end do

     call eos(q(1:nx,ID),e(1:nx),q(1:nx,IP),c(1:nx))

     !$acc loop vector reduction(max:cournox,cournoy)
     do i=1,nx
      cournox=max(cournox,c(i)+abs(q(i,IU)))
      cournoy=max(cournoy,c(i)+abs(q(i,IV)))
     enddo

  end do

  deallocate(q,e,c)

  !--MPI--!
  dt_local = courant_factor*dx/max(cournox,cournoy,smallc)
  call MPI_ALLREDUCE(dt_local,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm2d,code)
  !--MPI--!

end subroutine cmpdt


subroutine godunov(idim,dt)
  use hydro_commons
  use hydro_const
  use hydro_parameters
  use hydro_utils
  use hydro_work_space
  implicit none

  ! Dummy arguments
  integer(kind=prec_int), intent(in) :: idim
  real(kind=prec_real),   intent(in) :: dt
  ! Local variables
  integer(kind=prec_int) :: i,j,in
  real(kind=prec_real)   :: dtdx

  ! constant
  dtdx=dt/dx

  !$acc data present(uold(1:nx+4,1:ny+4,1:nvar))
  
  ! Update boundary conditions
  call make_boundary(idim)

  if (idim==1)then
     ! Allocate work space for 1D sweeps
     call allocate_work_space(imin,imax,nx+1)

#if defined (key_riemann_oldalg)
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))      &
 !$acc private (u(imin:imax,1:nvar), q(imin:imax,1:nvar), dq(imin:imax,1:nvar), &
 !$acc          qgdnv(1:nx+1,1:nvar), flux(1:nx+1,1:nvar),                      &
 !$acc          qxm(imin:imax,1:nvar), qxp(imin:imax,1:nvar),n_new(1),          &
 !$acc          qleft(1:nx+1,1:nvar) , qright(1:nx+1,1:nvar), c(imin:imax),     &
 !$acc          rl(1:nx+1), ul(1:nx+1), pl(1:nx+1), cl(1:nx+1), wl(1:nx+1),     &
 !$acc          rr(1:nx+1), ur(1:nx+1), pr(1:nx+1), cr(1:nx+1), wr(1:nx+1),     &
 !$acc          ro(1:nx+1), uo(1:nx+1), po(1:nx+1), co(1:nx+1), wo(1:nx+1),     &
 !$acc          rstar(1:nx+1), ustar(1:nx+1), pstar(1:nx+1), cstar(1:nx+1),     & 
 !$acc          scr(1:nx+1), delp(1:nx+1), pold(1:nx+1), ind(1:nx+1), ind2(1:nx+1), &
 !$acc          sgnm(1:nx+1), spin(1:nx+1), spout(1:nx+1),ushock(1:nx+1), frac(1:nx+1)) 
#elif defined (key_riemann_bigloop)
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))      &
 !$acc private (u(imin:imax,1:nvar), q(imin:imax,1:nvar), dq(imin:imax,1:nvar), &
 !$acc          qgdnv(1:nx+1,1:nvar), flux(1:nx+1,1:nvar),                      &
 !$acc          qxm(imin:imax,1:nvar), qxp(imin:imax,1:nvar),                   &
 !$acc          qleft(1:nx+1,1:nvar) , qright(1:nx+1,1:nvar), c(imin:imax),     &
 !$acc          rl(1:nx+1), ul(1:nx+1), pl(1:nx+1), cl(1:nx+1), wl(1:nx+1),     &
 !$acc          rr(1:nx+1), ur(1:nx+1), pr(1:nx+1), cr(1:nx+1), wr(1:nx+1),     &
 !$acc          ro(1:nx+1), uo(1:nx+1), po(1:nx+1), co(1:nx+1), wo(1:nx+1),     &
 !$acc          rstar(1:nx+1), ustar(1:nx+1), pstar(1:nx+1), cstar(1:nx+1),     & 
 !$acc          scr(1:nx+1), delp(1:nx+1), pold(1:nx+1), ind(1:nx+1),           &
 !$acc          sgnm(1:nx+1), spin(1:nx+1), spout(1:nx+1),ushock(1:nx+1), frac(1:nx+1))
#else
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))      &
 !$acc private (u(imin:imax,1:nvar), q(imin:imax,1:nvar), dq(imin:imax,1:nvar), &
 !$acc          qgdnv(1:nx+1,1:nvar), flux(1:nx+1,1:nvar),                      &
 !$acc          qxm(imin:imax,1:nvar), qxp(imin:imax,1:nvar),                   &
 !$acc          qleft(1:nx+1,1:nvar) , qright(1:nx+1,1:nvar), c(imin:imax)      )      
#endif

     !$acc loop gang
     do j=jmin+2,jmax-2
        ! Gather conservative variables
        !$acc loop vector
        do i=imin,imax
           u(i,ID)=uold(i,j,ID)
           u(i,IU)=uold(i,j,IU)
           u(i,IV)=uold(i,j,IV)
           u(i,IP)=uold(i,j,IP)
        end do

        if(nvar>4)then
           !$acc loop seq
           do in = 5,nvar
              !$acc loop vector
              do i=imin,imax
                 u(i,in)=uold(i,j,in)
              end do
           end do
        end if
 
        ! Convert to primitive variables
        call constoprim(u(:,:),q(:,:),c(:))

        ! Characteristic tracing
        call trace(q(:,:),dq(:,:),c(:),qxm(:,:),qxp(:,:),dtdx)

        !$acc loop seq
        do in = 1,nvar
           !$acc loop vector
           do i=1,nx+1
              qleft (i,in)=qxm(i+1,in)
              qright(i,in)=qxp(i+2,in)
           end do
        end do

        ! Solve Riemann problem at interfaces
#if defined (key_riemann_oldalg)
        call riemann(qleft(:,:),qright(:,:),qgdnv(:,:),   &
                     rl(:),ul(:),pl(:),cl(:),wl(:),       &
                     rr(:),ur(:),pr(:),cr(:),wr(:),       &
                     ro(:),uo(:),po(:),co(:),wo(:),       &
                     rstar(:),ustar(:),pstar(:),cstar(:), &
                     sgnm(:),spin(:),spout(:),ushock(:),  &
                     frac(:),scr(:),delp(:),pold(:),ind(:),ind2(:),n_new(:))
#elif defined (key_riemann_bigloop)
        call riemann_newalg_bigloop(qleft(:,:),qright(:,:),qgdnv(:,:), &
                                    rl(:),ul(:),pl(:),cl(:),wl(:),     &
                                    rr(:),ur(:),pr(:),cr(:),wr(:),     &
                                    ro(:),uo(:),po(:),co(:),wo(:),     &
                                    rstar(:),ustar(:),pstar(:),        &
                                    cstar(:),sgnm(:),spin(:),spout(:), &
                                    ushock(:),frac(:),scr(:),delp(:),  &
                                    pold(:),ind(:)                     )
#else
        call riemann_opt_mem(qleft(:,:),qright(:,:),qgdnv(:,:))
        ! renamned form of riemann_dev_PF2 for explicit purposes
#endif

        ! Compute fluxes
        call cmpflx(qgdnv(:,:),flux(:,:))

        ! Update conservative variables
        !$acc loop vector 
        do i=imin+2,imax-2
           uold(i,j,ID)=u(i,ID)+(flux(i-2,ID)-flux(i-1,ID))*dtdx
           uold(i,j,IU)=u(i,IU)+(flux(i-2,IU)-flux(i-1,IU))*dtdx
           uold(i,j,IV)=u(i,IV)+(flux(i-2,IV)-flux(i-1,IV))*dtdx
           uold(i,j,IP)=u(i,IP)+(flux(i-2,IP)-flux(i-1,IP))*dtdx
        end do

        if(nvar>4)then
           !$acc loop seq
           do in = 5,nvar
             !$acc loop vector
              do i=imin+2,imax-2
                 uold(i,j,in)=u(i,in)+(flux(i-2,in)-flux(i-1,in))*dtdx
              end do
           end do
        end if
     end do
     !$acc end parallel

     ! Deallocate work space
     call deallocate_work_space

  else

     ! Allocate work space for 1D sweeps
     call allocate_work_space(jmin,jmax,ny+1)

#if defined (key_riemann_oldalg)
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))          &
 !$acc private (u(jmin:jmax,1:nvar), qxm(jmin:jmax,1:nvar), qxp(jmin:jmax,1:nvar) , &
 !$acc          q(jmin:jmax,1:nvar)  , dq(jmin:jmax,1:nvar), c(jmin:jmax),          &
 !$acc          qleft(1:ny+1,1:nvar) , qright(1:ny+1,1:nvar),                       &
 !$acc          qgdnv(1:ny+1,1:nvar), flux(1:ny+1,1:nvar), n_new(1),               &
 !$acc          rr(1:ny+1), ur(1:ny+1), pr(1:ny+1), cr(1:ny+1), wr(1:ny+1),         &
 !$acc          ro(1:ny+1), uo(1:ny+1), po(1:ny+1), co(1:ny+1), wo(1:ny+1),         &
 !$acc          rl(1:ny+1), ul(1:ny+1), pl(1:ny+1), cl(1:ny+1), wl(1:ny+1),         &
 !$acc          rstar(1:ny+1), ustar(1:ny+1), pstar(1:ny+1), cstar(1:ny+1),         & 
 !$acc          scr(1:ny+1), delp(1:ny+1), pold(1:ny+1), ind(1:ny+1), ind2(1:ny+1), &
 !$acc          sgnm(1:ny+1), spin(1:ny+1), spout(1:ny+1), ushock(1:ny+1), frac(1:ny+1))
#elif defined (key_riemann_bigloop)
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))          &
 !$acc private (u(jmin:jmax,1:nvar), qxm(jmin:jmax,1:nvar), qxp(jmin:jmax,1:nvar) , &
 !$acc          q(jmin:jmax,1:nvar)  , dq(jmin:jmax,1:nvar) , c(jmin:jmax),         &
 !$acc          qleft(1:ny+1,1:nvar) , qright(1:ny+1,1:nvar),                       &
 !$acc          qgdnv(1:ny+1,1:nvar), flux(1:ny+1,1:nvar),                          &
 !$acc          rr(1:ny+1), ur(1:ny+1), pr(1:ny+1), cr(1:ny+1), wr(1:ny+1),         &
 !$acc          ro(1:ny+1), uo(1:ny+1), po(1:ny+1), co(1:ny+1), wo(1:ny+1),         &
 !$acc          rl(1:ny+1), ul(1:ny+1), pl(1:ny+1), cl(1:ny+1), wl(1:ny+1),         &
 !$acc          rstar(1:ny+1), ustar(1:ny+1), pstar(1:ny+1), cstar(1:ny+1),         & 
 !$acc          scr(1:ny+1), delp(1:ny+1), pold(1:ny+1), ind(1:ny+1),               &
 !$acc          sgnm(1:ny+1), spin(1:ny+1), spout(1:ny+1), ushock(1:ny+1), frac(1:ny+1))
#else
 !$acc parallel num_gangs(gang_number) present(uold(1:nx+4,1:ny+4,1:nvar))          &
 !$acc private (u(jmin:jmax,1:nvar), qxm(jmin:jmax,1:nvar), qxp(jmin:jmax,1:nvar) , &
 !$acc          q(jmin:jmax,1:nvar)  , dq(jmin:jmax,1:nvar), c(jmin:jmax),          &
 !$acc          qleft(1:ny+1,1:nvar) , qright(1:ny+1,1:nvar),                       &
 !$acc          qgdnv(1:ny+1,1:nvar), flux(1:ny+1,1:nvar)                           )
#endif

     !$acc loop gang 
     do i=imin+2,imax-2
        ! Gather conservative variables
        !$acc loop vector
        do j=jmin,jmax
           u(j,ID)=uold(i,j,ID)
           u(j,IU)=uold(i,j,IV)
           u(j,IV)=uold(i,j,IU)
           u(j,IP)=uold(i,j,IP)
        end do

        if(nvar>4)then
           !$acc loop seq
           do in = 5,nvar
              !$acc loop vector
              do j=jmin,jmax
                 u(j,in)=uold(i,j,in)
              end do
           end do
        end if

        ! Convert to primitive variables
        call constoprim(u(:,:),q(:,:),c(:))

        ! Characteristic tracing
        call trace(q(:,:),dq(:,:),c(:),qxm(:,:),qxp(:,:),dtdx)

        !$acc loop seq
        do in = 1, nvar
          !$acc loop vector
           do j = 1, ny+1
              qleft (j,in)=qxm(j+1,in)
              qright(j,in)=qxp(j+2,in)
           end do
        end do

        ! Solve Riemann problem at interfaces
#if defined (key_riemann_oldalg)
        call riemann(qleft(:,:),qright(:,:),qgdnv(:,:),   &
                     rl(:),ul(:),pl(:),cl(:),wl(:),       &
                     rr(:),ur(:),pr(:),cr(:),wr(:),       &
                     ro(:),uo(:),po(:),co(:),wo(:),       &
                     rstar(:),ustar(:),pstar(:),cstar(:), &
                     sgnm(:),spin(:),spout(:),ushock(:),  &
                     frac(:),scr(:),delp(:),pold(:),ind(:),ind2(:),n_new(:))
#elif defined (key_riemann_bigloop)
        call riemann_newalg_bigloop(qleft(:,:),qright(:,:),qgdnv(:,:), &
                                    rl(:),ul(:),pl(:),cl(:),wl(:),     &
                                    rr(:),ur(:),pr(:),cr(:),wr(:),     &
                                    ro(:),uo(:),po(:),co(:),wo(:),     &
                                    rstar(:),ustar(:),pstar(:),        &
                                    cstar(:),sgnm(:),spin(:),spout(:), &
                                    ushock(:),frac(:),scr(:),delp(:),  &
                                    pold(:),ind(:)                     )
#else
        call riemann_opt_mem(qleft(:,:),qright(:,:),qgdnv(:,:))
        ! renamned form of riemann_dev_PF2 for explicit purposes
#endif

        ! Compute fluxes
        call cmpflx(qgdnv(:,:),flux(:,:))

        ! Update conservative variables 
        !$acc loop vector
        do j=jmin+2,jmax-2
           uold(i,j,ID)=u(j,ID)+(flux(j-2,ID)-flux(j-1,ID))*dtdx
           uold(i,j,IU)=u(j,IV)+(flux(j-2,IV)-flux(j-1,IV))*dtdx
           uold(i,j,IV)=u(j,IU)+(flux(j-2,IU)-flux(j-1,IU))*dtdx
           uold(i,j,IP)=u(j,IP)+(flux(j-2,IP)-flux(j-1,IP))*dtdx
        end do
 
        if(nvar>4)then
           !$acc loop seq
           do in = 5,nvar
              !$acc loop vector
              do j=jmin+2,jmax-2
                 uold(i,j,in)=u(j,in)+(flux(j-2,in)-flux(j-1,in))*dtdx
              end do
           end do
        end if

      end do
      !$acc end parallel

     ! Deallocate work space
     call deallocate_work_space
  end if
  !$acc end data

contains

  subroutine allocate_work_space(ii1,ii2,ngrid)
    implicit none

    ! Dummy arguments
    integer(kind=prec_int), intent(in) :: ii1,ii2,ngrid
    allocate(u  (ii1:ii2,1:nvar))
    allocate(q  (ii1:ii2,1:nvar))
    allocate(dq (ii1:ii2,1:nvar))
    allocate(qxm(ii1:ii2,1:nvar))
    allocate(qxp(ii1:ii2,1:nvar))
    allocate(c  (ii1:ii2))
    allocate(qleft (1:ngrid,1:nvar))
    allocate(qright(1:ngrid,1:nvar))
    allocate(qgdnv (1:ngrid,1:nvar))
    allocate(flux  (1:ngrid,1:nvar))
#if defined (key_riemann_oldalg) 
    allocate(rl    (1:ngrid), ul   (1:ngrid), pl   (1:ngrid), cl    (1:ngrid))
    allocate(rr    (1:ngrid), ur   (1:ngrid), pr   (1:ngrid), cr    (1:ngrid))
    allocate(ro    (1:ngrid), uo   (1:ngrid), po   (1:ngrid), co    (1:ngrid))
    allocate(rstar (1:ngrid), ustar(1:ngrid), pstar(1:ngrid), cstar (1:ngrid))
    allocate(wl    (1:ngrid), wr   (1:ngrid), wo   (1:ngrid))
    allocate(sgnm  (1:ngrid), spin (1:ngrid), spout(1:ngrid), ushock(1:ngrid))
    allocate(frac  (1:ngrid), scr  (1:ngrid), delp (1:ngrid), pold  (1:ngrid))
    allocate(ind   (1:ngrid), ind2 (1:ngrid))
#elif defined (key_riemann_bigloop)
    allocate(rl    (1:ngrid), ul   (1:ngrid), pl   (1:ngrid), cl    (1:ngrid))
    allocate(rr    (1:ngrid), ur   (1:ngrid), pr   (1:ngrid), cr    (1:ngrid))
    allocate(ro    (1:ngrid), uo   (1:ngrid), po   (1:ngrid), co    (1:ngrid))
    allocate(rstar (1:ngrid), ustar(1:ngrid), pstar(1:ngrid), cstar (1:ngrid))
    allocate(wl    (1:ngrid), wr   (1:ngrid), wo   (1:ngrid))
    allocate(sgnm  (1:ngrid), spin (1:ngrid), spout(1:ngrid), ushock(1:ngrid))
    allocate(frac  (1:ngrid), scr  (1:ngrid), delp (1:ngrid), pold  (1:ngrid))
    allocate(ind   (1:ngrid) )
#endif
  end subroutine allocate_work_space

  subroutine deallocate_work_space
    deallocate(u,q,dq,qxm,qxp,c,qleft,qright,qgdnv,flux)
#if defined (key_riemann_oldalg)
    deallocate(rl,ul,pl,cl)
    deallocate(rr,ur,pr,cr)
    deallocate(ro,uo,po,co)
    deallocate(rstar,ustar,pstar,cstar)
    deallocate(wl,wr,wo)
    deallocate(sgnm,spin,spout,ushock)
    deallocate(frac,scr,delp,pold)
    deallocate(ind,ind2)
#elif defined (key_riemann_bigloop)
    deallocate(rl,ul,pl,cl)
    deallocate(rr,ur,pr,cr)
    deallocate(ro,uo,po,co)
    deallocate(rstar,ustar,pstar,cstar)
    deallocate(wl,wr,wo)
    deallocate(sgnm,spin,spout,ushock)
    deallocate(frac,scr,delp,pold,ind)
#endif 
  end subroutine deallocate_work_space

end subroutine godunov

end module hydro_principal
