!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module_hydro_utils.f90 --- 
!!!!
!! subroutine make_boundary(idim)
!! subroutine constoprim(u,q,c)
!! subroutine eos(rho,eint,p,c)
!! subroutine trace(q,dq,c,qxm,qxp,dtdx)
!! subroutine slope(q,dq)
!! subroutine riemann(qleft,qright,qgdnv, &
!!            & rl,ul,pl,cl,wl,rr,ur,pr,cr,wr,ro,uo,po,co,wo, &
!!            & rstar,ustar,pstar,cstar,sgnm,spin,spout, &
!!            & ushock,frac,scr,delp,pold,ind,ind2)
!! subroutine cmpflx(qgdnv,flux)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hydro_utils

  double precision,save,private :: tb

  interface
  function sync () bind(c,name="sync")
    use iso_c_binding, only: c_int
    integer(c_int) :: sync
  end function sync
  end interface

contains

subroutine make_boundary(idim)
  use hydro_mpi
  use hydro_commons
  use hydro_const
  use hydro_parameters
  implicit none

  ! Dummy arguments
  integer(kind=prec_int), intent(in) :: idim
  ! Local variables
  integer(kind=prec_int) :: ivar,i,i0,j,j0
  real(kind=prec_real)   :: sign
!!$ integer(kind=prec_int) :: ijet
!!$ real(kind=prec_real) :: djet,ujet,pjet

  integer(kind=prec_int) :: k,ik

  !$acc data present(uold(:,:,:), &
  !$acc              tabx1(1:2,3:ny+2,1:nvar), tabx2(1:2,3:ny+2,1:nvar), &
  !$acc              tabx3(1:2,3:ny+2,1:nvar), tabx4(1:2,3:ny+2,1:nvar), &
  !$acc              taby1(3:nx+2,1:2,1:nvar), taby2(3:nx+2,1:2,1:nvar), &
  !$acc              taby3(3:nx+2,1:2,1:nvar), taby4(3:nx+2,1:2,1:nvar)  )
  if(idim==1)then


     !--MPI--!
  !$acc parallel async(1)
  !$acc loop seq
  do k=1,nvar
   !$acc loop gang vector
   do j=3,ny+2
    !$acc loop seq
    do i=1,2 ! 3,nx
     tabx1(i,j,k) = uold(nx+i,j,k)
     tabx2(i,j,k) = uold(i,j,k)
     tabx3(i,j,k) = uold(i+2,j,k)
     tabx4(i,j,k) = uold(nx+i+2,j,k)
    enddo
   enddo
  enddo
  !$acc end parallel
 !$acc update host(tabx1(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update host(tabx2(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update host(tabx3(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update host(tabx4(1:2,3:ny+2,1:nvar)) async(1)
 !$acc wait(1)
   CALL MPI_SENDRECV(tabx1(1,3,1),2*nvar*ny,MPI_DOUBLE_PRECISION,voisins(EAST),etiquette, &
                     tabx2(1,3,1),2*nvar*ny,MPI_DOUBLE_PRECISION,voisins(WEST),etiquette, &
                     comm2d,MPI_STATUS_IGNORE,code)

  CALL MPI_SENDRECV( tabx3(1,3,1),2*nvar*ny,MPI_DOUBLE_PRECISION,voisins(WEST),etiquette, &
                     tabx4(1,3,1),2*nvar*ny,MPI_DOUBLE_PRECISION,voisins(EAST),etiquette, &
                     comm2d,MPI_STATUS_IGNORE,code)
 !$acc update device(tabx1(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update device(tabx2(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update device(tabx3(1:2,3:ny+2,1:nvar)) async(1)
 !$acc update device(tabx4(1:2,3:ny+2,1:nvar)) async(1)

  !$acc parallel async(1)
  !$acc loop seq
  do k=1,nvar
   !$acc loop gang vector
   do j=3,ny+2
    !$acc loop seq
    do i=1,2
     uold(nx+i,j,k)   = tabx1(i,j,k)
     uold(i,j,k)      = tabx2(i,j,k)
     uold(i+2,j,k)    = tabx3(i,j,k)
     uold(nx+i+2,j,k) = tabx4(i,j,k)
    enddo
   enddo
  enddo
  !$acc end parallel

  !$acc wait(1)

!     !--MPI--!

     !$acc parallel present(uold(:,:,:)) 
     ! Left boundary
     if (boundary_left>0) then
        !$acc loop seq
        do ivar=1,nvar
           !$acc loop seq
           do i=1,2           
              sign=1.0_prec_real
              if(boundary_left==1)then
                 i0=5-i
                 if(ivar==IU)sign=-1.0_prec_real
              else if(boundary_left==2)then
                 i0=3
              else
                 i0=nx+i
              end if
              !$acc loop
              do j=jmin+2,jmax-2
                 uold(i,j,ivar)=uold(i0,j,ivar)*sign
              end do
           end do
        end do
     end if

     ! Right boundary
     if (boundary_right>0) then
        !$acc loop seq
        do ivar=1,nvar
           !$acc loop seq
           do i=nx+3,nx+4
              sign=1.0_prec_real
              if(boundary_right==1)then
                 i0=2*nx+5-i
                 if(ivar==IU)sign=-1.0_prec_real
              else if(boundary_right==2)then
                 i0=nx+2
              else
                 i0=i-nx
              end if
              !$acc loop
              do j=jmin+2,jmax-2
                 uold(i,j,ivar)=uold(i0,j,ivar)*sign
              end do
           end do
        end do
     end if
    !$acc end parallel
  else

     !--MPI--!
  !$acc parallel async(3)
  !$acc loop seq
  do k=1,nvar
   !$acc loop seq
   do j=1,2
    !$acc loop gang vector
    do i=3,nx+2 ! 3,nx
     taby1(i,j,k)=uold(i,ny+j+2,k)
     taby2(i,j,k)=uold(i,j+2,k)
     taby3(i,j,k)=uold(i,ny+j,k)
     taby4(i,j,k)=uold(i,j,k)
    enddo
   enddo
  enddo
  !$acc end parallel
 !$acc update host(taby1(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update host(taby2(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update host(taby3(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update host(taby4(3:nx+2,1:2,1:nvar)) async(3)
 !$acc wait(3)
     ! Send to south and receive from north
  CALL MPI_SENDRECV( taby2(3,1,1),2*nvar*nx,MPI_DOUBLE_PRECISION,voisins(SOUTH),etiquette, &
                     taby1(3,1,1),2*nvar*nx,MPI_DOUBLE_PRECISION,voisins(NORTH),etiquette, &
                     comm2d,MPI_STATUS_IGNORE,code)
    ! Send to north and receive from south
  CALL MPI_SENDRECV( taby3(3,1,1),2*nvar*nx,MPI_DOUBLE_PRECISION,voisins(NORTH),etiquette, &
                     taby4(3,1,1),2*nvar*nx,MPI_DOUBLE_PRECISION,voisins(SOUTH),etiquette, &
                     comm2d,MPI_STATUS_IGNORE,code)
 !$acc update device(taby1(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update device(taby2(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update device(taby3(3:nx+2,1:2,1:nvar)) async(3)
 !$acc update device(taby4(3:nx+2,1:2,1:nvar)) async(3)

  !$acc parallel async(3)
  !$acc loop seq
  do k=1,nvar
   !$acc loop seq
   do j=1,2
    !$acc loop gang vector
    do i=3,nx+2
     uold(i,ny+j+2,k) = taby1(i,j,k)
     uold(i,j+2,k)    = taby2(i,j,k)
     uold(i,ny+j,k)   = taby3(i,j,k)
     uold(i,j,k)      = taby4(i,j,k)
    enddo
   enddo
  enddo
  !$acc end parallel
  !$acc wait(3)
     !--MPI--!

     !$acc parallel present(uold(:,:,:)) !async
     ! Lower boundary
     if (boundary_down>0) then
        !$acc loop seq
        do ivar=1,nvar
           !$acc loop seq
           do j=1,2           
              sign=1.0_prec_real
              if(boundary_down==1)then
                 j0=5-j
                 if(ivar==IV)sign=-1.0_prec_real
              else if(boundary_down==2)then
                 j0=3
              else
                 j0=ny+j
              end if
              !$acc loop
              do i=imin+2,imax-2
                 uold(i,j,ivar)=uold(i,j0,ivar)*sign
              end do
           end do
        end do
     endif

!!$        djet=1.0
!!$        ujet=300.
!!$        pjet=1.0
!!$        ijet=imax-20 
!!$        do ivar=1,nvar
!!$           do j=1,2           
!!$              do i=ijet,imax-2
!!$                 uold(i,j,1)=djet
!!$                 uold(i,j,2)=0.
!!$                 uold(i,j,3)=djet*ujet
!!$                 uold(i,j,4)=pjet/(gamma-1)+0.5*djet*ujet**2
!!$              end do
!!$           end do
!!$        end do

     ! Upper boundary
     if (boundary_up>0) then
        !$acc loop seq
        do ivar=1,nvar
           !$acc loop seq
           do j=ny+3,ny+4
              sign=1.0_prec_real
              if(boundary_up==1)then
                 j0=2*ny+5-j
                 if(ivar==IV)sign=-1.0_prec_real
              else if(boundary_up==2)then
                 j0=ny+2
              else
                 j0=j-ny
              end if
              !$acc loop
              do i=imin+2,imax-2
                 uold(i,j,ivar)=uold(i,j0,ivar)*sign
              end do
           end do
        end do
     endif
     !$acc end parallel
     !!$acc wait
  end if
  !$acc end data
end subroutine make_boundary


subroutine constoprim(u,q,c)
  use hydro_commons
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real), dimension(:,:), intent(in)   :: u
  real(kind=prec_real), dimension(:,:), intent(out)  :: q
  real(kind=prec_real), dimension(:)  , intent(out)  :: c
  ! Local variables
  integer(kind=prec_int) :: i,IN,ijmax,ijmin
  real(kind=prec_real)   :: eken
  real(kind=prec_real), dimension(1:size(c)) :: e

  ! ::: .....::::: convert to non-conservation form and compute
  ! ::: .....::::: pressure via equation of state
  !$acc data present(u(:,:),q(:,:),c(:))
  ijmin = 1
  ijmax=size(c)

  !$acc enter data create(e(1:ijmax))

  !$acc loop vector
  do i = ijmin, ijmax
     q(i,ID) = max(u(i,ID),smallr)
     q(i,IU) = u(i,IU)/q(i,ID)
     q(i,IV) = u(i,IV)/q(i,ID)
     eken = half*(q(i,IU)**2+q(i,IV)**2)
     q(i,IP) = u(i,IP)/q(i,ID) - eken     
  end do

  if (nvar > 4) then
     !$acc loop seq
     do IN = 5, nvar
        !$acc loop vector
        do i = ijmin, ijmax
           q(i,IN) = u(i,IN)/q(i,ID)
        end do
     end do
  end if
  !$acc loop vector
  do i = ijmin, ijmax
     e(i)=q(i,IP)
  end do

  call eos(q(1:ijmax,ID),e(1:ijmax),q(1:ijmax,IP),c(1:ijmax))
  !$acc exit data delete(e(1:ijmax))
  !$acc end data
end subroutine constoprim


subroutine eos(rho,eint,p,c)
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real), dimension(:), intent(in)  :: rho,eint
  real(kind=prec_real), dimension(:), intent(out) :: p,c
  ! Local variables
  integer(kind=prec_int) :: k
  real(kind=prec_real)   :: smallp
  !$acc data present(rho(:),eint(:),p(:),c(:))
  smallp = smallc**2/gamma
  !$acc loop  vector
  do k = 1, size(rho)
     p(k) = (gamma - one)*rho(k)*eint(k)
     p(k) = max(p(k),rho(k)*smallp)
     c(k) = sqrt(gamma*p(k)/rho(k))
  end do
  !$acc end data
end subroutine eos


subroutine trace(q,dq,c,qxm,qxp,dtdx)
  use hydro_commons
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real), intent(in) :: dtdx
  real(kind=prec_real), dimension(:,:), intent(in)  :: q
  real(kind=prec_real), dimension(:,:), intent(out) :: dq,qxm,qxp 
  real(kind=prec_real), dimension(:),   intent(in)  :: c 
  ! Local variables
  integer(kind=prec_int) :: i,IN,ijmin,ijmax,j
  real(kind=prec_real)   :: cc,csq,r,u,v,p,a
  real(kind=prec_real)   :: dr,du,dv,dp,da
  real(kind=prec_real)   :: alpham,alphap,alpha0r,alpha0v
  real(kind=prec_real)   :: spminus,spplus,spzero
  real(kind=prec_real)   :: apright,amright,azrright,azv1right,acmpright
  real(kind=prec_real)   :: apleft,amleft,azrleft,azv1left,acmpleft
  real(kind=prec_real)   :: zerol,zeror,project
  !$acc data present(q(:,:),dq(:,:),qxm(:,:),qxp(:,:),c(:))
  ijmin = 1
  ijmax = size(c)


  ! compute slopes
  !$acc loop seq
  do j=1,nvar
   !$acc loop vector
   do i=ijmin,ijmax
    dq(i,j) = zero
   enddo
  enddo
  if (iorder .ne. 1) then
     call slope(q(:,:),dq(:,:))
  endif

  if(scheme=='muscl')then ! MUSCL-Hancock method
     zerol=-100.0d0/dtdx
     zeror=+100.0d0/dtdx
     project=one
  endif
  if(scheme=='plmde')then ! standard PLMDE
     zerol=zero
     zeror=zero
     project=one  
  endif
  if(scheme=='collela')then ! Collela's method
     zerol=zero
     zeror=zero
     project=zero
  endif

  !$acc loop vector
  do i = ijmin+1, ijmax-1
     
     cc = c(i)
     csq = cc**2
     r = q(i,ID)
     u = q(i,IU)
     v = q(i,IV)
     p = q(i,IP)
     
     dr = dq(i,ID)
     du = dq(i,IU)
     dv = dq(i,IV)
     dp = dq(i,IP)
     
     alpham  = half*(dp/(r*cc) - du)*r/cc
     alphap  = half*(dp/(r*cc) + du)*r/cc
     alpha0r = dr - dp/csq
     alpha0v = dv
     
     ! Right state
     spminus = (u-cc)*dtdx+one
     spplus  = (u+cc)*dtdx+one
     spzero  =  u    *dtdx+one
     if((u-cc)>=zeror)spminus=project
     if((u+cc)>=zeror)spplus =project
     if( u    >=zeror)spzero =project
     apright   = -half*spplus *alphap
     amright   = -half*spminus*alpham
     azrright  = -half*spzero *alpha0r
     azv1right = -half*spzero *alpha0v     
     qxp(i,ID) = r + (apright + amright + azrright)
     qxp(i,IU) = u + (apright - amright           )*cc/r
     qxp(i,IV) = v + (azv1right                   )
     qxp(i,IP) = p + (apright + amright           )*csq
     
     ! Left state
     spminus = (u-cc)*dtdx-one
     spplus  = (u+cc)*dtdx-one
     spzero  =  u    *dtdx-one
     if((u-cc)<=zerol)spminus=-project
     if((u+cc)<=zerol)spplus =-project
     if( u    <=zerol)spzero =-project
     apleft   = -half*spplus *alphap
     amleft   = -half*spminus*alpham
     azrleft  = -half*spzero *alpha0r
     azv1left = -half*spzero *alpha0v
     qxm(i,ID) = r + (apleft + amleft + azrleft)
     qxm(i,IU) = u + (apleft - amleft          )*cc/r
     qxm(i,IV) = v + (azv1left                 )
     qxm(i,IP) = p + (apleft + amleft          )*csq
     
  end do

  if (nvar>4)then
     !$acc loop seq
     do IN = 5, nvar
        !$acc loop vector
        do i = ijmin+1, ijmax-1
           
           u  =  q(i,IU)
           a  =  q(i,IN)
           da = dq(i,IN)
           ! Right state
           spzero = u*dtdx+one
           if(u>=zeror)spzero=project
           acmpright = -half*spzero*da
           qxp(i,IN) = a + acmpright
           ! Left state
           spzero = u*dtdx-one
           if(u<=zerol)spzero=-project
           acmpleft = -half*spzero*da
           qxm(i,IN) = a + acmpleft
           
        end do
     end do
  end if
  !$acc end data
end subroutine trace


subroutine slope(q,dq)
  use hydro_commons
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real), dimension(:,:), intent(in)  :: q
  real(kind=prec_real), dimension(1:size(q,1),1:nvar), intent(out) :: dq
  ! Local arrays
  integer(kind=prec_int) :: i,n,ijmin,ijmax
  real(kind=prec_real)   :: dsgn,dlim,dcen,dlft,drgt,slop
  !$acc data present(q(:,:),dq(:,:))
  ijmin = 1
  ijmax = size(q,1)
  !$acc loop seq
  do n = 1, nvar
     !$acc loop vector
     do i = ijmin+1, ijmax-1        
        dlft = slope_type*(q(i  ,n) - q(i-1,n))
        drgt = slope_type*(q(i+1,n) - q(i  ,n))
        dcen = half*(dlft+drgt)/slope_type
        dsgn = sign(one, dcen)
        slop = min(abs(dlft),abs(drgt))
        dlim = slop
        if((dlft*drgt)<=zero)dlim=zero
        dq(i,n) = dsgn*min(dlim,abs(dcen))
     end do
  end do  
  !$acc end data
end subroutine slope


subroutine riemann(qleft,qright,qgdnv, &
                   rl,ul,pl,cl,wl,rr,ur,pr,cr,wr,ro,uo,po,co,wo, &
                   rstar,ustar,pstar,cstar,sgnm,spin,spout, &
                   ushock,frac,scr,delp,pold,ind,ind2,n_new)
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real),dimension(:,:), intent(in)  :: qleft,qright
  real(kind=prec_real),dimension(:,:), intent(out) :: qgdnv
  real(kind=prec_real),dimension(:),   intent(out) :: rl,ul,pl,cl,wl
  real(kind=prec_real),dimension(:),   intent(out) :: rr,ur,pr,cr,wr
  real(kind=prec_real),dimension(:),   intent(out) :: ro,uo,po,co,wo
  real(kind=prec_real),dimension(:),   intent(out) :: rstar,ustar,pstar,cstar
  real(kind=prec_real),dimension(:),   intent(out) :: sgnm,spin,spout,ushock
  real(kind=prec_real),dimension(:),   intent(out) :: frac,scr,delp,pold  
  integer(kind=prec_int),dimension(:), intent(out) :: ind,ind2,n_new

  ! Local variables
  real(kind=prec_real)   :: smallp,gamma6,ql,qr,usr,usl,wwl,wwr,smallpp
  !integer(kind=prec_int) :: i,j,in,n,iter,n_new,nface
  integer(kind=prec_int) :: i,j,in,n,iter,nface
  integer(kind=prec_int) :: acc_id

  !$acc data present(qleft(:,:),qright(:,:),qgdnv(:,:),     & 
  !$acc              rstar(:),ustar(:),pstar(:),cstar(:),   &
  !$acc              sgnm(:),spin(:),spout(:),ushock(:),    &
  !$acc              rl(:),ul(:),pl(:),cl(:),wl(:),         & 
  !$acc              rr(:),ur(:),pr(:),cr(:),wr(:),         &
  !$acc              ro(:),uo(:),po(:),co(:),wo(:),         &
  !$acc              frac(:),scr(:),delp(:),pold(:),ind(:), &
  !$acc              ind2(:),n_new(1)) 

  ! Constants
  nface=size(rl)
  smallp = smallc**2/gamma
  smallpp = smallr*smallp
  gamma6 = (gamma+one)/(two*gamma)

  ! Pressure, density and velocity
  !$acc loop
  do i=1,nface
     rl(i)=MAX(qleft (i,ID),smallr)
     ul(i)=    qleft (i,IU)
     pl(i)=MAX(qleft (i,IP),rl(i)*smallp)
     rr(i)=MAX(qright(i,ID),smallr)
     ur(i)=    qright(i,IU)
     pr(i)=MAX(qright(i,IP),rr(i)*smallp)
  end do

  ! Lagrangian sound speed
  !$acc loop
  do i=1,nface
     cl(i) = gamma*pl(i)*rl(i)
     cr(i) = gamma*pr(i)*rr(i)
  end do

  ! First guess
  !$acc loop
  do i=1,nface
   wl(i)    = sqrt(cl(i))
   wr(i)    = sqrt(cr(i))
   pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i)))/(wl(i)+wr(i))
   pstar(i) = MAX(pstar(i),0.d0)
   pold(i)  = pstar(i)
  enddo
  n = nface
  !$acc loop
  do i = 1,n
     ind(i)=i
  end do

  ! Newton-Raphson iterations to find pstar at the required accuracy
  !$acc loop seq
  do iter = 1,niter_riemann 
     !$acc loop
     do i=1,n
        wwl=sqrt(cl(ind(i))*(one+gamma6*(pold(i)-pl(ind(i)))/pl(ind(i))))
        wwr=sqrt(cr(ind(i))*(one+gamma6*(pold(i)-pr(ind(i)))/pr(ind(i))))
        ql=two*wwl**3/(wwl**2+cl(ind(i)))
        qr=two*wwr**3/(wwr**2+cr(ind(i)))
        usl=ul(ind(i))-(pold(i)-pl(ind(i)))/wwl
        usr=ur(ind(i))+(pold(i)-pr(ind(i)))/wwr
        delp(i)=MAX(qr*ql/(qr+ql)*(usl-usr),-pold(i))
     end do
     !$acc loop
     do i=1,n
        pold(i)=pold(i)+delp(i)
     end do
     ! Convergence indicator
     !$acc loop
     do i=1,n 
        uo(i)=ABS(delp(i)/(pold(i)+smallpp))
     end do
     n_new(1)=0
     !$acc loop 
     do i=1,n
        if(uo(i)>1.d-06)then
           !$acc atomic capture
           n_new(1)=n_new(1)+1
           acc_id  =n_new(1)
           !$acc end atomic
           ind2(acc_id)=ind (i)
           po  (acc_id)=pold(i)
        end if
     end do
     
     j=n_new(1)
     !$acc loop
     do i=1,n
        if(uo(i)<=1.d-06)then
           !$acc atomic capture
           n_new(1)=n_new(1)+1
           acc_id  =n_new(1)
           !$acc end atomic
           ind2(acc_id)=ind (i)
           po  (acc_id)=pold(i)
        end if
     end do

     !$acc loop
     do i=1,n
      ind (i)=ind2(i)
      pold(i)=po  (i)
     enddo
     n=j
  end do

  !$acc loop
  do i=1,nface
     pstar(ind(i))=pold(i)
  end do
  !$acc loop
  do i=1,nface
     wl(i) = sqrt(cl(i)*(one+gamma6*(pstar(i)-pl(i))/pl(i)))
     wr(i) = sqrt(cr(i)*(one+gamma6*(pstar(i)-pr(i))/pr(i)))
  end do
  !$acc loop
  do i=1,nface
     ustar(i) = half*(ul(i) + (pl(i)-pstar(i))/wl(i) + &
          &           ur(i) - (pr(i)-pstar(i))/wr(i) )
  end do
  !$acc loop
  do i=1,nface
     sgnm(i) = sign(one,ustar(i))
  end do
  !$acc loop
  do i=1,nface
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
     end if
  end do
  !$acc loop
  do i=1,nface
     co(i) = max(smallc,sqrt(abs(gamma*po(i)/ro(i))))
  end do
  !$acc loop
  do i=1,nface
     rstar(i) = ro(i)/(one+ro(i)*(po(i)-pstar(i))/wo(i)**2)
     rstar(i) = max(rstar(i),smallr)
  end do
  !$acc loop
  do i=1,nface
     cstar(i) = max(smallc,sqrt(abs(gamma*pstar(i)/rstar(i))))
  end do
  !$acc loop
  do i=1,nface
     spout (i) = co   (i)    - sgnm(i)*uo   (i)
     spin  (i) = cstar(i)    - sgnm(i)*ustar(i)
     ushock(i) = wo(i)/ro(i) - sgnm(i)*uo   (i)
  end do
  !$acc loop 
  do i=1,nface
     if(pstar(i)>=po(i))then
        spin (i) = ushock(i)
        spout(i) = ushock(i)
     end if
  end do
  !$acc loop 
  do i=1,nface
     scr(i) = MAX(spout(i)-spin(i),smallc+ABS(spout(i)+spin(i)))
  end do
  !$acc loop 
  do i=1,nface
     frac(i) = (one + (spout(i) + spin(i))/scr(i))*half
     frac(i) = max(zero,min(one,frac(i)))
  end do
  !$acc loop 
  do i=1,nface
     qgdnv(i,ID) = frac(i)*rstar(i) + (one - frac(i))*ro(i)
     qgdnv(i,IU) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
     qgdnv(i,IP) = frac(i)*pstar(i) + (one - frac(i))*po(i)
  end do
  !$acc loop 
  do i=1,nface
     if(spout(i)<zero)then
        qgdnv(i,ID) = ro(i)
        qgdnv(i,IU) = uo(i)
        qgdnv(i,IP) = po(i)
     end if
  end do
  !$acc loop 
  do i=1,nface
     if(spin(i)>zero)then
        qgdnv(i,ID) = rstar(i)
        qgdnv(i,IU) = ustar(i)
        qgdnv(i,IP) = pstar(i)
     end if
  end do

  ! transverse velocity
  !$acc loop 
  do i=1,nface
     if(sgnm(i)==one)then
        qgdnv(i,IV) = qleft (i,IV)
     else
        qgdnv(i,IV) = qright(i,IV)
     end if
  end do

  ! other passive variables
  if (nvar > 4) then
     !$acc loop seq
     do in = 4,nvar
        !$acc loop vector
        do i=1,nface
           if(sgnm(i)==one)then
              qgdnv(i,in) = qleft (i,in)
           else
              qgdnv(i,in) = qright(i,in)
           end if
        end do
     end do
  end if
  !$acc end data
end subroutine riemann

subroutine riemann_opt_mem(qleft,qright,qgdnv)
  use hydro_const
  use hydro_parameters
  implicit none
  !$acc routine vector

  ! Dummy arguments
  real(kind=prec_real),dimension(:,:), intent(in)  :: qleft,qright
  real(kind=prec_real),dimension(:,:), intent(out) :: qgdnv
  ! Local variables
  real(kind=prec_real)   :: smallp,gamma6,ql,qr,usr,usl,wwl,wwr,smallpp
  real(kind=prec_real)   :: rl,ul,pl,rr,ur,pr,cr,cl,wl,wr,pstar,pold,ustar,sgnm
  real(kind=prec_real)   :: delp,ro,uo,po,co,wo,rstar,cstar,spin,spout,ushock,frac,scr
  integer(kind=prec_int) :: i,j,in,n,iter,nface,ind
  !$acc data present(qleft(:,:),qright(:,:),qgdnv(:,:))

  ! Constants
  nface=size(qgdnv(:,1))
  smallp = smallc**2/gamma
  smallpp = smallr*smallp
  gamma6 = (gamma+one)/(two*gamma)

  ! Pressure, density and velocity

  !$acc loop  
  do i=1,nface
     rl=MAX(qleft (i,ID),smallr)
     ul=    qleft (i,IU)
     pl=MAX(qleft (i,IP),rl*smallp)
     rr=MAX(qright(i,ID),smallr)
     ur=    qright(i,IU)
     pr=MAX(qright(i,IP),rr*smallp)
     ! Lagrangian sound speed
     cl = gamma*pl*rl
     cr = gamma*pr*rr
     wl = sqrt(cl); wr = sqrt(cr)
     pstar = ((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr)
     pstar = MAX(pstar,0.d0)
     pold = pstar
     ind = 1
     iter = 1
     ! Newton-Raphson iterations to find pstar at the required accuracy
     do while(ind==1.and.iter<=niter_riemann)
        wwl=sqrt(cl*(one+gamma6*(pold-pl)/pl))
        wwr=sqrt(cr*(one+gamma6*(pold-pr)/pr))
        ql=two*wwl**3/(wwl**2+cl)
        qr=two*wwr**3/(wwr**2+cr)
        usl=ul-(pold-pl)/wwl
        usr=ur+(pold-pr)/wwr
        delp=MAX(qr*ql/(qr+ql)*(usl-usr),-pold)
        pold=pold+delp
        ! Convergence indicator
        uo=ABS(delp/(pold+smallpp))
        if (uo<=1.d-06) ind=0
        iter=iter+1
     end do
     pstar=pold
     wl = sqrt(cl*(one+gamma6*(pstar-pl)/pl))
     wr = sqrt(cr*(one+gamma6*(pstar-pr)/pr))
     ustar = half*(ul + (pl-pstar)/wl + &
          &           ur - (pr-pstar)/wr )
     sgnm = sign(one,ustar)
     if(sgnm==one)then
        ro = rl
        uo = ul
        po = pl
        wo = wl
     else
        ro = rr
        uo = ur
        po = pr
        wo = wr
     end if
     co = max(smallc,sqrt(abs(gamma*po/ro)))
     rstar = ro/(one+ro*(po-pstar)/wo**2)
     rstar = max(rstar,smallr)
     cstar = max(smallc,sqrt(abs(gamma*pstar/rstar)))
     spout  = co       - sgnm*uo
     spin   = cstar    - sgnm*ustar
     ushock = wo/ro - sgnm*uo
     if(pstar>=po)then
        spin  = ushock
        spout = ushock
     end if
     scr = MAX(spout-spin,smallc+ABS(spout+spin))
     frac = (one + (spout + spin)/scr)*half
     frac = max(zero,min(one,frac))
     qgdnv(i,ID) = frac*rstar + (one - frac)*ro
     qgdnv(i,IU) = frac*ustar + (one - frac)*uo
     qgdnv(i,IP) = frac*pstar + (one - frac)*po
     if(spin>zero)then
        qgdnv(i,ID) = rstar
        qgdnv(i,IU) = ustar
        qgdnv(i,IP) = pstar
     end if
  ! transverse velocity
     if(sgnm==one)then
        qgdnv(i,IV) = qleft (i,IV)
     else
        qgdnv(i,IV) = qright(i,IV)
     end if
  end do


! BUG, double the computing time if not commented while nvar=4
! even with -Mvect=nocond
  ! other passive variables
  !if (nvar > 4) then
  !   !$acc loop seq
  !   do in = 4,nvar
  !      !$acc loop vector
  !      do i=1,nface
  !         if(sgnm==one)then
  !            qgdnv(i,in) = qleft (i,in)
  !         else
  !            qgdnv(i,in) = qright(i,in)
  !         end if
  !      end do
  !   end do
  !end if
  !$acc end data
end subroutine riemann_opt_mem



subroutine cmpflx(qgdnv,flux)
  use hydro_parameters
  use hydro_const
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real), dimension(:,:), intent(in)  :: qgdnv
  real(kind=prec_real), dimension(:,:), intent(out) :: flux
  ! Local variables
  integer(kind=prec_int) :: i,IN,nface
  real(kind=prec_real)   :: entho,etot,ekin
  !$acc data present(qgdnv(:,:),flux(:,:))
  nface=size(flux,1)
  entho=one/(gamma-one)

  ! Compute fluxes
  !$acc loop 
  do i=1,nface

     ! Mass density
     flux(i,ID) = qgdnv(i,ID)*qgdnv(i,IU)

     ! Normal momentum
     flux(i,IU) = flux(i,ID)*qgdnv(i,IU)+qgdnv(i,IP)

     ! Transverse momentum 1
     flux(i,IV) = flux(i,ID)*qgdnv(i,IV)

     ! Total energy
     ekin = half*qgdnv(i,ID)*(qgdnv(i,IU)**2+qgdnv(i,IV)**2)
     etot = qgdnv(i,IP)*entho + ekin
     flux(i,IP) = qgdnv(i,IU)*(etot+qgdnv(i,IP))

     ! Other advected quantities
     if (nvar>4) then
        !$acc loop seq
        do IN = 5, nvar
           flux(i,IN) = flux(i,ID)*qgdnv(i,IN)
        end do
     end if

  end do
  !$acc end data
end subroutine cmpflx


subroutine riemann_newalg_f90(qleft,qright,qgdnv, &
                   rl,ul,pl,cl,wl,rr,ur,pr,cr,wr,ro,uo,po,co,wo, &
                   rstar,ustar,pstar,cstar,sgnm,spin,spout, &
                   ushock,frac,scr,delp,pold,ind)
  use hydro_const
  use hydro_parameters
  implicit none

  ! Dummy arguments
  real(kind=prec_real),dimension(:,:), intent(in)  :: qleft,qright
  real(kind=prec_real),dimension(:,:), intent(out) :: qgdnv
  real(kind=prec_real),dimension(:),   intent(out) :: rl,ul,pl,cl,wl
  real(kind=prec_real),dimension(:),   intent(out) :: rr,ur,pr,cr,wr
  real(kind=prec_real),dimension(:),   intent(out) :: ro,uo,po,co,wo
  real(kind=prec_real),dimension(:),   intent(out) :: rstar,ustar,pstar,cstar
  real(kind=prec_real),dimension(:),   intent(out) :: sgnm,spin,spout,ushock
  real(kind=prec_real),dimension(:),   intent(out) :: frac,scr,delp,pold  
  integer(kind=prec_int),dimension(:), intent(out) :: ind
  ! Local variables
  real(kind=prec_real)   :: smallp,gamma6,ql,qr,usr,usl,wwl,wwr,smallpp
  integer(kind=prec_int) :: i,j,in,n,iter,nface

  ! Constants
  nface=size(rl)
  smallp = smallc**2/gamma
  smallpp = smallr*smallp
  gamma6 = (gamma+one)/(two*gamma)

  ! Pressure, density and velocity
  rl(1:nface)=MAX(qleft (1:nface,ID),smallr)
  ul(1:nface)=    qleft (1:nface,IU)
  pl(1:nface)=MAX(qleft (1:nface,IP),rl(1:nface)*smallp)
  rr(1:nface)=MAX(qright(1:nface,ID),smallr)
  ur(1:nface)=    qright(1:nface,IU)
  pr(1:nface)=MAX(qright(1:nface,IP),rr(1:nface)*smallp)
  ! Lagrangian sound speed
  cl(1:nface) = gamma*pl(1:nface)*rl(1:nface)
  cr(1:nface) = gamma*pr(1:nface)*rr(1:nface)
  ! First guess
  wl(1:nface) = sqrt(cl(1:nface))
  wr(1:nface) = sqrt(cr(1:nface))
  pstar(1:nface) = ((wr(1:nface)*pl(1:nface)+wl(1:nface)*pr(1:nface))+ &
                   wl(1:nface)*wr(1:nface)*(ul(1:nface)-ur(1:nface)))/(wl(1:nface)+wr(1:nface))
  pstar(1:nface) = MAX(pstar(1:nface),0.d0)
  pold(1:nface) = pstar(1:nface)
  ind(1:nface) = 1
  ! Newton-Raphson iterations to find pstar at the required accuracy
  do iter = 1,niter_riemann 
     do i=1,nface
        if (ind(i)==1) then
           wwl=sqrt(cl(i)*(one+gamma6*(pold(i)-pl(i))/pl(i)))
           wwr=sqrt(cr(i)*(one+gamma6*(pold(i)-pr(i))/pr(i)))
           ql=two*wwl**3/(wwl**2+cl(i))
           qr=two*wwr**3/(wwr**2+cr(i))
           usl=ul(i)-(pold(i)-pl(i))/wwl
           usr=ur(i)+(pold(i)-pr(i))/wwr
           delp(i)=MAX(qr*ql/(qr+ql)*(usl-usr),-pold(i))
           pold(i)=pold(i)+delp(i)
           ! Convergence indicator
           uo(i)=ABS(delp(i)/(pold(i)+smallpp))
           if (uo(i)<=1.d-06) ind(i)=0
        end if
     end do
  end do
  pstar(1:nface)=pold(1:nface)
  wl(1:nface) = sqrt(cl(1:nface)*(one+gamma6*(pstar(1:nface)-pl(1:nface))/pl(1:nface)))
  wr(1:nface) = sqrt(cr(1:nface)*(one+gamma6*(pstar(1:nface)-pr(1:nface))/pr(1:nface)))
  ustar(1:nface) = half*(ul(1:nface) + (pl(1:nface)-pstar(1:nface))/wl(1:nface) + &
       &           ur(1:nface) - (pr(1:nface)-pstar(1:nface))/wr(1:nface) )
  sgnm(1:nface) = sign(one,ustar(1:nface))
  where(sgnm(1:nface)==one)
     ro(1:nface) = rl(1:nface)
     uo(1:nface) = ul(1:nface)
     po(1:nface) = pl(1:nface)
     wo(1:nface) = wl(1:nface)
  elsewhere
     ro(1:nface) = rr(1:nface)
     uo(1:nface) = ur(1:nface)
     po(1:nface) = pr(1:nface)
     wo(1:nface) = wr(1:nface)
  end where
  co(1:nface) = max(smallc,sqrt(abs(gamma*po(1:nface)/ro(1:nface))))
  rstar(1:nface) = ro(1:nface)/(one+ro(1:nface)*(po(1:nface)-pstar(1:nface))/wo(1:nface)**2)
  rstar(1:nface) = max(rstar(1:nface),smallr)
  cstar(1:nface) = max(smallc,sqrt(abs(gamma*pstar(1:nface)/rstar(1:nface))))
  spout (1:nface) = co   (1:nface)    - sgnm(1:nface)*uo   (1:nface)
  spin  (1:nface) = cstar(1:nface)    - sgnm(1:nface)*ustar(1:nface)
  ushock(1:nface) = wo(1:nface)/ro(1:nface) - sgnm(1:nface)*uo   (1:nface)
  where(pstar(1:nface)>=po(1:nface))
     spin (1:nface) = ushock(1:nface)
     spout(1:nface) = ushock(1:nface)
  end where
  scr(1:nface) = MAX(spout(1:nface)-spin(1:nface),smallc+ABS(spout(1:nface)+spin(1:nface)))
  frac(1:nface) = (one + (spout(1:nface) + spin(1:nface))/scr(1:nface))*half
  frac(1:nface) = max(zero,min(one,frac(1:nface)))
  qgdnv(1:nface,ID) = frac(1:nface)*rstar(1:nface) + (one - frac(1:nface))*ro(1:nface)
  qgdnv(1:nface,IU) = frac(1:nface)*ustar(1:nface) + (one - frac(1:nface))*uo(1:nface)
  qgdnv(1:nface,IP) = frac(1:nface)*pstar(1:nface) + (one - frac(1:nface))*po(1:nface)
  where(spout(1:nface)<zero)
     qgdnv(1:nface,ID) = ro(1:nface)
     qgdnv(1:nface,IU) = uo(1:nface)
     qgdnv(1:nface,IP) = po(1:nface)
  end where
  where(spin(1:nface)>zero)
     qgdnv(1:nface,ID) = rstar(1:nface)
     qgdnv(1:nface,IU) = ustar(1:nface)
     qgdnv(1:nface,IP) = pstar(1:nface)
  end where
  ! transverse velocity
  where(sgnm(1:nface)==one)
     qgdnv(1:nface,IV) = qleft (1:nface,IV)
  elsewhere
     qgdnv(1:nface,IV) = qright(1:nface,IV)
  end where

  ! other passive variables
  if (nvar > 4) then
     do in = 4,nvar
        do i=1,nface
           if(sgnm(i)==one)then
              qgdnv(i,in) = qleft (i,in)
           else
              qgdnv(i,in) = qright(i,in)
           end if
        end do
     end do
  end if
  
end subroutine riemann_newalg_f90


subroutine riemann_newalg_bigloop(qleft,qright,qgdnv, &
                   rl,ul,pl,cl,wl,rr,ur,pr,cr,wr,ro,uo,po,co,wo, &
                   rstar,ustar,pstar,cstar,sgnm,spin,spout, &
                   ushock,frac,scr,delp,pold,ind)
  use hydro_const
  use hydro_parameters
  !$acc routine vector
  implicit none

  ! Dummy arguments
  real(kind=prec_real),dimension(:,:), intent(in)  :: qleft,qright
  real(kind=prec_real),dimension(:,:), intent(out) :: qgdnv
  real(kind=prec_real),dimension(:),   intent(out) :: rl,ul,pl,cl,wl
  real(kind=prec_real),dimension(:),   intent(out) :: rr,ur,pr,cr,wr
  real(kind=prec_real),dimension(:),   intent(out) :: ro,uo,po,co,wo
  real(kind=prec_real),dimension(:),   intent(out) :: rstar,ustar,pstar,cstar
  real(kind=prec_real),dimension(:),   intent(out) :: sgnm,spin,spout,ushock
  real(kind=prec_real),dimension(:),   intent(out) :: frac,scr,delp,pold  
  integer(kind=prec_int),dimension(:), intent(out) :: ind
  ! Local variables
  real(kind=prec_real)   :: smallp,gamma6,ql,qr,usr,usl,wwl,wwr,smallpp
  integer(kind=prec_int) :: i,j,in,n,iter,nface
  !$acc data present(qleft(:,:),qright(:,:),qgdnv(:,:),   & 
  !$acc              rstar(:),ustar(:),pstar(:),cstar(:), &
  !$acc              sgnm(:),spin(:),spout(:),ushock(:),  &
  !$acc              rl(:),ul(:),pl(:),cl(:),wl(:),       & 
  !$acc              rr(:),ur(:),pr(:),cr(:),wr(:),       &
  !$acc              ro(:),uo(:),po(:),co(:),wo(:),       &
  !$acc              frac(:),scr(:),delp(:),pold(:),ind(:))

  ! Constants
  nface=size(rl)
  smallp = smallc**2/gamma
  smallpp = smallr*smallp
  gamma6 = (gamma+one)/(two*gamma)

  ! Pressure, density and velocity
  !$acc loop
  do i=1,nface
     rl(i)=MAX(qleft (i,ID),smallr)
     ul(i)=    qleft (i,IU)
     pl(i)=MAX(qleft (i,IP),rl(i)*smallp)
     rr(i)=MAX(qright(i,ID),smallr)
     ur(i)=    qright(i,IU)
     pr(i)=MAX(qright(i,IP),rr(i)*smallp)
     ! Lagrangian sound speed
     cl(i) = gamma*pl(i)*rl(i)
     cr(i) = gamma*pr(i)*rr(i)
     ! First guess
     wl(i) = sqrt(cl(i))
     wr(i) = sqrt(cr(i))
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i)))/(wl(i)+wr(i))
     pstar(i) = MAX(pstar(i),0.d0)
     pold(i) = pstar(i)
     ind(i) = 1
     ! Newton-Raphson iterations to find pstar at the required accuracy
     do iter = 1,niter_riemann 
        if (ind(i)==1) then
           wwl=sqrt(cl(i)*(one+gamma6*(pold(i)-pl(i))/pl(i)))
           wwr=sqrt(cr(i)*(one+gamma6*(pold(i)-pr(i))/pr(i)))
           ql=two*wwl**3/(wwl**2+cl(i))
           qr=two*wwr**3/(wwr**2+cr(i))
           usl=ul(i)-(pold(i)-pl(i))/wwl
           usr=ur(i)+(pold(i)-pr(i))/wwr
           delp(i)=MAX(qr*ql/(qr+ql)*(usl-usr),-pold(i))
           pold(i)=pold(i)+delp(i)
           ! Convergence indicator
           uo(i)=ABS(delp(i)/(pold(i)+smallpp))
           if (uo(i)<=1.d-06) ind(i)=0
        end if
     end do
     pstar(i)=pold(i)
     wl(i) = sqrt(cl(i)*(one+gamma6*(pstar(i)-pl(i))/pl(i)))
     wr(i) = sqrt(cr(i)*(one+gamma6*(pstar(i)-pr(i))/pr(i)))
     ustar(i) = half*(ul(i) + (pl(i)-pstar(i))/wl(i) + ur(i) - (pr(i)-pstar(i))/wr(i) )
     sgnm(i) = sign(one,ustar(i))
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
     end if
     co(i) = max(smallc,sqrt(abs(gamma*po(i)/ro(i))))
     rstar(i) = ro(i)/(one+ro(i)*(po(i)-pstar(i))/wo(i)**2)
     rstar(i) = max(rstar(i),smallr)
     cstar(i) = max(smallc,sqrt(abs(gamma*pstar(i)/rstar(i))))
     spout (i) = co   (i)    - sgnm(i)*uo   (i)
     spin  (i) = cstar(i)    - sgnm(i)*ustar(i)
     ushock(i) = wo(i)/ro(i) - sgnm(i)*uo   (i)
     if(pstar(i)>=po(i))then
        spin (i) = ushock(i)
        spout(i) = ushock(i)
     end if
     scr(i) = MAX(spout(i)-spin(i),smallc+ABS(spout(i)+spin(i)))
     frac(i) = (one + (spout(i) + spin(i))/scr(i))*half
     frac(i) = max(zero,min(one,frac(i)))
     qgdnv(i,ID) = frac(i)*rstar(i) + (one - frac(i))*ro(i)
     qgdnv(i,IU) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
     qgdnv(i,IP) = frac(i)*pstar(i) + (one - frac(i))*po(i)
     if (spout(i)<zero) then
        qgdnv(i,ID) = ro(i)
        qgdnv(i,IU) = uo(i)
        qgdnv(i,IP) = po(i)
     end if
     if (spin(i)>zero) then
        qgdnv(i,ID) = rstar(i)
        qgdnv(i,IU) = ustar(i)
        qgdnv(i,IP) = pstar(i)
     end if
     ! transverse velocity
     if (sgnm(i)==one) then
        qgdnv(i,IV) = qleft (i,IV)
     else
        qgdnv(i,IV) = qright(i,IV)
     end if
  end do

  ! other passive variables
  !if (nvar > 4) then
  !   do in = 4,nvar
  !      do i=1,nface
  !         if(sgnm(i)==one)then
  !            qgdnv(i,in) = qleft (i,in)
  !         else
  !            qgdnv(i,in) = qright(i,in)
  !         end if
  !      end do
  !   end do
  !end if
  !$acc end data
end subroutine riemann_newalg_bigloop


subroutine title(n,nchar)
  use hydro_precision
  implicit none

  integer(kind=prec_int), intent(in) :: n
  character(LEN=5), intent(out) :: nchar
  character(LEN=1) :: nchar1
  character(LEN=2) :: nchar2
  character(LEN=3) :: nchar3
  character(LEN=4) :: nchar4
  character(LEN=5) :: nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif
end subroutine title

subroutine timer_start
  use hydro_mpi
  implicit none

  !Write buffers to disk before writing files
  code = sync()

  call MPI_Barrier(MPI_COMM_WORLD,code)
  tb = MPI_Wtime()
end subroutine timer_start


subroutine timer_stop(filetype,writing,datasize)
  use hydro_mpi
  implicit none

  character(len=*)::filetype
  logical :: writing
  integer(kind=8)::datasize

  real(kind=8),parameter :: MiB=1024.d0*1024.d0, GiB=1024.d0*MiB, TiB=1024.d0*GiB
  real(kind=8) :: te, timeloc, tmin, tmax, throughput
  real(kind=8) :: timewithsync, tsync, tminsync, tmaxsync, throughputwithsync
  character(len=11) :: sizetext
  character(len=13) :: throughputtext ,throughputwithsynctext

  te = MPI_Wtime()
  timeloc = te-tb

  !Write buffers to disk (warning: all files are concerned)
  code = sync()
  tsync = MPI_Wtime()
  timewithsync = tsync-tb

  call MPI_Reduce(timeloc,tmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,code)
  call MPI_Reduce(timeloc,tmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,code)

  call MPI_Reduce(timewithsync,tminsync,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,code)
  call MPI_Reduce(timewithsync,tmaxsync,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,code)

  if (datasize<GiB) then
    write(sizetext,'(F8.3,"MiB")') datasize/MiB
  else if (datasize<TiB) then
    write(sizetext,'(F8.3,"GiB")') datasize/GiB
  else
    write(sizetext,'(F8.3,"TiB")') datasize/TiB
  end if

  throughput = datasize/tmax
  if (throughput<GiB) then
    write(throughputtext,'(F8.3,"MiB/s")') throughput/MiB
  else if (throughput<TiB) then
    write(throughputtext,'(F8.3,"GiB/s")') throughput/GiB
  else
    write(throughputtext,'(F8.3,"TiB/s")') throughput/TiB
  end if

  throughputwithsync = datasize/tmaxsync
  if (throughputwithsync<GiB) then
    write(throughputwithsynctext,'(F8.3,"MiB/s")') throughputwithsync/MiB
  else if (throughputwithsync<TiB) then
    write(throughputwithsynctext,'(F8.3,"GiB/s")') throughputwithsync/GiB
  else
    write(throughputwithsynctext,'(F8.3,"TiB/s")') throughputwithsync/TiB
  end if

  if (rang==0) then
    if(writing) then
      write(*,'(A," : Time to write:             max=",F9.3,"s min=",F9.3,"s")') trim(filetype),tmax,tmin
      write(*,'(A," : Time to write (with sync): max=",F9.3,"s min=",F9.3,"s")') trim(filetype),tmaxsync,tminsync
      print *,trim(filetype)," : Write throughput:             ",trim(throughputtext)," (data written: ",trim(sizetext),")"
      print *,trim(filetype)," : Write throughput (with sync): ",trim(throughputwithsynctext)," (data written: ",trim(sizetext),")"
    else
      write(*,'(A," : Time to read max=",F9.3,"s min=",F9.3,"s")') trim(filetype),tmax,tmin
      print *,trim(filetype)," : Read throughput: ",trim(throughputtext)," (data read: ",trim(sizetext),")"
    end if
  end if
end subroutine timer_stop

end module hydro_utils
