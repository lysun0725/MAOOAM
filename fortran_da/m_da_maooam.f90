module m_da_maooam
  use m_mt, only : eigen, inv, shuffle
  use m_mt, only : randn 
  implicit none

  public :: etkf
  public :: etkf4d
  public :: etkf_w, pf, pf2, MH_resampling
  public :: enkf_w, RSR_resampling_1D

  !integer,parameter :: nx = 20
  !integer,parameter :: nyo = 20
  !integer,parameter :: nx = 36
  !integer,parameter :: nyo = 36
  !integer,parameter :: kitermax = 500
  real(8),parameter :: epsG   = 1.d-2

  interface MH_resampling
    module procedure MH_resampling_1D
    !module procedure MH_resampling_1D_2
    module procedure MH_resampling_2D
  end interface

contains

subroutine etkf( nx, nyo, nn, xb, lyo, yo, erro, infl, xa, xam, fnum )
  implicit none
! passed args
  integer,    intent(in   ) :: nx	! dim of model space
  integer,    intent(in   ) :: nyo	! dim of obs space
  integer,    intent(in   ) :: nn	! dim of ens size
  integer,    intent(in   ) :: fnum
  real(8), intent(in   ) :: xb(nx,nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo)
  real(8), intent(in   ) :: erro(nyo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  integer :: nyo_use
  real(8),allocatable :: dyb(:,:)
  real(8),allocatable :: ybm(:)
  real(8),allocatable :: erro_use(:)
  real(8),allocatable :: yo_ybm(:)
  real(8),allocatable :: gain(:,:)
  real(8) :: xbm(nx)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn)
  real(8) :: dxb(nx,nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( dyb(nyo_use,nn) )
  allocate( ybm(nyo_use), yo_ybm(nyo_use), erro_use(nyo_use) )

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        dyb(k,:)    = xb(i,:)
        ybm(k)      = SUM(dyb(k,:))/nn
        dyb(k,:)    = dyb(k,:) - ybm(k)
        yo_ybm(k)   = yo(i) - ybm(k)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:))/nn
  enddo
 
  allocate( C(nn,nyo_use) )

! C=(Yb)^T * R^-1
  do j = 1, nyo_use
     C(:,j) = dyb(j,:)/(erro_use(j)**2)
  enddo
! invSPa = [(k-1)I+C*Yb]
  invSPa = MATMUL( C, dyb )
  do k = 1, nn
     invSPa(k,k) = invSPa(k,k)+(nn-1)/infl
  enddo
! let A=invSPa, where A=Q*V*Q^T
! so A^-1=Q*V^-1*Q^T
  call eigen( nn, invSPa, eigvect, eigval)
  !if ( ierr/=0 .or. np/=nn ) then
  !   write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
  !   xam = xbm
  !   xa = xb
  !   return
  !endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  ! wam = SPa * C *( yo -ybm )
  wam = MATMUL( C, yo_ybm )
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo

  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb(i,:) = xb(i,:)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(i,:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(i,:)*Wa(:,k))
     enddo
  enddo

  ! write gain matrix
  allocate(gain(nx,nyo))
  gain = MATMUL(SPa,C)
  gain = MATMUL(dxb,gain)
  !do i = 1,nx
  !  write(fnum,*) gain(i,:)
  !enddo

endsubroutine

subroutine etkf4d( nx, nyo, ntda, nn, xb, lyo, yo, erro, infl, xa, xam )
  implicit none
! passed args
  integer, intent(in) :: nx
  integer, intent(in) :: nyo
  integer, intent(in) :: ntda  ! Number of obs within the DA window
  integer,    intent(in   ) :: nn
  real(8), intent(in   ) :: xb(nx,nn,ntda)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo,ntda)
  real(8), intent(in   ) :: erro(nyo,ntda)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  integer :: nyo_use
  real(8),allocatable :: dyb(:,:,:)
  real(8),allocatable :: ybm(:,:)
  real(8),allocatable :: erro_use(:,:)
  real(8),allocatable :: yo_ybm(:,:)
  real(8) :: xbm(nx)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn)
  real(8) :: dxb(nn)
  integer :: ierr
  integer :: i, j, k, l, m

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i,ntda)
     endif
  enddo
  allocate( dyb(nyo_use,nn,ntda) )
  allocate( ybm(nyo_use,ntda), yo_ybm(nyo_use,ntda), erro_use(nyo_use,ntda) )

  do l = 1,ntda ! loop within the DA window
    k = 0
    do i = 1, nyo
       if ( lyo(i) ) then
          k           = k + 1
          dyb(k,:,l)    = xb(i,:,l)
          ybm(k,l)      = SUM(dyb(k,:,l))/nn
          dyb(k,:,l)    = dyb(k,:,l) - ybm(k,l)
          yo_ybm(k,l)   = yo(i,l) - ybm(k,l)
          erro_use(k,l) = erro(i,ntda)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
       endif
    enddo
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:,ntda))/nn
  enddo
 
  allocate( C(nn,nyo_use,ntda) )

! C=(Yb)^T * R^-1
  do l = 1, ntda
    do j = 1, nyo_use
      C(:,j,l) = dyb(j,:,l)/(erro_use(j,l)**2)
    enddo
  enddo
! invSPa = [(k-1)I+ \sum C*Yb]
  invSPa = 0.d0
  do l = 1, ntda
    invSPa = invSPa + MATMUL( C(:,:,l), dyb(:,:,l) )
  enddo
  do k = 1, nn
    invSPa(k,k) = invSPa(k,k)+(nn-1)/infl
  enddo

! let A=invSPa, where A=Q*V*Q^T
! so A^-1=Q*V^-1*Q^T
  call eigen( nn, invSPa, eigvect, eigval)
  !if ( ierr/=0 .or. np/=nn ) then
  !   write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
  !   xam = xbm
  !   xa = xb
  !   return
  !endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  
! wam = SPa * [sum C *( yo -ybm )]
  wam = 0.d0
  do l = 1, ntda
    wam = wam + MATMUL( C(:,:,l), yo_ybm(:,l) )
  enddo
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo

  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb = xb(i,:,ntda)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(:)*Wa(:,k))
     enddo
  enddo

endsubroutine

subroutine etkf_w( nx, nyo, nn, xb, w, lyo, yo, erro, infl, xa, xam, fnum )
  implicit none
! passed args
  integer,    intent(in   ) :: nx       ! dim of model space
  integer,    intent(in   ) :: nyo      ! dim of obs space
  integer,    intent(in   ) :: nn       ! dim of ens size
  integer,    intent(in   ) :: fnum
  real(8), intent(in   ) :: xb(nx,nn)
  real(8), intent(in   ) :: w(nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo)
  real(8), intent(in   ) :: erro(nyo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  integer :: nyo_use
  real(8),allocatable :: dyb(:,:)
  real(8),allocatable :: ybm(:)
  real(8),allocatable :: erro_use(:)
  real(8),allocatable :: yo_ybm(:)
  real(8),allocatable :: gain(:,:)
  real(8) :: xbm(nx)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn), neff
  real(8) :: dxb(nx,nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
  neff = 1/SUM(w*w) 
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( dyb(nyo_use,nn) )
  allocate( ybm(nyo_use), yo_ybm(nyo_use), erro_use(nyo_use) )

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        dyb(k,:)    = xb(i,:)
        ybm(k)      = SUM(dyb(k,:)*w)
        dyb(k,:)    = (dyb(k,:) - ybm(k))*SQRT(w*neff) 
        yo_ybm(k)   = yo(i) - ybm(k)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:)*w)
  enddo

  allocate( C(nn,nyo_use) )

! C=(Yb)^T * R^-1
  do j = 1, nyo_use
     C(:,j) = dyb(j,:)/(erro_use(j)**2)
  enddo
! invSPa = [(k-1)I+C*Yb]
  invSPa = MATMUL( C, dyb )
  do k = 1, nn
     invSPa(k,k) = invSPa(k,k)+(neff)/infl
  enddo
! let A=invSPa, where A=Q*V*Q^T
! so A^-1=Q*V^-1*Q^T
  call eigen( nn, invSPa, eigvect, eigval)
  !if ( ierr/=0 .or. np/=nn ) then
  !   write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
  !   xam = xbm
  !   xa = xb
  !   return
  !endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     !Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
     Wa(:,k) = eigvect(:,k)*SQRT( (neff)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  ! wam = SPa * C *( yo -ybm )
  wam = MATMUL( C, yo_ybm )
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo

  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb(i,:) = (xb(i,:)-xbm(i))*SQRT(w*neff)
     xam(i) = xbm(i) + SUM(dxb(i,:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(i,:)*Wa(:,k))
     enddo
  enddo
 
  ! write gain matrix
  allocate(gain(nx,nyo))
  gain = MATMUL(SPa,C)
  gain = MATMUL(dxb,gain)
  !do i = 1,nx
  !  write(fnum,*) gain(i,:)
  !enddo
endsubroutine

subroutine enkf_w( nx, nyo, nn, xb, w, lyo, yo, erro, infl, xa, xam, fnum )
  implicit none
! passed args
  integer,    intent(in   ) :: nx       ! dim of model space
  integer,    intent(in   ) :: nyo      ! dim of obs space
  integer,    intent(in   ) :: nn       ! dim of ens size
  integer,    intent(in   ) :: fnum
  real(8), intent(in   ) :: xb(nx,nn)
  real(8), intent(in   ) :: w(nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo)
  real(8), intent(in   ) :: erro(nyo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local var
  integer :: nyo_use
  real(8) :: Z(nx,nn)
  real(8) :: xbm(nx)
  real(8), allocatable :: Hz(:,:), yo_yb(:,:), erro_use(:)
  real(8), allocatable :: gain(:,:), eigvec(:,:), eigval(:), C(:,:)
  real(8), allocatable :: eta(:,:), d(:,:), Cinv(:,:)
  integer :: i,j,k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  
  allocate(Hz(nyo_use,nn), yo_yb(nyo_use,nn), erro_use(nyo_use))
  allocate(gain(nx,nyo_use),eigvec(nyo_use,nyo_use), eigval(nyo_use))
  allocate(C(nyo_use,nyo_use), eta(nyo_use,nn), d(nyo_use,nn))
  allocate(Cinv(nyo_use,nyo_use))
  xbm = MATMUL(xb,w)


  ! matrix of difference
  do i = 1, nx
    Z(i,:) = (xb(i,:) - xbm(i)) * sqrt(w)
  enddo 

  ! Kalman gain matrix
  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        Hz(k,:)     = Z(i,:)
        yo_yb(k,:)   = yo(i) - xb(i,:)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, Hz =", k, erro_use(k), yo(i), Hz(k,:)
        !pause
     endif
  enddo
  ! compute C=inv((Hz*Hz')+R)
  Cinv = MATMUL(Hz, TRANSPOSE(Hz))
  do k = 1, nyo_use 
     Cinv(k,k) = Cinv(k,k)+erro_use(k)**2/infl
  enddo
  call eigen( nyo_use, Cinv, eigvec, eigval) 
  do i = 1, nyo_use
     C(:,i) = eigvec(:,i)/eigval(i)
  enddo
 
  ! compute K = (Z*Hz')*inv((Hz*Hz')+R)
  gain = MATMUL(MATMUL(Z,TRANSPOSE(Hz)), C)

  ! inovation vectors
  ! define eta
  do j = 1, nn
    call randn(nyo_use,eta(:,j))
    eta(:,j) = eta(:,j) * erro_use * 0.1
  enddo
  d = yo_yb + eta
  xa = xb + MATMUL(gain,d)
  xam = SUM(xa,DIM=2)/nn
endsubroutine

subroutine pf(nyo, nx, nn, npart, xpart_b, w, lyo, yo, erro, w_a)
  implicit none
  integer,    intent(in   ) :: nyo      ! dim of obs space
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: nn       ! dim of ens size
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_b(nx, nn, npart)
  real(8),    intent(in   ) :: w(nn, npart)
  logical,    intent(in   ) :: lyo(nyo)
  real(8),    intent(in   ) :: yo(nyo)
  real(8),    intent(in   ) :: erro(nyo)
  real(8),    intent(out  ) :: w_a(nn, npart)
! local vars
  integer :: nyo_use
  real(8),allocatable :: erro_use(:)
  real(8),allocatable :: yo_use(:)
  real(8),allocatable :: hxpart_b(:,:,:)
  real(8),allocatable :: PFexp(:,:), cond_prob(:,:)
  real(8),allocatable :: yo_part(:)
  integer :: i, j, k, l, m
  real(8) :: infl

  allocate( PFexp(nn,npart), cond_prob(nn,npart) )
  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( yo_use(nyo_use), erro_use(nyo_use), hxpart_b(nyo_use,nn,npart) )
  allocate( yo_part(nyo_use) )  

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        yo_use(k)   = yo(i)
        erro_use(k) = erro(i)
        if (nyo > nx) then
          hxpart_b(k,:,:) = xpart_b(k,:,:)
        elseif (nyo == nx) then
          hxpart_b(k,:,:) = xpart_b(i,:,:)
        endif
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

  do j = 1, nn
    !$OMP PARALLEL DO
    do l = 1, npart
      ! PFexp(j,l) = 1/2 * (yo-hxpart)' * Rinv * (yo-hxpart) and Rinv = 1/(error_use**2)
      yo_part = (yo_use - hxpart_b(:,j,l))
      !IF (j==1 .AND. l==1)  PRINT *, "yo_part = ", yo_part
      PFexp(j,l) = 0.5 * SUM(yo_part*yo_part/erro_use/erro_use)
      cond_prob(j,l) = w(j,l) * EXP(-PFexp(j,l)) 
    enddo
    !$OMP END PARALLEL DO
  enddo

  infl = MINVAL(PFexp)
  if ( SUM(cond_prob) < 1.D-150 ) then
    write(6,*) "cond_prob=0, inflating: infl = ", infl
    do j = 1, nn
      !$OMP PARALLEL DO
      do l = 1, npart
        cond_prob(j,l) = w(j,l) * EXP(-PFexp(j,l) + infl)
      enddo
      !$OMP END PARALLEL DO
    enddo
  endif

  w_a = cond_prob / SUM(cond_prob)
  PRINT *,"DEBUG_PF: SUM(cond_prob) = ", SUM(cond_prob)

  deallocate(yo_use,erro_use,hxpart_b,yo_part,PFexp,cond_prob)
endsubroutine

subroutine pf2(nyo, nx, nn, npart, xpart_b, w, xam, lyo, yo, erro, std,  w_a)
  implicit none
  integer,    intent(in   ) :: nyo      ! dim of obs space
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: nn       ! dim of ens size
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_b(nx, nn, npart)
  real(8),    intent(in   ) :: xam(nx)
  real(8),    intent(in   ) :: w(nn, npart)
  logical,    intent(in   ) :: lyo(nyo)
  real(8),    intent(in   ) :: yo(nyo)
  real(8),    intent(in   ) :: erro(nyo)
  real(8),    intent(in   ) :: std(nyo)
  real(8),    intent(out  ) :: w_a(nn, npart)
! local vars
  integer :: nyo_use
  real(8),allocatable :: erro_use(:), std_use(:)
  real(8),allocatable :: yo_use(:)
  real(8),allocatable :: hxpart_b(:,:,:)
  real(8),allocatable :: PFexp(:,:), PFexp2(:,:), cond_prob(:,:)
  real(8),allocatable :: yo_part(:), mu_part(:)
  integer :: i, j, k, l, m
  real(8) :: infl

  allocate( PFexp(nn,npart), PFexp2(nn,npart), cond_prob(nn,npart) )
  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        !write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( yo_use(nyo_use), erro_use(nyo_use), hxpart_b(nyo_use,nn,npart) )
  allocate( yo_part(nyo_use), std_use(nyo_use), mu_part(nyo_use))

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        yo_use(k)   = yo(i)
        erro_use(k) = erro(i)
        std_use(k)  = std(i)
        if (nyo > nx) then
          hxpart_b(k,:,:) = xpart_b(k,:,:)
        elseif (nyo == nx) then
          hxpart_b(k,:,:) = xpart_b(i,:,:)
        endif
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

  do j = 1, nn
    !$OMP PARALLEL DO
    do l = 1, npart
      ! PFexp(j,l) = 1/2 * (yo-hxpart)' * Rinv * (yo-hxpart) and Rinv = 1/(error_use**2)
      yo_part = (yo_use - hxpart_b(:,j,l))
      mu_part = (hxpart_b(:,j,l) - xam)
      !IF (j==1 .AND. l==1)  PRINT *, "mu_part = ", mu_part
      PFexp(j,l) = 0.5 * SUM(yo_part*yo_part/erro_use/erro_use)
      !IF (j==1 .AND. l==1)  PRINT *, "PFexp = ", PFexp(j,l)
      PFexp2(j,l) = 0.5 * SUM(mu_part*mu_part/std_use/std_use)
      !IF (j==1 .AND. l==1)  PRINT *, "delta_PFexp = ", 0.5 * SUM(mu_part*mu_part/std_use/std_use)
      cond_prob(j,l) = w(j,l) * EXP(-PFexp(j,l))/EXP(-PFexp2(j,l))
    enddo
    !$OMP END PARALLEL DO
  enddo

  infl = MINVAL(PFexp-PFexp2)
  if ( SUM(cond_prob) < 1.D-150 ) then
    write(6,*) "cond_prob=0, inflating: infl = ", infl
    do j = 1, nn
      !$OMP PARALLEL DO
      do l = 1, npart
        cond_prob(j,l) = w(j,l) * EXP(-PFexp(j,l) + infl)
      enddo
      !$OMP END PARALLEL DO
    enddo
  endif

  w_a = cond_prob / SUM(cond_prob)
  PRINT *,"DEBUG_PF: SUM(cond_prob) = ", SUM(cond_prob)

  deallocate(yo_use,erro_use,hxpart_b,yo_part,PFexp,cond_prob)
endsubroutine


! resamples via Metropolis-Hastings
subroutine MH_resampling_1D_2(nx, npart, xpart_old, w_old, xpart_new, w_new)
  implicit none
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_old(nx,npart)
  real(8),    intent(in   ) :: w_old(npart)
  real(8),    intent(out  ) :: xpart_new(nx,npart)
  real(8),    intent(out  ) :: w_new(npart)
! local params
  integer :: b(npart), nn, Nfirst
  real(8) :: r
  integer :: i, j, k, l, m, n, cnt
  integer, allocatable :: seed(:)

  call system_clock(cnt)
  call random_seed(size=n)
  allocate(seed(n), source=cnt+37*[(i,i=0,n-1)])
  call random_seed(put=seed)

  b = (/ (nn, nn = 1,npart) /)
  xpart_new=0.d0
  r=0.d0

  do i = 2, npart
    if ( w_old(b(i)) > w_old(b(i-1)) ) then
      b(i) = b(i)
    else
      call random_number(r)
      if ( r < w_old(b(i))/w_old(b(i-1)) ) then
        b(i) = b(i)
      else
        b(i) = b(i-1)
      end if
    endif
  enddo

  !$OMP PARALLEL DO
  do k = 1, npart
    xpart_new(:,k) = xpart_old(:,b(k))
  enddo
  !$OMP END PARALLEL DO

  w_new = 1.d0/npart

endsubroutine

! resamples via Metropolis-Hastings
subroutine RSR_resampling_1D(nx, npart, xpart_old, w_old, xpart_new, w_new)
  implicit none
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_old(nx,npart)
  real(8),    intent(in   ) :: w_old(npart)
  real(8),    intent(out  ) :: xpart_new(nx,npart)
  real(8),    intent(out  ) :: w_new(npart)
! local params
  integer :: b(npart), nn, Nfirst
  real(8) :: r
  integer :: i,ii, j, jj, k, l, m, n, cnt
  integer, allocatable :: seed(:)

  call system_clock(cnt)
  call random_seed(size=n)
  allocate(seed(n), source=cnt+37*[(i,i=0,n-1)])
  call random_seed(put=seed)

  jj=0; call random_number(r); r=r/npart

  do i = 1, npart
    j = floor( (w_old(i)-r)*npart ) + 1
    do ii = 1, j
      jj=jj+1
      xpart_new(:,jj)=xpart_old(:,i);
    enddo
    r=r+real(j,kind=8)/real(npart,kind=8)-w_old(i)
  enddo

  w_new = 1.d0/npart
  deallocate(seed)
endsubroutine


! resamples via Metropolis-Hastings
subroutine MH_resampling_1D(nx, npart, xpart_old, w_old, xpart_new, w_new)
  implicit none
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_old(nx,npart)
  real(8),    intent(in   ) :: w_old(npart)
  real(8),    intent(out  ) :: xpart_new(nx,npart)
  real(8),    intent(out  ) :: w_new(npart)
! local params
  integer :: b(npart), nn, Nfirst
  real(8) :: r
  integer :: i, j, k, l, m, n, cnt
  integer, allocatable :: seed(:)

  call system_clock(cnt)
  call random_seed(size=n)
  allocate(seed(n), source=cnt+37*[(i,i=0,n-1)])
  call random_seed(put=seed)
  b = (/ (nn, nn = 1,npart) /); call shuffle(b)
  !print *, "debug: ", b
  xpart_new=0.d0
  Nfirst=max(nint(npart/2.0),1)
  j=npart; r=0.d0

  do i = 2, npart+Nfirst
    if (i>npart) then
      j=j-1
      if ( w_old(b(i-npart)) >= w_old(b(j)) ) then
        b(i-npart) = b(i-npart)
      else
        ! accept particle i with probability W(i)/W(b(i-1))
        call random_number(r) !; print *, "debug: ", r
        if ( r < w_old(b(i-npart))/w_old(b(j)) ) then
           !b(1) = 1
        else
          b(i-npart) = b(j)
        endif
      endif
    else
      if ( w_old(b(i)) >= w_old(b(i-1)) ) then !b(i-1) is index of previously accepted particle
        b(i) = b(i)
      else
        ! accept particle i with probability W(i)/W(b(i-1))
        call random_number(r) !; print *, "debug: ", r
        if ( r < w_old(b(i))/w_old(b(i-1)) ) then
          b(i) = b(i)
        else
          b(i) = b(i-1)
        endif
      endif
    endif
  enddo

  !$OMP PARALLEL DO
  do k = 1, npart
    xpart_new(:,k) = xpart_old(:,b(k))
  enddo
  !$OMP END PARALLEL DO

  w_new = 1.d0/npart
  deallocate(seed)

endsubroutine

subroutine MH_resampling_2D(nx, nn, npart, xpart_old, w_old, xpart_new, w_new)
  implicit none
  integer,    intent(in   ) :: nx       ! dim of each particle 
  integer,    intent(in   ) :: nn       ! num of ensemble members 
  integer,    intent(in   ) :: npart    ! num of particles 
  real(8),    intent(in   ) :: xpart_old(nx,nn,npart)
  real(8),    intent(in   ) :: w_old(nn,npart)
  real(8),    intent(out  ) :: xpart_new(nx,nn,npart)
  real(8),    intent(out  ) :: w_new(nn,npart)
! local params
  real(8), allocatable :: xpart_old_long(:,:), w_old_long(:)
  integer :: b(npart*nn), Nfirst
  real(8) :: r
  integer :: i, j, k, l, m, n, cnt
  integer, allocatable :: seed(:)

  call system_clock(cnt)
  call random_seed(size=n)
  allocate(seed(n), source=cnt+37*[(i,i=0,n-1)])
  call random_seed(put=seed)

  ! reshape the xpart_old and w_old
  allocate(xpart_old_long(nx,nn*npart),w_old_long(nn*npart))
  do k = 1, nn
    !$OMP PARALLEL DO
    do l = 1, npart
      xpart_old_long(:,(k-1)*npart+l) = xpart_old(:,k,l)
      w_old_long((k-1)*npart+l) = w_old(k,l) 
    enddo
    !$OMP END PARALLEL DO
  enddo

  b = (/ (m, m = 1,npart*nn) /); call shuffle(b)
  xpart_new=0.d0
  Nfirst=max(nint(npart*nn/25.0),1)
  j=npart*nn; r=0.d0

  do i = 2, npart*nn+Nfirst
    if (i>npart*nn) then
      j=j-1
      if ( w_old_long(b(i-npart*nn)) >= w_old_long(b(j)) ) then
        b(i-npart*nn) = b(i-npart*nn)
      else
        ! accept particle i with probability W(i)/W(b(i-1))
        call random_number(r)
        if ( r < w_old_long(b(i-npart*nn))/w_old_long(b(j)) ) then
           !b(1) = 1
        else
          b(i-npart*nn) = b(j)
        endif
      endif
    else
      if ( w_old_long(b(i)) >= w_old_long(b(i-1)) ) then !b(i-1) is index of previously accepted particle
        b(i) = b(i)
      else
        ! accept particle i with probability W(i)/W(b(i-1))
        call random_number(r)
        if ( r < w_old_long(b(i))/w_old_long(b(i-1)) ) then
          b(i) = b(i)
        else
          b(i) = b(i-1);
        endif
      endif
    endif
  enddo

  do k = 1, nn
    !$OMP PARALLEL DO
    do l = 1, npart
      xpart_new(:,k,l) = xpart_old_long(:,b((k-1)*npart+l))
    enddo
    !$OMP END PARALLEL DO
  enddo

  w_new = 1.d0/npart/nn


  deallocate(xpart_old_long,w_old_long)
endsubroutine


endmodule
