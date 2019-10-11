module m_da_maooam
  use m_mt, only : eigen, inv
  implicit none

  private 
  public :: etkf
  public :: etkf2

  !integer,parameter :: nx = 20
  !integer,parameter :: nyo = 20
  !integer,parameter :: nx = 36
  !integer,parameter :: nyo = 36
  !integer,parameter :: kitermax = 500
  real(8),parameter :: epsG   = 1.d-2

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


subroutine etkf2( nx, nyo, nn, xb, lyo, yo, erro, infl, xa, xam )
  implicit none
! passed args
  integer, intent(in) :: nx
  integer, intent(in) :: nyo
  integer,    intent(in   ) :: nn
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
  real(8) :: xbm(nx)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn)
  real(8) :: dxb(nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
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
     dxb = xb(i,:)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(:)*Wa(:,k))
     enddo
  enddo

endsubroutine




endmodule
