module m_da_maooam_phys
  use m_mt, only : eigen, inv
  use proj_spec_phys, only: compute_ocn_uvel, compute_ocn_vvel
  use params, only : nobs_dr,tw_da
  use params, only : oms, ndim, natm, noc,n,dr_num 
  implicit none

  public :: read_dimension_phys, read_dimension_drf
  public :: read_obs_phys
  public :: etkf_phys,etkf4d_phys

  real(8),parameter :: epsG   = 1.d-2

contains

subroutine obsop(nx,nn,nyo,ncol_yo,xb,yo,dyb,yo_ybm,error)
  use params, only: obs_atm_u, obs_atm_v, obs_atm_t
  use params, only: obs_ocn_u, obs_ocn_v, obs_ocn_t
  use params, only: obs_drf_x, obs_drf_y, obs_drf_t
  implicit none

  integer,    intent(in   ) :: nx
  integer,    intent(in   ) :: nn
  integer,    intent(in   ) :: nyo
  integer,    intent(in   ) :: ncol_yo    ! columns of yobs_*_phys.dat file.
  real(8),    intent(in   ) :: xb(nx,nn)
  real(8),    intent(in   ) :: yo(nyo,ncol_yo)
  real(8),    intent(out  ) :: dyb(nyo,nn)
  real(8),    intent(out  ) :: yo_ybm(nyo)
  real(8),    intent(out  ) :: error(nyo)
  integer :: i,j,k,l,m, Ho, Po
  real(8) :: oelem, opos(2), R(nx), oerr_ou, oerr_ov 
  real(8) :: oerr_dx, oerr_dy
  real(8),    allocatable :: hxb(:,:)
  real(8),    allocatable :: ybm(:)

  allocate(hxb(nyo,nn)); hxb = 0.d0
  allocate(ybm(nyo)); ybm = 0.d0; yo_ybm = 0.d0
  oerr_ou=0.d0; oerr_ov=0.d0; oerr_dx=0.d0; oerr_dy=0.d0

  ! Example of yobs_drf_phys.dat: ncol_yo = 4
  ! oelem   olon    olat    odata
  open(103,file="fort.202",action="read",form="formatted",access="sequential")
  R = 0.d0
  do j = 1, ndim
    read(103,"(10000(D24.17,1x))") R(j)
  enddo
  close(103)
  !print*, "debug R=", R
  !print*, "debug oerr_ou=", oerr_ou

  !compute observation error for ocn u and v, atms u and v are TBD
  do m = 1,noc
    Ho=oms(m,1); Po=oms(m,2)
    oerr_ou = oerr_ou + (R(m+2*natm)**2)*(Po**2)/noc
    oerr_ov = oerr_ov + (R(m+2*natm)**2)*((Ho*n/2)**2)/noc
  enddo

  oerr_ou=SQRT(oerr_ou); oerr_ov=SQRT(oerr_ov)
  !print*, "oerr_ou = ", oerr_ou, "oerr_ov = ", oerr_ov
  oerr_dx=oerr_ou*tw_da*100; oerr_dy = oerr_ov*tw_da*100

  do i = 1, nyo
    ! Read the position
    oelem = yo(i,1); opos(1) = yo(i,2); opos(2) = yo(i,3)
    ! print *, "i = ", i, "oelem = ", oelem, "opos = ", opos
    ! For each ensemble member
    do j = 1, nn
      select case (INT(oelem))
        case (obs_ocn_u)
          call compute_ocn_uvel(xb(1:nx,j),opos,hxb(i,j))
          error(i)=oerr_ou
        case (obs_ocn_v) 
          call compute_ocn_vvel(xb(1:nx,j),opos,hxb(i,j))
          error(i)=oerr_ov
        case (obs_drf_x)
          error(i)=oerr_dx
        case (obs_drf_y)
          error(i)=oerr_dy
      end select
      !print*, "j = ", j, "hxb=",hxb(i,j),"oerror=", error(i)
    end do

    ! Compute yo_ybm, dyb and ybm
    dyb(i,:) = hxb(i,:)
    ybm(i) = SUM(dyb(i,:))/nn
    dyb(i,:) = hxb(i,:) - ybm(i)
    yo_ybm(i) = yo(i,4) - ybm(i)
  end do

  
endsubroutine

subroutine etkf_phys( nx, nyo, nn, ncol_yo, xb, yo, infl, xa, xam, fnum )
  implicit none
! passed args
  integer,    intent(in   ) :: nx	! dim of model space
  integer,    intent(in   ) :: nyo	! dim of obs space
  integer,    intent(in   ) :: nn	! dim of ens size
  integer,    intent(in   ) :: ncol_yo	! dim of yo
  integer,    intent(in   ) :: fnum
  real(8), intent(in   ) :: xb(nx,nn)
  real(8), intent(in   ) :: yo(nyo,ncol_yo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  real(8),allocatable :: dyb(:,:)
  real(8),allocatable :: yo_ybm(:), error(:)
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

  allocate( dyb(nyo,nn) )
  allocate( yo_ybm(nyo), error(nyo) ) ! LUYU TBD: replace erro

  call obsop(nx,nn,nyo,ncol_yo,xb,yo,dyb,yo_ybm,error)
  ! print *, "debug: dyb(1,:) = ", dyb(1,:)
  ! print *, "debug: yo(1,4) = ", yo(1,4)
  ! print *, "debug: yo_ybm(1) = ", yo_ybm(1)
  ! print *, "debug: error(last_two) = ", error(nyo-1:nyo)

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:))/nn
  enddo
 
  allocate( C(nn,nyo) )

! C=(Yb)^T * R^-1
  do j = 1, nyo
     C(:,j) = dyb(j,:)/(error(j)**2)
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


subroutine etkf4d_phys( nx, nyo, ntda, nn, ncol_yo, xb, yo, infl, xa, xam )
  implicit none
! passed args
  integer, intent(in) :: nx
  integer, intent(in) :: nyo
  integer, intent(in) :: ntda  ! Number of obs within the DA window
  integer, intent(in) :: nn
  integer, intent(in) :: ncol_yo
  real(8), intent(in   ) :: xb(nx,nn,ntda)
  real(8), intent(in   ) :: yo(nyo,ncol_yo,ntda)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  real(8),allocatable :: dyb(:,:,:)
  real(8),allocatable :: ybm(:,:)
  real(8),allocatable :: error(:,:)
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

  allocate( dyb(nyo,nn,ntda) )
  allocate( ybm(nyo,ntda), yo_ybm(nyo,ntda), error(nyo,ntda) )

  do l = 1, ntda
    call obsop(nx,nn,nyo,ncol_yo,xb(:,:,l),yo(:,:,l),dyb(:,:,l),yo_ybm(:,l),error(:,l))
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:,ntda))/nn
  enddo
 
  allocate( C(nn,nyo,ntda) )

! C=(Yb)^T * R^-1
  do l = 1, ntda
    do j = 1, nyo
      C(:,j,l) = dyb(j,:,l)/(error(j,l)**2)
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

subroutine pf_phys
endsubroutine

subroutine hybrid_enkf_pf_phys()
endsubroutine


subroutine read_dimension_phys(do_atm_obs, do_ocn_obs, do_drf_vel, cnt_obs_phys)
  implicit none
  logical, intent(in) :: do_atm_obs
  logical, intent(in) :: do_ocn_obs
  logical, intent(in) :: do_drf_vel
  integer, intent(inout) :: cnt_obs_phys
  integer :: io
  real(8) :: wk(5),t_old, t_new 
  character(len=5) :: pathname

  io=0
  WRITE(pathname,'("d",I0.3,"/")') dr_num
  if (do_atm_obs) then
    open(801,file=pathname//'yobs_atm_phys.dat',action="read",access="sequential")
    read(801,*,IOSTAT=io) wk
    if (io==0) cnt_obs_phys = cnt_obs_phys+1; t_old = wk(1) 

    do while (io == 0)
      read(801,*,IOSTAT=io) wk 
      t_new = wk(1)
      if (io==0  .and. t_new == t_old) then
        cnt_obs_phys = cnt_obs_phys+1; t_old = t_new
      else
        exit
      endif
    enddo    
    close(801)
  endif
  print*, "read_dim_phys: after atm_phys cnt_obs = ", cnt_obs_phys

  if (do_ocn_obs) then
    open(802,file=pathname//'yobs_ocn_phys.dat',action="read",access="sequential")
    read(802,*,IOSTAT=io) wk  
    if (io==0) cnt_obs_phys = cnt_obs_phys+1; t_old = wk(1)
    
    do while (io == 0)
      read(802,*,IOSTAT=io) wk
      t_new = wk(1)
      if (io==0  .and. t_new == t_old) then
        cnt_obs_phys = cnt_obs_phys+1; t_old = t_new
      else
        exit
      endif
    enddo
    close(802)
  endif
  print*, "read_dim_phys: after ocn_phys cnt_obs = ", cnt_obs_phys
 
  if (do_drf_vel) then
    open(803,file=pathname//'yobs_drf_vel.dat',action="read",access="sequential")
    read(803,*, IOSTAT=io) wk
    if (io==0) cnt_obs_phys = cnt_obs_phys+1; t_old = wk(1)
    
    do while (io == 0)
      read(803,*, IOSTAT=io) wk
      t_new = wk(1)
      if (io==0  .and. t_new == t_old) then
        cnt_obs_phys = cnt_obs_phys+1; t_old = t_new
      else
        exit
      endif
    enddo
    close(803)
  endif
  print*, "read_dim_phys: after drf_velo cnt_obs = ", cnt_obs_phys
 
endsubroutine

subroutine read_dimension_drf(do_drf_pos,cnt_obs_drf)
  implicit none
  logical, intent(in) :: do_drf_pos
  integer, intent(inout) :: cnt_obs_drf

  if (do_drf_pos) then
    cnt_obs_drf = nobs_dr
  endif
endsubroutine

subroutine read_obs_phys(do_atm_obs, do_ocn_obs, do_drf_vel, cnt_obs_phys, nt, yobs)
  implicit none
  logical, intent(in) :: do_atm_obs
  logical, intent(in) :: do_ocn_obs
  logical, intent(in) :: do_drf_vel
  integer, intent(in) :: cnt_obs_phys
  integer, intent(in) :: nt
  real(8), intent(inout) :: yobs(cnt_obs_phys,4,nt)
  real(8) :: wk(5)
  integer :: i,j,k,io
  integer :: cnt_atm, cnt_ocn
  character(len=5) :: pathname
  
  io = 0; cnt_atm = 0; cnt_ocn = 0
  print *, "cnt_obs_phys", cnt_obs_phys

  WRITE(pathname,'("d",I0.3,"/")') dr_num
  if (do_atm_obs) then
  endif

  if (do_ocn_obs) then
  endif

  if (do_drf_vel) then
    open(803,file=pathname//'yobs_drf_vel.dat',action="read",access="sequential")
    do i = 1, nt
      do j = cnt_atm+cnt_ocn+1, cnt_obs_phys
        read(803,*) wk
        yobs(j,1,i) = wk(2) !oelem
        yobs(j,2,i) = wk(3) !olon
        yobs(j,3,i) = wk(4) !olat
        yobs(j,4,i) = wk(5) !odat
      enddo
    enddo
    close(803)
    print*, "shape(yobs) = ", SHAPE(yobs)
    print *, "yobs(cnt_obs_phys,:,end)=", yobs(cnt_obs_phys,:,nt)
  endif
endsubroutine
endmodule
