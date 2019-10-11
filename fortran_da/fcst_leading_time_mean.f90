PROGRAM fcst_leading_time
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, infl, ini_err, t_lead
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf
  USE m_mt,     only: randn
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,n
  INTEGER :: nt, nseed
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL(KIND=8) :: t, t_tmp
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=32) ::filename
  CHARACTER(LEN=33) :: ensaname
  CHARACTER(LEN=3) :: expname  
 
  CALL init_aotensor    ! Compute the tensor
  CALL init_integrator  ! Initialize the integrator

  nt = INT(t_run/dt)

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(Xens(ndim,ens_num))
  ALLOCATE(R(ndim),luse(ndim)); luse = .true.
  ALLOCATE(Xbm(ndim),Xam(ndim))
  ALLOCATE(std_ens(ndim))

!----------------------------------------------------------------

  PRINT*, 'Model MAOOAM v1.3 for ETKF'
  PRINT*, 'Loading information...'

! set initial ensembles
   expname='wcp'
   ALLOCATE(wk(ndim+1))
   WRITE(ensaname,'("Xam_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,t_run_da, INT(ini_err)
   OPEN(200,file=ensaname)
   DO i = 1, ens_num
     READ(200,*) wk
   END DO  
   Xam = wk(2:ndim+1)
   CLOSE(200)

  t=0.D0
  t_up=dt/t_run_da*100.D0

  IF (writeout) THEN
    WRITE(filename,'("Xam_lead_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(104,file=filename)
  END IF

  ! run the DA cycle
  DO WHILE (t<t_lead)
  
    ! generate forecast
    t_tmp = 0.d0
    X(1:ndim) = Xam
    X(0) = 1.d0
    CALL step(X,t_tmp,dt,Xnew)
    Xam = Xnew(1:ndim)

    t = t + dt

    IF (mod(t,tw)<dt) THEN      
      ! Write Xam
      IF (writeout) WRITE(104,*) t, Xam
    END IF

    IF (mod(t/t_run_da*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)

  ENDDO

  IF (writeout) THEN
    CLOSE(104)
  END IF
!  x(0,2:nt) = 0.0d0
!  do n = 1, nt-1
!     do k = 1, nens
!        ! ensemble forecast
!        call step(Xens(:,k), t, dt, Xens(:,k))
!        ! H(x)
!        Yens(:,k) = Xens(1:ndim,k)
!     enddo

     ! luse
     ! (1:natm_maooam): atm streamfunction
     ! (natm_maooam+1:2*natm_maooam): atm temp
     ! (2*natm_maooam:2*natm_maooam+nocn_maooam): ocn streamfunction
     ! (2*natm_maooam+1:end): ocn temp
     !
!     call etkf( ndim, &
!                ndim, &
!                nens, &
!                Xens(1:ndim,:), &
!                luse, &
!                yobs(:,n+1), &
!                R, &
!                infl, &
!                Xens(1:ndim,:), &
!                x(1:ndim,n+1))

!  enddo

!  do n = 1, nt
!     write(lout_etkf,"(10000(D24.17,1x))") (x(k,n),k=1,ndim)
!  enddo

  PRINT*, 'Evolution finished.'

END PROGRAM fcst_leading_time
