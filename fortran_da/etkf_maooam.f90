PROGRAM etkf_maooam
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, nobs, infl, ini_err
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf
  USE m_mt,     only: randn
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: X_free
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,n, cnt_obs
  INTEGER :: nt, nseed
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL(KIND=8) :: t, t_tmp
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=32) ::filename
  CHARACTER(LEN=33) :: gainname 
  CHARACTER(LEN=24) :: freename
  CHARACTER(LEN=33) :: sprdname
  CHARACTER(LEN=33) :: ensaname
  CHARACTER(LEN=3) :: expname  
 
  CALL init_aotensor    ! Compute the tensor
  CALL init_integrator  ! Initialize the integrator

  nt = INT(t_run/dt)

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(Xens(ndim,ens_num),Xens_a(ndim,ens_num))
  ALLOCATE(R(ndim),luse(ndim)); luse = .true.
  ALLOCATE(Xbm(ndim),Xam(ndim))
  ALLOCATE(std_ens(ndim))

!----------------------------------------------------------------

  PRINT*, 'Model MAOOAM v1.3 for ETKF'
  PRINT*, 'Loading information...'

  OPEN(100,file='nature.dat',action="read",access="sequential")
  READ(100,*) X(0:ndim)
  X(0) = 1.D0
  CLOSE(100)

! load obs
  PRINT*, 'Start loading obs...'
  nt = INT(t_run/tw)+2
  ALLOCATE(yobs(ndim,nt))
  yobs = 0.0
  If (nobs == 2*natm) THEN
    OPEN(102,file='yobs_atm.dat',action="read",access="sequential")
    expname="atm"  
  ELSEIF (nobs ==  2*noc) THEN
    OPEN(102,file='yobs_ocn.dat',action="read",access="sequential")
    expname="ocn"
  ELSEIF (nobs == ndim) THEN
    OPEN(102,file='yobs_atm.dat',action="read",access="sequential")
    OPEN(107,file='yobs_ocn.dat',action="read",access="sequential")
    expname="cpl"    
  ENDIF
  ALLOCATE(wk(nobs+1))
  DO i=1,nt
     IF (nobs==2*natm) THEN
       READ(102,*) wk
       yobs(1:2*natm,i) = wk(2:nobs+1)
     ELSEIF (nobs==2*noc) THEN
       READ(102,*) wk
       yobs(2*natm+1:ndim,i) = wk(2:nobs+1)
     ELSEIF (nobs == ndim) THEN
       READ(102,*) wk(1:2*natm+1)
       yobs(1:2*natm,i) = wk(2:2*natm+1)
       READ(107,*) wk(1:2*noc+1)
       yobs(2*natm+1:ndim,i) = wk(2:2*noc+1)
     ENDIF
  ENDDO
  PRINT*, "yobs(:,end) = ", yobs(:,nt)
  PRINT*, "nt = ", nt
  CLOSE(102)
  IF (nobs==ndim) CLOSE(107)
  PRINT*, 'Read R matrix...'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO j = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(j)
    R(j) = R(j) ! Try 10% climatology 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

! set up the luse for obs space
  If (nobs == 2*natm) THEN
    luse(1:2*natm) = .true.
    luse(2*natm+1:ndim) = .false.  
  ELSEIF (nobs ==  2*noc) THEN
    luse(1:2*natm) = .false.
    luse(2*natm+1:ndim) = .true.    
  ENDIF

! set initial ensembles
  ALLOCATE(err(ndim))
  !CALL RANDOM_SEED(size=nseed)
  !ALLOCATE(seed(nseed))
  !seed(1) = 3333 
  !CALL RANDOM_SEED(put=seed) 
  DO k = 1,ens_num
    CALL randn(ndim,err)
    IF (nobs == 2*natm) THEN
      !Xens(:,k) = X(1:ndim) + DOT_PRODUCT(err,R)*ini_err*0.1
      Xens(:,k) = X(1:ndim) + (err*R)*ini_err
      print*, k, (err*R)*ini_err
      !Xens(:,k) = X(1:ndim) + (err)*0.0001*ini_err
    ELSEIF (nobs == 2*noc) THEN
      Xens(:,k) = X(1:ndim) + (err*R)
    ELSEIF (nobs == ndim) THEN
      Xens(:,k) = X(1:ndim) + (err*R)*ini_err
    ENDIF
    !PRINT*, k, Xens(:,k)
    PRINT*, k, err
  ENDDO

  t=0.D0
  t_up=dt/t_run_da*100.D0
  cnt_obs=1

  IF (writeout) THEN
    WRITE(filename,'("Xam_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(104,file=filename)
    WRITE(gainname,'("gain_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname, INT(ens_num), infl,tw_da, INT(ini_err)
    !OPEN(105,file=gainname)
    WRITE(sprdname,'("sprd_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname, INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(106,file=sprdname)
  END IF

  ! run the DA cycle
  DO WHILE (t<t_run_da)
  
    ! generate forecast
    DO k = 1, ens_num
      t_tmp = 0.d0
      X(1:ndim) = Xens(:,k)
      CALL step(X,t_tmp,dt,Xnew)
      Xens(:,k) = Xnew(1:ndim)
    ENDDO

    t = t + dt

    ! generate analysis
    IF (mod(t,tw_da)<dt) THEN

      CALL etkf( ndim, ndim, ens_num, Xens, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)

      Xens = Xens_a
    ELSE
      DO i = 1, ndim
        Xam(i) = SUM(Xens(i,:))/ens_num
      END DO
    END IF

    IF (mod(t,tw)<dt) THEN
    ! Write ensemble spread
      DO i = 1, ndim
        std_ens(i) = 0.d0
        DO j = 1, ens_num
          std_ens(i) = std_ens(i) + (Xens(i,j)-Xam(i))*(Xens(i,j)-Xam(i))
        END DO
        std_ens(i) = SQRT(std_ens(i)/(ens_num-1))
      END DO
      IF (writeout) WRITE(106,*) t,std_ens
      
      ! Write Xam
      IF (writeout) WRITE(104,*) t, Xam

      ! Updat cnt_obs
      cnt_obs = cnt_obs + 1
    END IF

    IF (mod(t/t_run_da*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)

  ENDDO

  IF (writeout) THEN
    CLOSE(104)
    !CLOSE(105)
    !CLOSE(106)
  END IF

  IF (.false.) THEN
    WRITE(ensaname,'("Xens_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,t_run_da, INT(ini_err)
    OPEN(200,file=ensaname)
    DO j = 1, ens_num
      WRITE(200,*) j, Xens(:,j)
    END DO
    CLOSE(200)
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



END PROGRAM etkf_maooam
