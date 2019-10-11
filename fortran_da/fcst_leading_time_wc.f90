PROGRAM etkf_maooam
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, nobs, infl, ini_err,t_lead
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf
  USE m_mt,     only: randn
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xhist
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_lead_true
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: X_lead, X_lead_new
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: err_lead_ave, err_lead_sum
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: X_free
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a, Xens_aa, Xens_ao
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam, Xam_a, Xam_o
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,n, cnt_obs, cnt_lead, cnt_da=0
  INTEGER :: nt, nseed, nt_lead
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL(KIND=8) :: t, t_tmp, t_tmp_lead=0.D0, t_tmp_lead2=0.D0
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=32) ::filename
  CHARACTER(LEN=35) :: gainname
  CHARACTER(LEN=35) :: gainname2 
  CHARACTER(LEN=24) :: freename
  CHARACTER(LEN=33) :: sprdname
  CHARACTER(LEN=33) :: ensaname
  CHARACTER(LEN=3) :: expname  
 
  CALL init_aotensor    ! Compute the tensor
  CALL init_integrator  ! Initialize the integrator

  nt = INT(t_run/tw)+2
  nt_lead = INT(t_lead/tw)+2

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(Xhist(nt,0:ndim))
  ALLOCATE(X_lead_true(0:ndim))
  ALLOCATE(X_lead(0:ndim,ens_num),X_lead_new(0:ndim,ens_num))
  ALLOCATE(err_lead_ave(nt_lead,ndim),err_lead_sum(nt_lead,ndim)); err_lead_sum=0.d0
  ALLOCATE(Xens(ndim,ens_num),Xens_a(ndim,ens_num))
  ALLOCATE(Xens_aa(2*natm,ens_num),Xens_ao(2*noc,ens_num))
  ALLOCATE(Xam_a(2*natm),Xam_o(2*noc))
  ALLOCATE(R(ndim),luse(ndim)); luse = .true.
  ALLOCATE(Xbm(ndim),Xam(ndim))
  ALLOCATE(std_ens(ndim))

!----------------------------------------------------------------

  PRINT*, 'Model MAOOAM v1.3 for ETKF'
  PRINT*, 'Loading information...'

  OPEN(100,file='nature.dat',action="read",access="sequential")
  DO j=1,nt  
    READ(100,*) Xhist(j,0:ndim)
    Xhist(j,0) = 1.D0
  ENDDO
  X=Xhist(1,:)
  CLOSE(100)

! load obs
  PRINT*, 'Start loading obs...'
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
    expname="wcp"    
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
  CALL RANDOM_SEED(size=nseed)
  ALLOCATE(seed(nseed))
  seed(1) = 4444
  CALL RANDOM_SEED(put=seed) 
  DO k = 1,ens_num
    CALL randn(ndim,err)
    IF (nobs == 2*natm) THEN
      Xens(:,k) = X(1:ndim) + DOT_PRODUCT(err,R)*ini_err*0.1
    ELSEIF (nobs == 2*noc) THEN
      Xens(:,k) = X(1:ndim) + DOT_PRODUCT(err,R)
    ELSEIF (nobs == ndim) THEN
      Xens(:,k) = X(1:ndim) + DOT_PRODUCT(err,R)*ini_err*0.1
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
    WRITE(gainname,'("gain_etkf_a_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname, INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(105,file=gainname)
    WRITE(gainname2,'("gain_etkf_o_",A3,"_",I2.2,"_",F3.1,E6.1"_",I2.2,".dat")') expname, INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(107,file=gainname2)
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

      CALL etkf( 2*natm, 2*natm, ens_num, Xens(1:2*natm,:), &
                 luse(1:2*natm), yobs(1:2*natm,cnt_obs+1), R(1:2*natm), infl,&
                 Xens_aa, Xam_a, 105)
      
      CALL etkf( 2*noc, 2*noc, ens_num, Xens(2*natm+1:ndim,:), &
                 luse(2*natm+1:ndim), yobs(2*natm+1:ndim,cnt_obs+1), &
                 R(2*natm+1:ndim), infl, Xens_ao, Xam_o, 107)

      Xens_a(1:2*natm,:) = Xens_aa
      Xam(1:2*natm) = Xam_a
      Xens_a(2*natm+1:ndim,:) = Xens_ao
      Xam(2*natm+1:ndim) = Xam_o
     
      Xens = Xens_a
    
      ! Start doing fcst leading time
      t_tmp_lead = 0.D0
      cnt_lead = 1
      cnt_da = cnt_da + 1
      X_lead(1:ndim,:) = Xens_a; X_lead(0,:)=1.D0; X_lead_new(1:ndim,:) = 0.D0; X_lead_new(0,:)=1.D0
      X_lead_true(0:ndim) = Xhist(cnt_obs+1,0:ndim)

      DO k = 1, ens_num
        err_lead_sum(1,:) = err_lead_sum(1,:) + ABS(X_lead(1:ndim,k)-X_lead_true(1:ndim))
      ENDDO
      err_lead_ave(1,:) = err_lead_sum(1,:)/cnt_da/ens_num
      
      DO WHILE (t_tmp_lead<t_lead)
        DO k = 1, ens_num
          t_tmp_lead2 = t_tmp_lead
          CALL step(X_lead(:,k),t_tmp_lead2,dt,X_lead_new(:,k))
          X_lead(:,k)=X_lead_new(:,k)
        ENDDO
        t_tmp_lead = t_tmp_lead2
        IF (mod(t_tmp_lead,tw)<dt) THEN
          cnt_lead = cnt_lead + 1
          X_lead_true(0:ndim) = Xhist(cnt_obs+cnt_lead,0:ndim)
          DO k = 1,ens_num
            err_lead_sum(cnt_lead,:) = err_lead_sum(cnt_lead,:) + ABS(X_lead(1:ndim,k)-X_lead_true(1:ndim))
          ENDDO
          err_lead_ave(cnt_lead,:) = err_lead_sum(cnt_lead,:)/cnt_da/ens_num
        END IF
      ENDDO

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
    CLOSE(105)
    CLOSE(106)
  END IF

  IF (writeout) THEN
    WRITE(ensaname,'("Xerr_lead_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,t_run_da, INT(ini_err)
    OPEN(200,file=ensaname)
    DO j = 1, cnt_lead
      WRITE(200,*) j, err_lead_ave(j,:)
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


  PRINT*, 'cnt_da = ',cnt_da
  PRINT*, 'Evolution finished.'



END PROGRAM etkf_maooam
