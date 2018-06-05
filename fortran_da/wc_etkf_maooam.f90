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
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a, Xens_aa, Xens_ao
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam, Xam_a, Xam_o
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,n, cnt_obs
  INTEGER :: nt
  REAL(KIND=8) :: t, t_tmp
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=32) ::filename
  CHARACTER(LEN=35) :: gainname
  CHARACTER(LEN=35) :: gainname2
  CHARACTER(LEN=24) :: freename
  CHARACTER(LEN=33) :: sprdname
  CHARACTER(LEN=3) :: expname  
 
  CALL init_aotensor    ! Compute the tensor
  CALL init_integrator  ! Initialize the integrator

  nt = INT(t_run/dt)

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
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
  READ(100,*) X(0:ndim)
  X(0) = 1.D0
  CLOSE(100)

! load obs
  PRINT*, 'Start loading obs...'
  nt = INT(t_run/tw)+2
  ALLOCATE(yobs(ndim,nt))
  yobs = 0.0

  IF (nobs == ndim) THEN
    OPEN(102,file='yobs_atm.dat',action="read",access="sequential")
    OPEN(108,file='yobs_ocn.dat',action="read",access="sequential")
    !OPEN(108,file='nature.dat',action="read",access="sequential")
    expname="wcp"
  ELSE
    PRINT*, "ERROR: wrong nobs."    
  ENDIF
  ALLOCATE(wk(nobs+1))
  DO i=1,nt
    IF (nobs == ndim) THEN
       READ(102,*) wk(1:2*natm+1)
       yobs(1:2*natm,i) = wk(2:2*natm+1)
       READ(108,*) wk(1:2*noc+1)
       yobs(2*natm+1:ndim,i) = wk(2:2*noc+1)
       !READ(108,*) wk(1:ndim+1)
       !yobs(1:ndim,i) = wk(2:ndim+1)
     ELSE
       PRINT *, "ERROR: wrong nobs"
     ENDIF
  ENDDO
  PRINT*, "yobs(:,end) = ", yobs(:,nt)
  PRINT*, "nt = ", nt
  CLOSE(102)
  CLOSE(108)
  IF (nobs==ndim) CLOSE(107)
  PRINT*, 'Read R matrix...'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO j = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(j)
    R(j) = R(j) * 0.1 ! Try 10% climatology 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

! set initial ensembles
  ALLOCATE(err(ndim))
  DO k = 1,ens_num
    CALL randn(ndim,err)
    IF (nobs == ndim) THEN
      Xens(:,k) = X(1:ndim) + DOT_PRODUCT(err,R)*ini_err
    ELSE
      PRINT *, "ERROR: wrong nobs"
    ENDIF
    PRINT*, k, Xens(:,k)
  ENDDO

! set initial ensembles
!  IF (tw_solo .lt. 1) THEN
!     write(freename,'("freerun_atm",F0.3,".dat")') tw_solo
!  ELSE
!     write(freename,'("freerun_atm",I4.4,".dat")') INT(tw_solo)
!  ENDIF
!  OPEN(106,file=freename,action="read",access="sequential")
!  READ(106,*) X_free(1,0:ndim)
!  X_free(1,0) = 1.D0
!  print *, X_free(1,0:ndim)
!  CLOSE(106)
!  t=0.D0
!  DO n = 1, nt-1
!     call step(X_free(n,:),t,dt,X_free(n+1,:))
!  enddo
!  print *, "Start assigning ensemble values......"
!  do n = 1, ens_num
!     k = (nt/ens_num)*n
!     print *, "k = ", k
!     Xens(1:2*natm,n) = X_free(k,1:2*natm)
!     print*, "ensemble: mem, step=", n, k, nt
!  enddo

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
      
      ! DA on atm components
      CALL etkf( 2*natm, 2*natm, ens_num, Xens(1:2*natm,:),&
                 luse(1:2*natm), yobs(1:2*natm,cnt_obs+1), R(1:2*natm), infl,&
                 Xens_aa, Xam_a, 105)

      ! DA on ocn components
      CALL etkf( 2*noc, 2*noc, ens_num, Xens(2*natm+1:ndim,:), &
                 luse(2*natm+1:ndim), yobs(2*natm+1:ndim,cnt_obs+1), &
                 R(2*natm+1:ndim), infl, Xens_ao, Xam_o, 107)

      Xens_a(1:2*natm,:) = Xens_aa
      Xam(1:2*natm) = Xam_a
      Xens_a(2*natm+1:ndim,:) = Xens_ao
      Xam(2*natm+1:ndim) = Xam_o

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
    CLOSE(105)
    CLOSE(107)
    CLOSE(106)
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
