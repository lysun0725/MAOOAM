PROGRAM etkf_maooam
  USE params, only: ndim, dt, tw, t_run,writeout, tw_solo
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, nobs, infl
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf
  USE m_mt,     only: randn
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_solo
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_solo_new
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: X_free
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Yens
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: X_ocn
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,n, cnt_obs
  INTEGER :: nt
  REAL(KIND=8) :: t, t_tmp
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=25) ::filename
  CHARACTER(LEN=26) :: gainname 
  CHARACTER(LEN=24) :: freename
  CHARACTER(LEN=26) :: sprdname
  
 
  CALL init_aotensor    ! Compute the tensor
  CALL init_integrator  ! Initialize the integrator

  nt = INT(t_run/dt)

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(X_solo(0:ndim),X_solo_new(0:ndim))
  ALLOCATE(Xens(2*natm,ens_num),Yens(nobs,ens_num),Xens_a(2*natm,ens_num))
  ALLOCATE(X_ocn(2*noc))
  ALLOCATE(R(nobs))
  ALLOCATE(luse(nobs)); luse=.true.
  ALLOCATE(Xbm(2*natm),Xam(2*natm))
  ALLOCATE(X_free(nt,ndim))
  ALLOCATE(std_ens(2*natm))

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
  ALLOCATE(yobs(nobs,nt))
  OPEN(102,file='yobs.dat',action="read",access="sequential")
  ALLOCATE(wk(nobs+1))
  DO i=1,nt
     READ(102,*) wk
     yobs(:,i) = wk(2:nobs+1)
  ENDDO
  PRINT*, "yobs(:,end) = ", yobs(:,nt)
  PRINT*, "nt = ", nt
  CLOSE(102)
  PRINT*, 'Read R matrix...'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO j = 1, ndim
    IF (j .gt. 2*natm) CYCLE
    READ(103,"(10000(D24.17,1x))") R(j)
    R(j) = R(j) * 0.1 ! Try 10% climatology 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

! set initial ensembles
  ALLOCATE(err(2*natm))
  DO k = 1,ens_num
    CALL randn(2*natm,err)
    Xens(:,k) = X(1:2*natm) + DOT_PRODUCT(err,R)
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
  X_ocn = X(2*natm+1:ndim)
  cnt_obs=1

  IF (writeout) THEN
    WRITE(filename,'("Xam_etkf_",I2.2,"_",F3.1,E6.1,".dat")') INT(ens_num), infl,tw_da
    OPEN(104,file=filename)
    WRITE(gainname,'("gain_etkf_",I2.2,"_",F3.1,E6.1,".dat")') INT(ens_num), infl,tw_da
    OPEN(105,file=gainname)
    WRITE(sprdname,'("sprd_etkf_",I2.2,"_",F3.1,E6.1,".dat")') INT(ens_num), infl,tw_da
    OPEN(106,file=sprdname)
  END IF

  ! run the DA cycle
  DO WHILE (t<t_run_da)
  
    ! generate forecast
    DO k = 1, ens_num
      t_tmp = 0.d0
      X_solo = (/1.D0, Xens(:,k), X_ocn/)
      CALL step(X_solo,t_tmp,dt,X_solo_new)
      Xens(:,k) = X_solo_new(1:2*natm)
      Yens(:,k) = Xens(:,k)
    ENDDO

    t_tmp=0.d0
    ! generate the driving signal
    CALL step(X,t_tmp,dt,Xnew)
    X=Xnew
    t = t + dt

    ! generate analysis
    IF (mod(t,tw_da)<dt) THEN

      CALL etkf( 2*natm, nobs, ens_num, Xens, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)

      Xens = Xens_a
    ELSE
      DO i = 1, 2*natm
        Xam(i) = SUM(Xens(i,:))/ens_num
      END DO
    END IF

    IF (mod(t,tw)<dt) THEN
    ! Write ensemble spread
      DO i = 1, 2*natm
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

    ! update the driving signal
    IF (mod(t,tw_solo)<dt) THEN
      X_ocn = X(2*natm+1:ndim)
    END IF

    IF (mod(t/t_run_da*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)

  ENDDO

  IF (writeout) THEN
    CLOSE(104)
    CLOSE(105)
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
