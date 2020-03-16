PROGRAM etkf_maooam_phys
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, infl, ini_err, spin_up
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  !USE m_da_maooam, only: etkf
  USE m_mt,     only: randn

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size, oms, n
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr

  USE params, only: do_atm_obs, do_ocn_obs, do_drf_pos, do_drf_vel
  USE m_da_maooam_phys, only : read_dimension_phys, read_dimension_drf
  USE m_da_maooam_phys, only : read_obs_phys
  USE m_da_maooam_phys, only : etkf_phys 
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_dr,X_dummy
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew_dr
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a
  REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R, R_phys
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,nn,m, cnt_obs
  INTEGER :: cnt_obs_phys, cnt_obs_drf
  INTEGER :: Ho, Po, total
  INTEGER :: nt, nseed, nt_spinup
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL(KIND=8) :: t, t_tmp
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk, wk_dr
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: err_dr
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: std_ens
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=37) ::filename
  CHARACTER(LEN=38) :: gainname 
  CHARACTER(LEN=24) :: freename
  CHARACTER(LEN=38) :: sprdname
  CHARACTER(LEN=38) :: ensaname
  CHARACTER(LEN=3) :: expname  
  CHARACTER(LEN=5) :: pathname
 
  CALL init_aotensor    ! Compute the tensor
  IF (do_drf_pos) THEN
    CALL init_integrator_dr
  END IF
  CALL init_integrator  ! Initialize the integrator

  IF (do_drf_pos) THEN
    total = ndim+ndim_dr
  ELSE
    total = ndim
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim)); ALLOCATE(X_dr(ndim_dr),Xnew_dr(ndim_dr))
  ALLOCATE(X_dummy(0:total))
  ALLOCATE(Xens(total,ens_num),Xens_a(total,ens_num)); Xens=0.d0
  ALLOCATE(R(total))
  ALLOCATE(Xbm(total),Xam(total))
  ALLOCATE(std_ens(total))

  nt = INT(t_run/tw)
!----------------------------------------------------------------

  PRINT*, 'Model MAOOAM v1.3 for ETKF'
  PRINT*, 'Loading information...'
  
  t=0.d0
  WRITE(pathname,'("d",I0.3,"/")') dr_num
  OPEN(100,file=pathname//'nature.dat',action="read",access="sequential")
  DO WHILE (t - spin_up <dt)
     READ(100,*) X_dummy(0:total)
     IF (t==0.d0) THEN
       X(0)=1.D0 
       X(1:ndim)=X_dummy(1:ndim) ! Assign true fluid IC
     ENDIF
     IF (do_drf_pos .AND. ABS(t-spin_up)<dt) THEN
       X_dr=X_dummy(ndim+1:ndim+ndim_dr) ! Assign true drifter IC at t=spin_up
       PRINT *, "DEBUG: t=",t,"X_dr=",X_dr
     ENDIF
     t=t+tw
  ENDDO
  CLOSE(100)
  DEALLOCATE(X_dummy)

  ! load obs
  PRINT*, 'Start loading obs...'
  ! step 1: read dimension
  cnt_obs_phys = 0; cnt_obs_drf = 0;
  CALL read_dimension_phys(do_atm_obs, do_ocn_obs, do_drf_vel, cnt_obs_phys)
  CALL read_dimension_drf(do_drf_pos, cnt_obs_drf)
  ALLOCATE(yobs(cnt_obs_phys+cnt_obs_drf,4,nt),R_phys(cnt_obs_phys+cnt_obs_drf))  
 
  ! step 2: read observation 
  CALL read_obs_phys(do_atm_obs, do_ocn_obs, do_drf_vel, cnt_obs_phys, nt, yobs(1:cnt_obs_phys,:,:))
  !CALL read_obs_drf(do_drf_pos, cnt_obs_drf, nt, yobs)

  PRINT*, 'Read R matrix...'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  R = 0.d0
  DO j = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(j)
    R(j) = R(j) ! Try 10% climatology 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

! set initial ensembles
  ALLOCATE(err(ndim),err_dr(ndim))
  !CALL RANDOM_SEED(size=nseed)
  !ALLOCATE(seed(nseed))
  !seed(1) = 3333 
  !CALL RANDOM_SEED(put=seed) 
  DO k = 1,ens_num
    CALL randn(ndim,err)
    Xens(1:ndim,k) = X(1:ndim) + (err*R)*ini_err
    !PRINT*, k, err

    IF (do_drf_pos .AND. spin_up==0.d0) THEN
      Xens(ndim+1:ndim+ndim_dr,k) = X_dr(1:ndim_dr) !+ (err_dr*R(ndim+1:ndim+ndim_dr))
    END IF   
  ENDDO
  DEALLOCATE(err)

! set initial ensemble for drifters

  t=0.D0
  t_up=dt/t_run_da*100.D0
  cnt_obs=1

  expname='vel'

  IF (writeout) THEN
    WRITE(filename,'("Xam_etkf_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') &
& INT(dr_num),expname,INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(104,file=filename)
    WRITE(gainname,'("gain_etkf_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') & 
& INT(dr_num),expname, INT(ens_num), infl,tw_da, INT(ini_err)
    !OPEN(105,file=gainname) !output the gain matrix for computing Lyapunov Exponents
    WRITE(sprdname,'("sprd_etkf_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') &
& INT(dr_num),expname, INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(106,file=sprdname)
  END IF

  ! run the DA cycle
  PRINT*, "Start Running DA"
  DO WHILE (t<t_run_da)
  
    ! generate forecast
    DO k = 1, ens_num
      t_tmp = t 
      IF (do_drf_pos .AND. t-spin_up>-dt ) THEN
        X(1:ndim) = Xens(1:ndim,k); X_dr = Xens(ndim+1:ndim+ndim_dr,k)
        CALL step_dr(X,X_dr,t_tmp,dt,Xnew,Xnew_dr)
        Xens(1:ndim,k) = Xnew(1:ndim); Xens(ndim+1:ndim+ndim_dr,k) = Xnew_dr
      ELSE
        X(1:ndim) = Xens(1:ndim,k)
        CALL step(X,t_tmp,dt,Xnew)
        Xens(1:ndim,k) = Xnew(1:ndim)
      END IF
    ENDDO

    t = t + dt

    ! deploy drifters within fluid ensemble members
    IF (do_drf_pos .AND. ABS(t-spin_up)<dt) THEN
      DO k = 1, ens_num
        Xens(ndim+1:ndim+ndim_dr,k) = X_dr(1:ndim_dr)
      ENDDO
    ENDIF

    ! generate analysis
    IF (mod(t,tw_da)<dt .AND. t-spin_up>dt) THEN
      CALL etkf_phys( total, cnt_obs_phys, ens_num, 4, Xens, yobs(:,:, cnt_obs), infl, Xens_a, Xam, 105)

      Xens = Xens_a
    ELSE
      DO i = 1, total
        Xam(i) = SUM(Xens(i,:))/ens_num
      END DO
    END IF

    IF (mod(t,tw)<dt) THEN
    ! Write ensemble spread
      DO i = 1, total
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
    CLOSE(106)
  END IF

  !IF (.false.) THEN
  !  WRITE(ensaname,'("Xens_etkf_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') expname,INT(ens_num), infl,t_run_da, INT(ini_err)
  !  OPEN(200,file=ensaname)
  !  DO j = 1, ens_num
  !    WRITE(200,*) j, Xens(:,j)
  !  END DO
  !  CLOSE(200)
  !END IF
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



END PROGRAM etkf_maooam_phys
