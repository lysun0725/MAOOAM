PROGRAM hybrid_maooam
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, nobs, infl, ini_err, spin_up
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf, etkf_w, enkf_w,pf, pf2, MH_resampling, RSR_resampling_1D
  USE m_mt,     only: randn

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size, oms, n, nobs_dr
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr
  
  USE params, only: part_num, do_hybrid
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_dr,X_dummy
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xnew
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xnew_dr
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens ! Driving ensemble (1:2*natm,ens_num)
  REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Xdr_part 
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a
  REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Xdr_part_a 
  REAL(KIND=8),DIMENSION(:,:), ALLOCATABLE :: W, W_a
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: W_tilde, W_tilde_a, W_dummy
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: yobs
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: R
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xbm
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Xam
  LOGICAL, DIMENSION(:),ALLOCATABLE :: luse

  INTEGER :: i,j,k,l,nn,m, cnt_obs
  INTEGER :: Ho, Po, total
  INTEGER :: nt, nseed, nt_spinup,nsize,time
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
  REAL(KIND=8) :: N_eff  

  REAL(KIND=8) :: spin_up_etkf
  spin_up_etkf = 200.0 ! 20 days 

  CALL init_aotensor    ! Compute the tensor
  IF (do_drifter) THEN
    CALL init_integrator_dr
  END IF
  CALL init_integrator  ! Initialize the integrator

  IF (do_drifter) THEN
    total = ndim+ndim_dr
  ELSE
    total = ndim
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim)); ALLOCATE(X_dr(ndim_dr),Xnew_dr(ndim_dr,part_num))
  ALLOCATE(X_dummy(0:total))
  ALLOCATE(Xens(total,ens_num),Xens_a(total,ens_num)); Xens=0.d0
  ALLOCATE(Xdr_part(ndim_dr,ens_num,part_num)); Xdr_part=0.d0
  ALLOCATE(Xdr_part_a(ndim_dr,ens_num,part_num)); Xdr_part_a=0.d0
  ALLOCATE(R(total),luse(total)); luse = .true.
  ALLOCATE(Xbm(total),Xam(total))
  ALLOCATE(std_ens(total))
  ALLOCATE(W_tilde(ens_num),W_tilde_a(ens_num),W_dummy(ens_num))
  ALLOCATE(W(ens_num,part_num),W_a(ens_num,part_num))
  W = 1.d0/ens_num/part_num; W_a = 1.d0/ens_num/part_num
  W_tilde = 1.d0/ens_num; W_tilde_a = 1.d0/ens_num

  nt = INT(t_run/tw)+1
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
     IF (do_drifter .AND. ABS(t-spin_up)<dt) THEN
       X_dr=X_dummy(ndim+1:ndim+ndim_dr) ! Assign true drifter IC at t=spin_up
       PRINT *, "DEBUG: t=",t,"X_dr=",X_dr
     ENDIF
     t=t+tw
  ENDDO
  CLOSE(100)
  DEALLOCATE(X_dummy)

! load obs
  PRINT*, 'Start loading obs...'
  ALLOCATE(yobs(total,nt))
  yobs = 0.0
  IF (nobs == 2*natm) THEN
    OPEN(102,file=pathname//'yobs_atm.dat',action="read",access="sequential")
    expname="atm"  
  ELSEIF (nobs ==  2*noc) THEN
    OPEN(102,file=pathname//'yobs_ocn.dat',action="read",access="sequential")
    expname="ocn"
  ELSEIF (nobs == ndim) THEN
    OPEN(102,file=pathname//'yobs_atm.dat',action="read",access="sequential")
    OPEN(107,file=pathname//'yobs_ocn.dat',action="read",access="sequential")
    expname="cpl"    
  ELSEIF (nobs == 0 .AND. do_drifter) THEN
    expname='drf'
  ELSE
    PRINT*, "ERROR: WRONG NOBS INPUT IN DA_PARAMS.NML"
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
  IF (nobs > 0) CLOSE(102)
  IF (nobs==ndim) CLOSE(107)
  DEALLOCATE(wk)

  !load drifter obs
  PRINT*, 'Start loading drifter obs...'
  IF (do_drifter) THEN
    OPEN(108,file=pathname//'yobs_drf.dat',action="read",access="sequential")
    ALLOCATE(wk_dr(nobs_dr+1)) 
    DO i=1,nt
      READ(108,*) wk_dr
      yobs(ndim+1:ndim+ndim_dr,i) = wk_dr(2:nobs_dr+1)
    ENDDO
    PRINT*, "yobs(dr,end) = ", yobs(ndim+1:ndim+ndim_dr,nt)
    CLOSE(108)
    DEALLOCATE(wk_dr)
  END IF

  PRINT*, 'Read R matrix...'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  R = 0.d0
  DO j = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(j)
    R(j) = R(j) ! Try 10% climatology 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

  IF (do_drifter) THEN
      ! Define R_dr
      DO m = 1,noc
        Ho=oms(m,1); Po=oms(m,2)
        R(ndim+1) = R(ndim+1) + (R(m+2*natm)**2)*(Po**2)/noc
        R(ndim+2) = R(ndim+2) + (R(m+2*natm)**2)*((Ho*n/2)**2)/noc
      ENDDO

      R(ndim+1)=SQRT(R(ndim+1))*10; R(ndim+2)=SQRT(R(ndim+2))*10
      !R(ndim+1)=1.D-2; R(ndim+2)=1.D-3
      PRINT*, "R(ndim+1) = ", R(ndim+1), "R(ndim+2) = ", R(ndim+2)

      DO nn = 2,dr_num
        R(ndim+(nn-1)*dr_size+1) = R(ndim+1)
        R(ndim+(nn-1)*dr_size+2) = R(ndim+2)
      ENDDO
  END IF

! set up the luse for obs space
  IF (nobs == 2*natm) THEN
    luse(1:2*natm) = .true.
    luse(2*natm+1:ndim) = .false.  
  ELSEIF (nobs ==  2*noc) THEN
    luse(1:2*natm) = .false.
    luse(2*natm+1:ndim) = .true.    
  ELSEIF (nobs == 0 .AND. do_drifter) THEN
    luse(1:ndim) = .false.
    luse(ndim+1:ndim+ndim_dr) = .true.
  ENDIF

! set initial ensembles
  ALLOCATE(err(ndim),err_dr(ndim_dr))
  call system_clock(time)
  call random_seed(size=nsize)
  allocate(seed(nsize), source=time+37*[(i,i=0,nsize-1)])
  call random_seed(put=seed)

  !CALL RANDOM_SEED(size=nseed)
  !ALLOCATE(seed(nseed))
  !seed(1) = 3333 
  !CALL RANDOM_SEED(put=seed) 
  DO k = 1,ens_num
    CALL randn(ndim,err)
    Xens(1:ndim,k) = X(1:ndim) + (err*R)*ini_err
    !PRINT*, k, err

    !IF (do_drifter .AND. spin_up==0.d0) THEN
    IF (do_drifter) THEN
      Xens(ndim+1:ndim+ndim_dr,k) = X_dr(1:ndim_dr) !+ (err_dr*R(ndim+1:ndim+ndim_dr))
      ! initialize particles
      !$OMP PARALLEL DO
      DO l = 1, part_num
        CALL randn(ndim_dr,err_dr)
        Xdr_part(:,k,l) = X_dr(1:ndim_dr) + err_dr*R(ndim+1:ndim+ndim_dr)*0.001
      ENDDO   
      !$OMP END PARALLEL DO
    END IF   
  ENDDO

! set initial ensemble for drifters

  t=0.D0
  t_up=dt/t_run_da*100.D0
  cnt_obs=1

  IF (writeout) THEN
    WRITE(filename,'("Xam_hybr_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') &
& INT(dr_num),expname,INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(104,file=filename)
    WRITE(gainname,'("gain_hybr_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') & 
& INT(dr_num),expname, INT(ens_num), infl,tw_da, INT(ini_err)
    !OPEN(105,file=gainname) !output the gain matrix for computing Lyapunov Exponents
    WRITE(sprdname,'("sprd_hybr_d",I0.3,"_",A3,"_",I2.2,"_",F3.1,E6.1,"_",I2.2,".dat")') &
& INT(dr_num),expname, INT(ens_num), infl,tw_da, INT(ini_err)
    OPEN(106,file=sprdname)
  END IF

  ! run the DA cycle
  PRINT*, "Start Running DA"
  DO WHILE (t<t_run_da)
  
    ! generate forecast
    DO k = 1, ens_num
      t_tmp = t 
      IF (do_drifter .AND. t-spin_up>-dt ) THEN
        X(1:ndim) = Xens(1:ndim,k)
        CALL step_dr(X,Xdr_part(:,k,:),t_tmp,dt,Xnew,Xnew_dr)
        Xens(1:ndim,k) = Xnew(1:ndim); Xdr_part(:,k,:)=Xnew_dr
        W_tilde(k) = SUM(W(k,:))
        IF (W_tilde(k) < 1.d-150) THEN
          !$OMP PARALLEL DO
          DO m = 1, ndim_dr
            Xens(ndim+m,k) = SUM(Xdr_part(m,k,:))/part_num
          ENDDO
          !$OMP END PARALLEL DO
        ELSE
          Xens(ndim+1:ndim+ndim_dr,k) = 1.d0/(W_tilde(k))*MATMUL(Xdr_part(:,k,:),W(k,:))
        END IF

      ELSE
        X(1:ndim) = Xens(1:ndim,k)
        CALL step(X,t_tmp,dt,Xnew)
        Xens(1:ndim,k) = Xnew(1:ndim)
      END IF
    ENDDO

    t = t + dt

    ! deploy drifters within fluid ensemble members
    IF (do_drifter .AND. ABS(t-spin_up)<dt) THEN
      DO k = 1, ens_num
        !$OMP PARALLEL DO
        DO l = 1, part_num
          CALL randn(ndim_dr,err_dr)
          Xdr_part(:,k,l) = X_dr(1:ndim_dr) + err_dr*R(ndim+1:ndim+ndim_dr)*0.001
        ENDDO
        !$OMP END PARALLEL DO
      ENDDO
    ENDIF

    IF (mod(t,tw_da)<dt .AND. t-spin_up_etkf<=dt ) THEN
      CALL etkf( total, total, ens_num, Xens, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)
      Xens = Xens_a
      DO k = 1, ens_num
        !$OMP PARALLEL DO
        DO j = 1, part_num
          CALL randn(ndim_dr,err_dr)
          Xdr_part(:,k,j) = Xens(ndim+1:ndim+ndim_dr,k) + err_dr*R(ndim+1:ndim+ndim_dr)*0.0001
        ENDDO
        !$OMP END PARALLEL DO
      ENDDO
    ELSE IF ( mod(t,tw_da)<dt .AND. t-spin_up_etkf<=dt ) THEN
      DO i = 1, total
        Xam(i) = SUM(Xens(i,:))/ens_num
      END DO
    ENDIF

    ! generate analysis
    IF (mod(t,tw_da)<dt .AND. t-spin_up>dt .AND. t-spin_up_etkf>dt) THEN
      N_eff = 1.d0/SUM(W*W)
      !PRINT *, "DEBUG : Neff", 1.d0/SUM(W*W)
      IF ( N_eff > 1.1d0*ens_num*part_num ) THEN
        !CALL etkf_w( total, total, ens_num, Xens, W_tilde, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)
        CALL pf(total, ndim_dr, ens_num, part_num, Xdr_part, W, luse, yobs(:,cnt_obs+1), R, W_a)
        W = W_a
        W_tilde = SUM(W,DIM=2)
        !Xens = Xens_a
        Xam = MATMUL(Xens, W_tilde)
      ELSE 
        !CALL etkf_w( total, total, ens_num, Xens, W_tilde, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)
        !W=1.d0/ens_num/part_num; W_tilde = SUM(W,DIM=2)
        !PRINT*, "Xbm_dr = ", Xbm(ndim+1:ndim+ndim_dr)
        Xbm = SUM(Xens,DIM=2)/ens_num
        CALL etkf( total, total, ens_num, Xens, luse, yobs(:,cnt_obs+1), R, infl, Xens_a, Xam, 105)
        CALL pf(total, ndim_dr, ens_num, part_num, Xdr_part, W, luse, yobs(:,cnt_obs+1), R, W_a)
        !CALL pf2(total, ndim_dr, ens_num, part_num, Xdr_part, W, Xam(ndim+1:ndim+ndim_dr), &
        !& luse, yobs(:,cnt_obs+1), R, R, W_a)
        !PRINT *, "DEBUG N_eff", 1.d0/SUM(W_a*W_a)
        W_tilde_a = SUM(W_a,DIM=2)
        W_tilde = W_tilde_a
        W = W_a
       ! PRINT *, "DEBUG W_tilde", W_tilde

        !CALL MH_resampling(total,ens_num,Xens_a,W_tilde,Xens,W_dummy)
        !Xens(1:2*natm,:) = Xens_a(1:2*natm,:)
        !Xens(1+ndim:ndim+ndim_dr,:) = Xens_a(ndim+1:ndim+ndim_dr,:)
        !CALL MH_resampling(2*noc+ndim_dr,ens_num,Xens_a(2*natm+1:ndim+ndim_dr,:),W_tilde,Xens(2*natm+1:ndim+ndim_dr,:),W_dummy)
        !CALL RSR_resampling_1D(2*noc+ndim_dr,ens_num,Xens_a(2*natm+1:ndim+ndim_dr,:),W_tilde,Xens(2*natm+1:ndim+ndim_dr,:),W_dummy)
        !CALL RSR_resampling_1D(ndim+ndim_dr,ens_num,Xens_a,W_tilde,Xens,W_dummy)
        !std_ens = 0.1*R
        !Xens(2*natm+1:ndim,:) = Xens_a(2*natm+1:ndim,:)
        !Xens(1+ndim:ndim+ndim_dr,:) = Xens_a(ndim+1:ndim+ndim_dr,:)
        !CALL MH_resampling(2*natm,ens_num,Xens_a(1:2*natm,:),W_tilde,Xens(1:2*natm,:),W_dummy)
        !CALL RSR_resampling_1D(2*natm,ens_num,Xens_a(1:2*natm,:),W_tilde,Xens(1:2*natm,:),W_dummy)
        !!$OMP PARALLEL DO
        !DO k = 1, ens_num
        !  CALL randn(ndim,err)
        !  Xens(1:2*natm,k) = Xens(1:2*natm,k) + err(1:2*natm)*std_ens(1:2*natm)*0.01
        !  Xens(2*natm+1:ndim,k) = Xens(2*natm+1:ndim,k) + err(2*natm+1:ndim)*std_ens(2*natm+1:ndim)*0.01
        !ENDDO
        !!$OMP END PARALLEL DO
        !Xam = MATMUL(Xens_a, W_tilde)

        !CALL MH_resampling(ndim_dr,ens_num,part_num,Xdr_part,W_a,Xdr_part_a,W)
        DO k = 1, ens_num
          W_a(k,:) = W_a(k,:)/SUM(W_a(k,:))
          CALL MH_resampling(ndim_dr,part_num,Xdr_part(:,k,:),W_a(k,:),Xdr_part_a(:,k,:),W(k,:))
          !$OMP PARALLEL DO
          DO j = 1, part_num
            CALL randn(ndim_dr,err_dr)
            Xdr_part(:,k,j) = Xdr_part_a(:,k,j) - Xens(ndim+1:ndim+ndim_dr,k) &
               & + Xens_a(ndim+1:ndim+ndim_dr,k)
            !Xdr_part(:,k,j) = Xdr_part_a(:,k,j) - Xbm(ndim+1:ndim+ndim_dr) &
            !  & + yobs(ndim+1:ndim+ndim_dr,cnt_obs+1)
            !PRINT *, "DEBUG: yobs_dr", yobs(ndim+1:ndim+ndim_dr,cnt_obs+1)
            !Xdr_part(:,k,j) = Xdr_part(:,k,j) + err_dr*R(ndim+1:ndim+ndim_dr)*0.0001
            Xdr_part(:,k,j) = Xdr_part(:,k,j) + err_dr*std_ens(ndim+1:ndim+ndim_dr)*0.01
            !Xdr_part(:,k,j) = Xam(ndim+1:ndim+ndim_dr) + err_dr*std_ens(ndim+1:ndim+ndim_dr)*0.1
        !!    Xdr_part(:,k,j) = Xdr_part(:,k,j) + err_dr*R(ndim+1:ndim+ndim_dr)*0.0001 
          ENDDO
          !$OMP END PARALLEL DO
        ENDDO
        Xens = Xens_a
        !Xdr_part = Xdr_part_a
        W = 1.d0/ens_num/part_num; W_tilde = SUM(W,DIM=2)
        !Xam(ndim+1:ndim+ndim_dr) = MATMUL(Xens(ndim+1:ndim+ndim_dr,:), W_tilde) 
         
      END IF
      !PRINT *, "DEBUG MH_resampling Xens(1,:): ", Xens(1,:)
    ELSE
      ! update Xam
      W_tilde = SUM(W,DIM=2)
      Xam(1:ndim) = MATMUL(Xens(1:ndim,:), W_tilde)
      !$OMP PARALLEL DO
      DO i = 1, ndim_dr
        Xam(ndim+i) = SUM(Xdr_part(i,:,:)*W)
      ENDDO
      !$OMP END PARALLEL DO
    END IF

    IF (mod(t,tw)<dt) THEN
    ! Write ensemble spread
      DO i = 1, total
        std_ens(i) = 0.d0
        DO j = 1, ens_num
          std_ens(i) = std_ens(i) + ((Xens(i,j)-Xam(i))*SQRT(W_tilde(j)))*((Xens(i,j)-Xam(i))*SQRT(W_tilde(j)))
        END DO
        std_ens(i) = SQRT(std_ens(i))
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



END PROGRAM hybrid_maooam
