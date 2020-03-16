PROGRAM nature_obs_maooam

  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE params, only: natm, noc ! LUYU: add for solo 
  USE params, only: tw_da
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_mt, only: randn

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size, oms, n
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_dr    !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew_dr !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R_dr
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err_dr
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y_obs_atm, y_obs_ocn
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y_obs_drf
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=5) :: pathname 
  INTEGER :: i,j,k,nn,m,Ho,Po

  PRINT*, 'Model MAOOAM v1.3'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition
  IF (do_drifter) THEN
    CALL load_IC_dr       ! LUYU
  END IF

  IF (do_drifter) THEN
    CALL init_integrator_dr
  ELSE
    CALL init_integrator  ! Initialize the integrator 
  END IF

  t_up=dt/t_trans*100.D0

  ALLOCATE(X(0:ndim),Xnew(0:ndim)); X=IC
  IF (do_drifter) THEN
    ALLOCATE(X_dr(1:ndim_dr),Xnew_dr(1:ndim_dr)); X_dr=IC_DR
  END IF

  PRINT*, 'Starting the transient time evolution...'
  PRINT*, 'IC_DR=',IC_DR
  DO WHILE (t<t_trans)

     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat

  IF (writeout) THEN
      WRITE(pathname,'("d",I0.3,"/")') dr_num
      OPEN(100,file=pathname//'nature.dat')
      OPEN(102,file=pathname//'yobs_atm.dat')
      OPEN(104,file=pathname//'yobs_ocn.dat')
      IF (do_drifter) OPEN(105,file=pathname//'yobs_drf.dat')
  END IF
  
  t=0.D0
  t_up=dt/t_run*100.D0

  ! Read R matrix Climatology
  ALLOCATE(R(ndim))
  ALLOCATE(err(ndim))
  ALLOCATE(y_obs_atm(2*natm),y_obs_ocn(2*noc))
  IF (do_drifter) THEN
      ALLOCATE(R_dr(ndim_dr),err_dr(ndim_dr),y_obs_drf(ndim_dr))
  END IF
  PRINT*, 'Read R matrix....'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO nn = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(nn) 
    print*, "n, R(n)=", nn, R(nn)
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

  IF (writeout) THEN
      IF (do_drifter) THEN
        WRITE(100,*) t,X(1:ndim),X_dr
      ELSE
        WRITE(100,*) t,X(1:ndim)
      END IF
      CALL randn(ndim,err)
      DO nn = 1,2*natm
        y_obs_atm(nn) = X(nn) + err(nn)*R(nn)*0.1 ! Use 10% climatology.
        !print*, y_obs(n)
      END DO
      DO nn = 1,2*noc
        y_obs_ocn(nn) = X(nn+2*natm) + err(nn+2*natm)*R(nn+2*natm)*0.1 ! Use 10% climatology.
        !print*, y_obs(n)
      END DO
      WRITE(102,*) t,y_obs_atm(1:2*natm)
      WRITE(104,*) t,y_obs_ocn(1:2*noc)
  END IF

  R_dr = 0.d0
  IF (do_drifter) THEN
      ! Define R_dr
      DO m = 1,noc
        Ho=oms(m,1); Po=oms(m,2)
        R_dr(1) = R_dr(1) + (R(m+2*natm)**2)*(Po**2)/noc
        R_dr(2) = R_dr(2) + (R(m+2*natm)**2)*((Ho*n/2)**2)/noc 
      ENDDO

      R_dr(1)=SQRT(R_dr(1)); R_dr(2)=SQRT(R_dr(2))
      PRINT*, "R_dr(1) = ", R_dr(1), "R_dr(2) = ", R_dr(2)
   
      DO nn = 2,dr_num
        R_dr((nn-1)*dr_size+1) = R_dr(1)
        R_dr((nn-1)*dr_size+2) = R_dr(2)
      ENDDO

      !write the initial observation of the drifter
      CALL randn(ndim_dr,err_dr)
      DO nn = 1, ndim_dr
        y_obs_drf(nn) = X_dr(nn) + err_dr(nn)*R_dr(nn)*0.1
        PRINT*, "X_dr =", X_dr(nn), "y_obs_drf =", y_obs_drf(nn)
      ENDDO
      WRITE(105,*) t,y_obs_drf(1:ndim_dr)
  END IF

  DO WHILE (t<t_run)
     IF (do_drifter) THEN
       CALL step_dr(X,X_dr,t,dt,Xnew,Xnew_dr)
       X=Xnew; X_dr=Xnew_dr
     ELSE
       CALL step(X,t,dt,Xnew)
       X=Xnew
     END IF

     IF (mod(t,tw)<dt) THEN 
        IF (writeout) THEN
          IF (do_drifter) THEN
            WRITE(100,*) t,X(1:ndim),X_dr
          ELSE 
            WRITE(100,*) t,X(1:ndim)
          END IF
        END IF
        CALL acc(X)

        CALL randn(ndim,err)
        DO nn = 1,2*natm
          y_obs_atm(nn) = X(nn) + err(nn)*R(nn)*0.1 ! Use 10% climatology
        END DO
        DO nn = 1,2*noc
          y_obs_ocn(nn) = X(nn+2*natm) + err(nn+2*natm)*R(nn+2*natm)*0.1 ! Use 10% climatology
        END DO
        WRITE(102,*) t,y_obs_atm(1:2*natm)      
        WRITE(104,*) t,y_obs_ocn(1:2*noc)      
        
        IF (do_drifter) THEN
          CALL randn(ndim_dr,err_dr)
          DO nn = 1, ndim_dr
            y_obs_drf(nn) = X_dr(nn) + err_dr(nn)*R_dr(nn)*0.1
          ENDDO
          WRITE(105,*) t,y_obs_drf(1:ndim_dr)
        END IF
     END IF

     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) THEN
      CLOSE(100)
      CLOSE(102)
      CLOSE(104)
      IF (do_drifter) CLOSE(105)
  END IF

END PROGRAM nature_obs_maooam
