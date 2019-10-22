PROGRAM freerun_maooam
  
  USE params, only: ndim, dt, tw, t_run,writeout, tw_solo
  USE params, only: natm, noc
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_mt,     only: randn

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size, oms, n
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X           !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_dr        !< LUYU: add drifter states
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_free      !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_free_dr   !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew_free   !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew_free_dr!< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R_dr
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err_dr
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up
  CHARACTER(len=24) :: freename

  INTEGER :: i,j,k,nn, m
  INTEGER :: Ho, Po

  PRINT*, 'Model MAOOAM v1.3 for free run'
  PRINT*, 'Loading information...'

  CALL init_aotensor
  IF (do_drifter) THEN
    CALL init_integrator_dr
  ELSE
    CALL init_integrator
  END IF

  ALLOCATE(X(0:ndim),X_free(0:ndim),Xnew_free(0:ndim))
  IF (do_drifter) THEN
    ALLOCATE(X_dr(ndim_dr),X_free_dr(ndim_dr),Xnew_free_dr(ndim_dr))
  END IF 

  OPEN(100,file='nature.dat',action="read",access="sequential")
  READ(100,*) X(0:ndim),X_dr(1:ndim_dr)
  X(0) = 1.D0
  PRINT*, X(0:ndim), X_dr
  CLOSE(100)

  IF (writeout) THEN
      write(freename,'("freerun.dat")')
      OPEN(101,file=freename)
  END IF

  ! Read R matrix Climatology
  ALLOCATE(R(ndim))
  ALLOCATE(err(ndim))
  IF (do_drifter) THEN
      ALLOCATE(R_dr(ndim_dr),err_dr(ndim_dr))
  END IF
  PRINT*, 'Read R matrix....'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO nn = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(nn) 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

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
        R_dr((nn-1)*dr_num+1) = R_dr(1)
        R_dr((nn-1)*dr_num+2) = R_dr(2)
      ENDDO
  END IF

  ! Add noise to X_atm and initialize X_solo and X_ocn
  CALL randn(ndim,err)
  DO nn = 1,ndim
     X_free(nn) = X(nn) + err(nn)*R(nn)
  END DO
  X_free(0) = X(0)
  IF (do_drifter) THEN
    CALL randn(ndim_dr,err_dr)
    DO nn = 1, ndim_dr
      X_free_dr(nn) = X_dr(nn) + err_dr(nn)*R_dr(nn)*0.1
      PRINT*, "X_dr =", X_dr(nn), "X_free_dr =", X_free_dr(nn)
    END DO
  END IF

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0
  t_up=dt/t_run*100.D0

  IF (writeout) THEN
    IF (do_drifter) THEN
      WRITE(101,*) t,X_free(1:ndim),X_free_dr
    ELSE
      WRITE(101,*) t,X_free(1:ndim)
    END IF
  END IF

  DO WHILE (t<t_run)

     IF (do_drifter) THEN
       CALL step_dr(X_free,X_free_dr,t,dt,Xnew_free,Xnew_free_dr)
       X_free = Xnew_free; X_free_dr = Xnew_free_dr
     ELSE
       CALL step(X_free, t, dt, Xnew_free)
       X_free = Xnew_free
     END IF

     IF (mod(t,tw)<dt) THEN
        IF (writeout) THEN
          IF (do_drifter) THEN
            WRITE(101,*) t,X_free(1:ndim),X_free_dr
          ELSE
            WRITE(101,*) t,X_free(1:ndim)
          END IF
          CALL acc(X_free)
        END IF       
     END IF

     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) THEN
    CLOSE(101)
  END IF
END PROGRAM
