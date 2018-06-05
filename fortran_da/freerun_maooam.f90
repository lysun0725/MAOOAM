PROGRAM freerun_maooam
  
  USE params, only: ndim, dt, tw, t_run,writeout, tw_solo
  USE params, only: natm, noc
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_mt,     only: randn
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_atm
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_ocn
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_solo  !< Combine X_atm and ocn states from X
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_solo_new !< UPdated state solo variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_solo=0.D0
  REAL(KIND=8) :: t_up
  CHARACTER(len=24) :: soloname

  INTEGER :: i,j,k,n

  PRINT*, 'Model MAOOAM v1.3 for free run'
  PRINT*, 'Loading information...'

  CALL init_aotensor
  CALL init_integrator

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(X_solo(0:ndim),X_solo_new(0:ndim))
  ALLOCATE(X_atm(2*natm),X_ocn(2*noc))

  OPEN(100,file='nature.dat',action="read",access="sequential")
  READ(100,*) X(0:ndim)
  X(0) = 1.D0
  CLOSE(100)

  IF (writeout) THEN
      IF (tw_solo .lt. 1) THEN
          write(soloname,'("freerun_atm",F0.3,".dat")') tw_solo
      ELSE
          write(soloname,'("freerun_atm",I4.4,".dat")') INT(tw_solo)
      ENDIF
      OPEN(101,file=soloname)
  END IF

  ! Read R matrix Climatology
  ALLOCATE(R(ndim))
  ALLOCATE(err(ndim))
  PRINT*, 'Read R matrix....'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO n = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(n) 
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

  ! Add noise to X_atm and initialize X_solo and X_ocn
  CALL randn(ndim,err)
  DO n = 1,2*natm
     X_atm(n) = X(n) + err(n)*R(n)
  END DO
  X_ocn = X(2*natm+1:ndim)
  X_solo = (/ X(0), X_atm, X_ocn/) 

    PRINT*, 'Starting the time evolution...'
  PRINT*, 'tw_solo = ', tw_solo

  CALL init_stat
  
  t=0.D0
  t_solo=0.D0
  t_up=dt/t_run*100.D0

  IF (writeout) THEN
    WRITE(101,*) t,X_solo(1:ndim)
  END IF

  DO WHILE (t<t_run)

     CALL step(X_solo, t_solo, dt, X_solo_new)
     X_atm = X_solo_new(1:2*natm)
     X_solo = (/X_solo_new(0), X_atm, X_ocn/)

     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(101,*) t_solo,X_solo(1:ndim)
        CALL acc(X_solo)      
     END IF

     ! LUYU: update forced data
     IF (mod(t,tw_solo)<dt) THEN
        X_ocn = X(2*natm+1:ndim)
        X_solo(2*natm+1:ndim) = X_ocn ! Update the ocn part in X_sol
     END IF

     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) THEN
    CLOSE(101)
  END IF
END PROGRAM
