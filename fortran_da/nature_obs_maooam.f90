PROGRAM nature_obs_maooam

  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, tw_solo
  USE params, only: natm, noc ! LUYU: add for solo 
  USE params, only: tw_da
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_mt, only: randn
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y_obs
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up

  INTEGER :: i,j,k,n

  PRINT*, 'Model MAOOAM v1.3'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  t_up=dt/t_trans*100.D0

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  X=IC

  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)

     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat

  IF (writeout) THEN
      OPEN(100,file='nature.dat')
      OPEN(102,file='yobs.dat')
  END IF
  
  t=0.D0
  t_up=dt/t_run*100.D0

  ! Read R matrix Climatology
  ALLOCATE(R(ndim))
  ALLOCATE(err(ndim))
  ALLOCATE(y_obs(2*natm))
  PRINT*, 'Read R matrix....'
  OPEN(103,file="fort.202",action="read",form="formatted",access="sequential")
  DO n = 1, ndim
    READ(103,"(10000(D24.17,1x))") R(n) 
    !print*, "n, R(n)=", n, R(n)
  ENDDO
  CLOSE(103)
  PRINT*, 'Finish reading R matrix...'

  IF (writeout) THEN
      WRITE(100,*) t,X(1:ndim)
      CALL randn(ndim,err)
      DO n = 1,2*natm
        y_obs(n) = X(n) + err(n)*R(n)*0.1 ! Use 10% climatology.
        !print*, y_obs(n)
      END DO
      WRITE(102,*) t,y_obs(1:2*natm)
  END IF

  DO WHILE (t<t_run)

     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t,tw)<dt) THEN 
        IF (writeout) WRITE(100,*) t,X(1:ndim)
        CALL acc(X)

        CALL randn(ndim,err)
        DO n = 1,2*natm
          y_obs(n) = X(n) + err(n)*R(n)*0.1 ! Use 10% climatology
        END DO
        WRITE(102,*) t,y_obs(1:2*natm)      
     END IF

     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) THEN
      CLOSE(100)
      CLOSE(102)
  END IF

END PROGRAM nature_obs_maooam
