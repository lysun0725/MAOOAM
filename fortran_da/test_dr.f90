
! test_dr.f90
!
!> Test for drifter forecast model of MAOOAM
!
!> 2019 Luyu Sun


PROGRAM test_dr
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, tw_solo
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_dr    !< LUYU: add drifter states
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew_dr !< LUYU: add updated drifter state variables
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up

  PRINT*, 'Model MAOOAM v1.3 with drifter simulator'
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
 
  IF (writeout) THEN
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim)); X=IC
  IF (do_drifter) THEN
    ALLOCATE(X_dr(ndim_dr),Xnew_dr(ndim_dr)); X_dr=IC_DR
  END IF

  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution...'
 
  CALL init_stat

  t=0.D0
  t_up=dt/t_run*100.D0

  OPEN(11,file='evol_field_dr.dat')
  IF (writeout) THEN
    IF (do_drifter) THEN
      WRITE(11,*) t, X(1:ndim), X_dr
    ELSE
      WRITE(11,*) t, X(1:ndim) 
    END IF
  END IF

  DO WHILE (t<t_run)
    IF (do_drifter) THEN
      CALL step_dr(X,X_dr,t,dt,Xnew,Xnew_dr)
      X=Xnew; X_dr=Xnew_dr
    ELSE
      CALL step(X,t,dt,Xnew)
      X=Xnew
    END IF


    IF (mod(t,tw)<dt .AND. writeout) THEN
      IF (do_drifter) THEN
        WRITE(11,*) t,X(1:ndim),X_dr
      ELSE
        WRITE(11,*) t,X(1:ndim)
      END IF
    END IF

    IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'
  CLOSE(11)
END PROGRAM test_dr
