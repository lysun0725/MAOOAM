
!  maooam.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere 
!> model MAOOAM.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam 
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, tw_solo
  USE params, only: natm, noc ! LUYU: add for solo 
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_atm
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_ocn
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_solo  !< Combine X_atm and ocn states from X
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X_solo_new !< UPdated state solo variable
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_solo=0.D0
  REAL(KIND=8) :: t_up

  PRINT*, 'Model MAOOAM v1.3'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
      OPEN(10,file='evol_field.dat')
      OPEN(11,file='evol_field_solo.dat')
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  ! LUYU: allocate X_solo, X_atm
  ALLOCATE(X_solo(0:ndim),X_solo_new(0:ndim))
  ALLOCATE(X_atm(2*natm),X_ocn(2*noc))

  X=IC

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

  ! LUYU: initialization of X_solo and X_atm
  X_atm = X(1:2*natm)
  X_ocn = X(2*natm+1:ndim)
  X_solo = (/ X(0), X_atm, X_ocn/)

  PRINT*, 'tw_solo = ', tw_solo

  IF (writeout) THEN
      WRITE(10,*) t,X(1:ndim)
      WRITE(11,*) t,X_solo(1:ndim)
  END IF

  DO WHILE (t<t_run)

     CALL step(X_solo, t_solo, dt, X_solo_new)
     X_atm = X_solo_new(1:2*natm)
     X_solo = (/X_solo_new(0), X_atm, X_ocn/)

     CALL step(X,t,dt,Xnew)
     X=Xnew

     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        CALL acc(X)
 
        IF (writeout) WRITE(11,*) t_solo,X_solo(1:ndim)
        !CALL acc(X_solo)      
     END IF

     ! LUYU: update forced data
     IF (mod(t,tw_solo)<dt) X_ocn = X(2*natm+1:ndim)

     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) THEN
      CLOSE(10)
      CLOSE(11)
  END IF

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam 
