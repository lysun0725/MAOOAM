PROGRAM obs_phys
  USE params, only: ndim, dt, tw, t_run,writeout
  USE params, only: natm, noc
  USE params, only: t_run_da, tw_da, ens_num, nobs, infl, ini_err, spin_up
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  USE m_da_maooam, only: etkf
  USE m_mt,     only: randn

  USE params, only: do_drifter, ndim_dr, dr_num, dr_size, oms, n, nobs_dr
  USE IC_def, only: load_IC_dr, IC_DR
  USE integrator_dr, only: init_integrator_dr,step_dr

  USE params, only: obs_atm_u, obs_atm_v, obs_atm_t
  USE params, only: obs_ocn_u, obs_ocn_v, obs_ocn_t
  USE params, only: obs_drf_x, obs_drf_y, obs_drf_t
  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X ! Driving system
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_dr
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_dr_old
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X_dr_vel
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: wk
  REAL(KIND=8) :: t
  INTEGER :: i,j,k,l,m,nn
  CHARACTER(LEN=21) :: drfname = 'yobs_drf_vel.dat.orig'
  CHARACTER(LEN=16) :: drfname2 = 'yobs_drf_vel.dat'
  CHARACTER(LEN=16) :: atmname = 'yobs_atm_phy.dat'
  CHARACTER(LEN=16) :: ocnname = 'yobs_ocn_phy.dat'
  CHARACTER(LEN=5) :: pathname
  CALL init_aotensor

  ALLOCATE(X(ndim),X_dr(ndim_dr),X_dr_old(ndim_dr),wk(ndim+ndim_dr+1))
  ALLOCATE(X_dr_vel(ndim_dr))

  WRITE(pathname,'("d",I0.3,"/")') dr_num
  OPEN(100,file=pathname//'nature.dat',action="read",access="sequential")
  ! Read initial condition when t=0.d0
  READ(100,*) wk(1:ndim+ndim_dr+1)
  t = wk(1); X = wk(2:ndim+1); X_dr_old = wk(ndim+2:ndim+ndim_dr+1)
  OPEN(300,file=drfname,action="write")
  DO WHILE (t+tw-t_run <= dt)
    READ(100,*) wk
    t = wk(1); X = wk(2:ndim+1); X_dr = wk(ndim+2:ndim+ndim_dr+1)
    
    IF (dr_size == 2) THEN
      X_dr_vel = (X_dr - X_dr_old)/tw
      DO i = 1, dr_num
        WRITE(300,*) t, obs_ocn_u, X_dr((i-1)*dr_size+1), X_dr((i-1)*dr_size+2), X_dr_vel((i-1)*dr_size+1)
        WRITE(300,*) t, obs_ocn_v, X_dr((i-1)*dr_size+1), X_dr((i-1)*dr_size+2), X_dr_vel((i-1)*dr_size+2)
      END DO
      X_dr_old = X_dr
    END IF
  END DO
  CLOSE(300)
  CLOSE(100)

  wk = 0.d0; t=0.d0; X_dr_old=0.d0; X_dr=0.d0; X_dr_vel=0.d0
  OPEN(200,file=pathname//'yobs_drf.dat',action="read",access="sequential")
  ! Read initial condition when t=0.d0
  READ(200,*) wk(1:ndim_dr+1)
  t = wk(1); X_dr_old = wk(2:ndim_dr+1)
  OPEN(400,file=pathname//drfname2,action="write")
  DO WHILE (t+tw-t_run <= dt)
    READ(200,*) wk(1:ndim_dr+1)
    t = wk(1);   X_dr = wk(2:ndim_dr+1)
    
    IF (dr_size == 2) THEN
      X_dr_vel = (X_dr - X_dr_old)/tw
      DO i = 1, dr_num
        WRITE(400,*) t, obs_ocn_u, X_dr((i-1)*dr_size+1), X_dr((i-1)*dr_size+2), X_dr_vel((i-1)*dr_size+1)
        WRITE(400,*) t, obs_ocn_v, X_dr((i-1)*dr_size+1), X_dr((i-1)*dr_size+2), X_dr_vel((i-1)*dr_size+2)
      END DO
      X_dr_old = X_dr
    END IF
  END DO
  CLOSE(400)
  CLOSE(200)  
END PROGRAM obs_phys
