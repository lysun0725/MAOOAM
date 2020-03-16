
! integrator.f90
!
!>  Module with the integration routines.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------

MODULE integrator_dr
  USE params, only: ndim
  USE params, only: oms, natm, noc, dr_num, dr_size, n, ndim_dr
  USE params, only: part_num
  USE tensor, only:sparse_mul3
  USE aotensor_def, only: aotensor
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0 !< Buffer to hold tendencies at the initial position
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1 !< Buffer to hold tendencies at the intermediate position

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1_dr !< Buffer to hold the intermediate position (Heun algorithm)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0_dr !< Buffer to hold tendencies at the initial position
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1_dr !< Buffer to hold tendencies at the intermediate position


  PUBLIC :: init_integrator_dr, step_dr

  INTERFACE step_dr
    module procedure step_dr_1D
    module procedure step_dr_2D
  END INTERFACE step_dr
CONTAINS
  
  !> Routine to initialise the integration buffers.
  SUBROUTINE init_integrator_dr
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim) ,STAT=AllocStat)
    IF (part_num > 1) THEN
      ALLOCATE(buf_y1_dr(1:ndim_dr*part_num),buf_f0_dr(1:ndim_dr*part_num),buf_f1_dr(1:ndim_dr*part_num) ,STAT=AllocStat)
    ELSE
      ALLOCATE(buf_y1_dr(1:ndim_dr),buf_f0_dr(1:ndim_dr),buf_f1_dr(1:ndim_dr) ,STAT=AllocStat)
    ENDIF
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_integrator_dr
  
  !> Routine computing the tendencies of the model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies

  SUBROUTINE tendencies_dr(t,ndr,y,ydr,res_dr)
    REAL(KIND=8), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: ndr
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(1:ndr), INTENT(IN) :: ydr
    REAL(KIND=8), DIMENSION(1:ndr), INTENT(OUT) :: res_dr
   
    CALL compute_vel(ndr,y,ydr,res_dr)
  END SUBROUTINE

  !> Routine to perform an integration step (Heun algorithm). The incremented time is returned.
  !> @param y Initial point.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE step_dr_1D(y,ydr,t,dt,res,res_dr)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(1:ndim_dr), INTENT(IN) :: ydr ! drifter initial states
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    REAL(KIND=8), DIMENSION(1:ndim_dr), INTENT(OUT) :: res_dr ! drifter final states
    
    CALL tendencies(t,y,buf_f0)
    CALL tendencies_dr(t,ndim_dr,y,ydr,buf_f0_dr)
    buf_y1 = y+dt*buf_f0
    buf_y1_dr = ydr + dt*buf_f0_dr

    CALL tendencies(t+dt,buf_y1,buf_f1)
    CALL tendencies_dr(t+dt,ndim_dr,buf_y1,buf_y1_dr,buf_f1_dr)
    res=y+0.5*(buf_f0+buf_f1)*dt
    res_dr=ydr+0.5*(buf_f0_dr+buf_f1_dr)*dt
    t=t+dt
  END SUBROUTINE step_dr_1D

  SUBROUTINE step_dr_2D(y,ydr,t,dt,res,res_dr)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(1:ndim_dr,part_num), INTENT(IN) :: ydr ! drifter initial states
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    REAL(KIND=8), DIMENSION(1:ndim_dr, part_num), INTENT(OUT) :: res_dr ! drifter final states
    REAL(KIND=8), ALLOCATABLE :: ydr_long(:), res_dr_long(:)
    INTEGER :: i,j,k

    ALLOCATE(ydr_long(ndim_dr*part_num))
    ALLOCATE(res_dr_long(ndim_dr*part_num))

    ! initialize  ydr_long
    !$OMP PARALLEL DO
    DO i = 1, part_num
      ydr_long((i-1)*ndim_dr+1:i*ndim_dr) = ydr(:,i)
    ENDDO 
    !$OMP END PARALLEL DO

    CALL tendencies(t,y,buf_f0)
    CALL tendencies_dr(t,ndim_dr*part_num,y,ydr_long,buf_f0_dr)
    buf_y1 = y+dt*buf_f0
    buf_y1_dr = ydr_long + dt*buf_f0_dr

    CALL tendencies(t+dt,buf_y1,buf_f1)
    CALL tendencies_dr(t+dt,ndim_dr*part_num,buf_y1,buf_y1_dr,buf_f1_dr)
    res=y+0.5*(buf_f0+buf_f1)*dt
    res_dr_long=ydr_long+0.5*(buf_f0_dr+buf_f1_dr)*dt

    ! reshape res_dr_long
    !$OMP PARALLEL DO
    DO i = 1, part_num
      res_dr(:,i) = res_dr_long((i-1)*ndim_dr+1:i*ndim_dr)
    ENDDO
    !$OMP END PARALLEL DO

    t=t+dt
  END SUBROUTINE step_dr_2D

  SUBROUTINE compute_vel(ndr,y,ydr,res_dr)
    INTEGER, INTENT(IN) :: ndr
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(1:ndr), INTENT(IN) :: ydr ! drifter initial states
    REAL(KIND=8), DIMENSION(1:ndr), INTENT(OUT) :: res_dr ! drifter velocity 
    REAL(KIND=8), DIMENSION(1:noc) :: psi_node
    INTEGER :: Ho, Po
    INTEGER :: i, j, k, cnt
    REAL(KIND=8) :: u, v, xp, yp
   
    cnt = 0
    psi_node = y(2*natm+1:2*natm+noc)
    
    !DO i=1,dr_num ! loop for the drifter number
    !  u=0.0; v=0.0
    !  xp=ydr((i-1)*dr_size+1)
    !  yp=ydr((i-1)*dr_size+2)
    !  DO j=1,noc
    !    Ho=oms(j,1); Po=oms(j,2)
    !    u=u-psi_node(j)*2*sin(Ho*n/2*xp)*cos(Po*yp)*Po
    !    v=v+psi_node(j)*2*cos(Ho*n/2*xp)*sin(Po*yp)*Ho*n/2        
    !  ENDDO
    !  res_dr((i-1)*dr_size+1)=u
    !  res_dr((i-1)*dr_size+2)=v
    !ENDDO 

    DO WHILE (cnt < ndr)
      u=0.0; v=0.0
      xp=ydr(cnt+1)
      yp=ydr(cnt+2)
      DO j=1,noc
        Ho=oms(j,1); Po=oms(j,2)
        u=u-psi_node(j)*2*sin(Ho*n/2*xp)*cos(Po*yp)*Po
        v=v+psi_node(j)*2*cos(Ho*n/2*xp)*sin(Po*yp)*Ho*n/2
      ENDDO
      res_dr(cnt+1)=u
      res_dr(cnt+2)=v
      cnt = cnt + dr_size
    ENDDO
     
  END SUBROUTINE compute_vel

END MODULE integrator_dr
