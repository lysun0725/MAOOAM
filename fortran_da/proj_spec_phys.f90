MODULE proj_spec_phys
  USE params, only: ndim
  USE params, only: oms, natm, noc, dr_num, dr_size, n, ndim_dr
  USE tensor, only:sparse_mul3
  USE aotensor_def, only: aotensor
  IMPLICIT NONE

  PUBLIC :: compute_ocn_uvel, compute_ocn_vvel

CONTAINS

SUBROUTINE compute_ocn_uvel(spec,pos,phys_u)
  REAL(KIND=8), DIMENSION(1:ndim), INTENT(IN) :: spec
  REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: pos !  states
  REAL(KIND=8), INTENT(OUT) :: phys_u ! drifter velocity 
  REAL(KIND=8), DIMENSION(1:noc) :: psi_node
  INTEGER :: Ho, Po
  INTEGER :: i, j, k
  REAL(KIND=8) :: u, xp, yp

  psi_node = spec(2*natm+1:2*natm+noc)

  u=0.0
  xp=pos(1); yp=pos(2)
  DO j=1,noc
    Ho=oms(j,1); Po=oms(j,2)
    u=u-psi_node(j)*2*sin(Ho*n/2*xp)*cos(Po*yp)*Po
  ENDDO
  phys_u=u

END SUBROUTINE compute_ocn_uvel

SUBROUTINE compute_ocn_vvel(spec,pos,phys_v)
  REAL(KIND=8), DIMENSION(1:ndim), INTENT(IN) :: spec
  REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: pos !  states
  REAL(KIND=8), INTENT(OUT) :: phys_v ! drifter velocity 
  REAL(KIND=8), DIMENSION(1:noc) :: psi_node
  INTEGER :: Ho, Po
  INTEGER :: i, j, k
  REAL(KIND=8) :: v, xp, yp

  psi_node = spec(2*natm+1:2*natm+noc)
  v=0.0
  xp=pos(1); yp=pos(2)
  DO j=1,noc
    Ho=oms(j,1); Po=oms(j,2)
    v=v+psi_node(j)*2*cos(Ho*n/2*xp)*sin(Po*yp)*Ho*n/2
  ENDDO
  phys_v=v

END SUBROUTINE compute_ocn_vvel

SUBROUTINE compute_ocn_temp(spec,pos,phys_t)
  REAL(KIND=8), DIMENSION(1:ndim), INTENT(IN) :: spec
  REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: pos !  states
  REAL(KIND=8), INTENT(OUT) :: phys_t ! drifter temperature 
  REAL(KIND=8), DIMENSION(1:noc) :: psi_node
  INTEGER :: Ho, Po
  INTEGER :: i, j, k
  REAL(KIND=8) :: t, xp, yp

  psi_node = spec(2*natm+noc+1:2*natm+2*noc)
  t=0.0
  xp=pos(1); yp=pos(2)
  DO j=1,noc
    Ho=oms(j,1); Po=oms(j,2)
    t=t+psi_node(j)*2*sin(Ho*n/2*xp)*sin(Po*yp)
  ENDDO
  phys_t=t

END SUBROUTINE compute_ocn_temp

END MODULE proj_spec_phys
