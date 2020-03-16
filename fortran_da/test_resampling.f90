PROGRAM test_resampling
  USE m_da_maooam, only: MH_resampling, RSR_resampling_1D

  IMPLICIT NONE

  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Xens_a, Xens
  REAL(KIND=8),DIMENSION(:), ALLOCATABLE :: W_tilde, W_tilde2, W_dummy
  INTEGER :: total, ens_num

  total = 1; ens_num = 6 
  ALLOCATE(Xens_a(total,ens_num),Xens(total,ens_num))
  ALLOCATE(W_tilde(ens_num),W_tilde2(ens_num),W_dummy(ens_num))

  W_tilde = (/0.1, 0.1, 0.2, 0.1, 0.1, 0.4/)
  W_tilde2 = (/0.2, 0.1, 0.3, 0.05, 0.05, 0.5/)
  print *, "W_tilde = ", W_tilde
  print *, "W_tilde2 = ", W_tilde2
  print *, "MAX ARRAY = ", MAX(W_tilde,W_tilde2)

  Xens_a = reshape((/ 1, 2, 3, 4, 5, 6 /), shape(Xens_a))
  

  !CALL MH_resampling(total,ens_num,Xens_a,W_tilde,Xens,W_dummy)
  CALL RSR_resampling_1D(total,ens_num,Xens_a,W_tilde,Xens,W_dummy)

  PRINT*, "TEST RESAMPLING: Xens = ", Xens

END PROGRAM test_resampling
