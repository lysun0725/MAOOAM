module m_mt
  implicit none

  integer,parameter :: r_single = 4
  integer,parameter :: r_kind   = kind(0.d0)
  integer,parameter :: i_kind   = 4

  public :: EIGEN
  public :: EYE
  public :: QR
  public :: INV
  public :: RENORM
  public :: SVD
  public :: RANDN    ! Guassian noise N(0,1) generateor
  public :: SHUFFLE

  interface RENORM
    module procedure RENORM1d
    module procedure RENORM2d
  endinterface

  interface RANDN
    module procedure RANDN1d
    module procedure RANDN2d
  endinterface

contains 

subroutine EIGEN( n, A, eigvect, eigval)
!
! check "dsyev" in the LAPACK
!
implicit none
! passed args
  integer,   intent(in   ) :: n
  real(8),intent(in   ) :: A(n,n)
  real(8),intent(  out) :: eigvect(n,n)
  real(8),intent(  out) :: eigval(n)

! local vars
  real(8)         :: r8A(n,n)
  real(8)         :: r8eigvect(n,n)
  real(8)         :: r8eigval(n)
  integer           :: lwk
  integer,parameter :: lwkmax = 3000 
  real(8)         :: wk(lwkmax)
  integer           :: i, ierr

  ierr = 0
  r8A      = REAL(A, kind(0.0d0))
  r8eigvect = 0.0d0
  r8eigval = 0.0d0
  wk       = 0.0d0

  if ( 3*n-1 > lwkmax ) stop "ERROR in EIGEN: increase lwkmax"
! get the optimal work space
  lwk = -1
  call dsyev( 'V', 'U', n, r8A, n, r8eigval, wk, lwk, ierr )
  if (ierr/=0) stop "ERROR in EIGEN: dsyev"
  lwk = MIN( lwkmax, INT(wk(1)) )
! solve eigenproblem
  call dsyev( 'V', 'U', n, r8A, n, r8eigval, wk, lwk, ierr )
  if ( ierr /= 0 ) stop "ERROR in EIGEN: dsyev"
! put the largest eigenvalue as the 1st column, which brings the largest variance
  do i = 1, n
     eigval(i)    = r8eigval(n+1-i)
     eigvect(:,i) = r8A(:,n+1-i)
  enddo

endsubroutine


!
! generate random 1d-array
!
subroutine RANDN1d(n, A)
  implicit none
 
  integer,intent(in) :: n
  real(r_kind),intent(out) :: A(n)

  integer :: i
  real(r_kind) :: am

  do i = 1, n
     A(i) = real(rnorm(),r_kind)
  enddo
  am = sum(A(:))/n
  A(:) = A(:)-am

endsubroutine

!
! random 2D-array
!
subroutine RANDN2d(n,m,A)
  implicit none

  integer,intent(in) :: n, m
  real(r_kind),intent(out) :: A(n,m)

  real(r_kind) :: Am
  integer :: i, j

  do i = 1, n
     do j = 1, m
        A(i,j) = real(rnorm(),r_kind)
     enddo
  enddo
  Am = sum(A(:,:))/(n*m)
  A = A - Am

endsubroutine

!
! SVD: A= U * S * V^t
!
subroutine SVD(n,m,A,U,S,VT)
  implicit none

  integer,intent(in) :: n, m
  real(r_kind),intent(in)  :: A(m,n)
  real(r_kind),intent(out) :: U(m,m), VT(n,n)
  real(r_kind),intent(inout) :: S(:)

  character :: jobu, jobvt
  integer :: lda, ldu, ldvt
  integer :: info
  integer :: lwork
  real(r_kind) :: bufA(m,n)
  real(r_kind),allocatable :: work(:)

  bufA=A
  lda=m; ldu=m; ldvt=n
  lwork=3*min(m,n) + max ( max(m,n), 2*min(m,n) )
  allocate( work(lwork) )
  jobu='A'; jobvt='A'
  call dgesvd( jobu, jobvt, m, n, bufA, lda, S, U, ldu, VT, ldvt, work, lwork, info)
  if (info/=0) stop 'svd: error in dgesvd'
  deallocate( work )

endsubroutine

!
! identity matrix
!
subroutine EYE(n,A)
  implicit none

  integer,intent(in) :: n
  real(r_kind),intent(out) :: A(n,n)

  integer :: i

  A=0.d0
  do i = 1, n
     A(i,i) = 1.d0
  enddo

endsubroutine

!
! A random upper triangle square matrix
!
subroutine RANDUTRI(n,A)
  implicit none
  integer,intent(in) :: n
  real(r_kind),intent(out) :: A(n,n)

  integer :: i, j

  A=0.d0
  do i = 1, n
     do j = i,n
        A(i,j) = REAL(rnorm(),8)
     enddo
  enddo

  call RENORM(n,n,A)

endsubroutine


subroutine RENORM1d(n,A)
  implicit none

  integer,intent(in) :: n
  real(r_kind),intent(inout) :: A(n)

  real(r_kind) :: rn2

  rn2 = sqrt(sum(A(:)*A(:)))
  A(:) = A(:)/rn2

endsubroutine



!
! renormalize each column of the general matrix
!
subroutine RENORM2d(n,m,A)
  implicit none

  integer,intent(in) :: n, m
  real(r_kind),intent(inout) :: A(n,m)

  integer :: i
  real(r_kind) :: rn2

  do i = 1, m
     rn2 = sqrt(sum(A(:,i)*A(:,i)))
     A(:,i) = A(:,i)/rn2
  enddo

endsubroutine

!
! inverse of a general square matrix
!
subroutine INV(n,A,invA)
  implicit none

  integer,intent(in) :: n
  real(r_kind),intent(in) :: A(n,n)
  real(r_kind),intent(out) :: invA(n,n)

  integer :: info, ipiv(n)
  integer :: lwork, lda
  real(r_kind) :: work(n)

  invA = A
  lda  = n; lwork = n
  call dgetrf(n,n,invA,lda,ipiv,info)
  if (info/=0) stop "ERROR in INV: dgetrf"
  call dgetri(n,invA,lda,ipiv,work,lwork,info)
  if (info/=0) stop "ERROR in INV: dgetri"
endsubroutine

!
! QR decomposition of a square matrix
!
subroutine QR(n,A,Q,R)
  implicit none

  integer,intent(in) :: n
  real(r_kind),intent(in) :: A(n,n) 
  real(r_kind),intent(out) :: Q(n,n), R(n,n)

  integer :: ierr
  real(r_kind) :: work(n), bufA(n,n), tau(n)
  integer :: lda
  integer :: i

  ierr=0
  lda = n
  bufA = A
  call dgeqrf(n,n,bufA,lda,tau,work,n,ierr)
  if (ierr/=0) stop "ERROR in QR: dgetqrf"
  R=0.d0
  do i = 1, n
     R(i,i:n) = bufA(i,i:n)
  enddo
  Q=bufA
  call dorgqr(n,n,n,Q,lda,tau,work,n,ierr)
  if (ierr/=0) stop "ERROR in QR: qorgqr"

endsubroutine

!
! random number generator
!
subroutine rnorm1d( nn, r1d )
  implicit none

  integer,intent(in) :: nn
  real(r_single),intent(inout) :: r1d(nn)
  
  integer :: i

  do i =1 , nn
     r1d(i) = rnorm()
  enddo

endsubroutine 

real(r_single) function rnorm() result( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

implicit none

! Local variables

real(r_single)            :: u, sum
real(r_single), save      :: v, sln
logicaL, save           :: second = .false.
real(r_single), parameter :: one = 1.0_r_single, vsmall = TINY( one )

if (second) then
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

else
! First call; generate a pair of random normals

  second = .true.
  do
    call random_number( u )
    call random_number( v )
    u = scale( u, 1 ) - one
    v = scale( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) exit
  end do
  sln = sqrt(- scale( log(sum), 1 ) / sum)
  fn_val = u*sln
end if

return
end function rnorm

subroutine iran(l,n,ran_int)
! generate an array of N random integers between 0 and L-1.
 real(r_single) :: rnd
 integer(i_kind) :: L,i,N
 integer(i_kind) :: ran_int(N)
 do i = 1,n
    call random_number(rnd)
    ran_int(i) = nint(float(l-1)*rnd)
 end do
end subroutine iran

subroutine set_random_seed ( iseed , myrank)
!
!*******************************************************************************
!
!! SET_RANDOM_SEED initializes the FORTRAN 90 random number generator.
!
!
!  Discussion:
!
!    If ISEED is nonzero, then that value is used to construct a seed.
!
!    If ISEED is zero, then the seed is determined by calling the date 
!    and time routine.  Thus, if the code is run at different times, 
!    different seed values will be set.
!
!  Parameters:
!
!    Input, integer ISEED, is nonzero for a user seed, or 0 if the
!    seed should be determined by this routine.
!
  implicit none
!
  integer(i_kind), intent(in), optional :: myrank
  integer(i_kind) date_time(8)
  logical, parameter :: debug = .false.
  integer(i_kind) i,j,k
  integer(i_kind) iseed
  integer(i_kind), allocatable :: seed(:)

!
!  Initialize the random seed routine.
!
  call random_seed
!
!  Request the size of a typical seed.
!  (On the DEC ALPHA, K is returned as 2.)
!
  call random_seed ( size = k )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SET_RANDOM_SEED:'
    write ( *, '(a,i6)' ) '  Random seed size is K = ', k
  end if
!
!  Set up space for a seed vector.
!
  allocate ( seed(k) )

  if ( iseed /= 0 ) then

    seed(1:k) = iseed

  else
!
!  Make up a "random" value based on date and time information.
!
    call date_and_time ( values = date_time )

    do i = 1, k

      seed(i) = 0

      do j = 1, 8
        if (present(myrank)) then
        seed(i) = seed(i) + ( j + i ) * date_time(j) + myrank * 100
        else
        seed(i) = seed(i) + ( j + i ) * date_time(j)
        endif
        seed(i) = ishftc ( seed(i), 4 * ( j - 1 ) )
      end do

    end do

  end if

  if  ( debug ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SET_RANDOM_SEED:'
    write ( *, '(a)' ) '  The random seed vector:'
    write ( *, '(a)' ) ' '

    do i = 1, k
      write ( *, '(i12)' ) seed(i)
    end do

  end if
!
!  Send this random value back to the RANDOM_SEED routine, to be
!  used as the seed of the random number generator.
!
  call random_seed ( put = seed(1:k) )

  deallocate ( seed )

end subroutine set_random_seed


subroutine SHUFFLE(a)
    integer, intent(inout) :: a(:)
    integer :: i, randpos, temp
    real :: r

    do i = size(a), 2, -1
      call random_number(r)
      randpos = int(r * i) + 1
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do

end subroutine SHUFFLE

endmodule

!
! test SVD
!
!program main
!  use m_mt
!  implicit none
!
!  integer :: n, m
!  real(8),allocatable :: A(:,:), U(:,:),S(:),VT(:,:)
!
!  integer :: i, j
!
!  m=3; n=4
!  allocate(A(m,n),U(m,m),S(max(m,n)),VT(n,n))
!  A(1,:)=(/1.d0, 2.d0, 3.d0, 4.d0/)
!  A(2,:)=(/3.d0, 4.d0, 5.d0, 6.d0/)
!  A(3,:)=(/6.d0, 5.d0, 4.d0, 1.d0/)
!  call SVD(n,m,A,U,S,VT)
!  print*, "A="
!  do i = 1, m
!     write(*,*) (A(i,j),j=1,n)
!  enddo
! 
!  print*, "U="
!  do i = 1, m
!     write(*,*) (U(i,j),j=1,m)
!  enddo
!
!  print*, "U=", (S(i),i=1,max(m,n))
!
!  print*, "VT="
!  do i = 1, n
!     write(*,*) (VT(i,j),j=1,n)
!  enddo
!endprogram

!
! test randn
!
!program main
!  use m_mt, only :randn
!  implicit none
!
!  real(8) :: a(100)
!  integer :: i , k 
!
!  do i = 1, 3
!    call randn(100,a)
!    do k = 1 , 5
!       write(*,*) a(k)
!    enddo
!    pause "i"
!  enddo
! 
!endprogram

