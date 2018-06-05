!
! a wrapper for MAOOAM
!
module m_maooam
  use params,           only : dt, t_trans, t_run, ndim
  use integrator,       only : step
  use tl_ad_integrator, only : get_tlm
  implicit none

  private

  public :: dt, dt_day, t_trans, t_run
  public :: ndim, nxa, nya, nxo, nyo
  public :: get_res_maooam, init_maooam
  public :: step, get_tlm
  public :: read_x0_maooam, write_le_maooam

  integer,save :: nxa, nya, nxo, nyo
  real(8),save :: dt_day
  character(80),save :: cshared


contains

subroutine get_res_maooam(nxa_in,nya_in,nxo_in,nyo_in)
  implicit none

  integer,intent(out) :: nxa_in, nya_in, nxo_in, nyo_in

  character(80) :: cin

  if (command_argument_count()<4) stop "should provide nxa_in, nya_in, nxo_in, nyo_in"
  call get_command_argument(1,cin); read(cin,"(I)") nxa_in
  call get_command_argument(2,cin); read(cin,"(I)") nya_in
  call get_command_argument(3,cin); read(cin,"(I)") nxo_in
  call get_command_argument(4,cin); read(cin,"(I)") nyo_in

endsubroutine

!
! should call first before anything except get_maooam_res()
!
subroutine init_maooam(nxa_in,nya_in,nxo_in,nyo_in,ndim_out)
  use params,           only: init_res, ndim, fstrout
  use aotensor_def,     only: init_aotensor
  use integrator,       only: init_integrator
  use tl_ad_tensor,     only: init_tltensor, init_adtensor
  use tl_ad_integrator, only: init_tl_ad_integrator
  implicit none

  integer,intent(in) :: nxa_in, nya_in, nxo_in, nyo_in
  integer,intent(out) :: ndim_out

  PRINT*, 'Model MAOOAM v1.2'
  PRINT*, 'Loading information...'
  call init_res(nxa_in,nya_in,nxo_in,nyo_in)
  call init_aotensor    ! Compute the tensor
  call init_tltensor
  call init_adtensor

  call init_integrator  ! Initialize the integrator
  call init_tl_ad_integrator

  call fstrout(nxa_in,nya_in,nxo_in,nyo_in,cshared)
  dt_day = dt*1.0/(1.032d-04)/(86400.0d0)
  nxa = nxa_in; nya = nya_in
  nxo = nxo_in; nyo = nyo_in

  ndim_out = ndim

endsubroutine

!
! read initial condition
! NOTE: use read_x0 after calling init_maooam
!
subroutine read_x0_maooam(ndim,X0)
  USE params, only: fstrout
  implicit none

  integer,intent(in) ::ndim
  real(kind=8),intent(out) :: X0(0:ndim)
  
  character(80) :: cin
  integer :: n

  cin="x0."//TRIM(cshared)//".bin"
  write(*,*) "read initial conditions from file:", TRIM(cin)
  open(200,file=TRIM(cin),action="read",form="unformatted",access="sequential")
  read(200) n
  if (n/=ndim) then
     write(*,*) "dimension mistach: ndim, ndim_from_file=", ndim, n
     stop 1
  endif
  do n = 0, ndim
     read(200) X0(n)
  enddo; close(200)
  do n = 0, ndim
     write(*,*) n, X0(n)
  enddo
  X0(0) = 1.d0
  !pause "conitue?"
endsubroutine

subroutine write_le_maooam(ctype,ndim,LE) 
  implicit none

  character(*),intent(in) :: ctype
  integer,     intent(in) :: ndim
  real(8),     intent(in) :: LE(ndim)
  
  character(80):: cout
  integer :: n

  cout=TRIM(ctype)//TRIM(cshared)//".txt"
  WRITE(*,*) TRIM(ctype), " written into file:", TRIM(cout)
  OPEN(11,file=TRIM(cout),action="write")
  DO n = 1, ndim
     WRITE(11,*) n, LE(n)
     PRINT*, "LV: l, val=", n, LE(n)
  ENDDO; CLOSE(11)
endsubroutine

!
! get the value of atmos physical vars given the location x', y'
!
subroutine s2g_a(ndim,s,x,y,atyp,aval)
  use params,          only : natm, f0, RR, n0 => n
  use inprod_analytic, only : awavenum
  implicit none

  integer,intent(in)  :: ndim
  real(8),intent(in)  :: s(0:ndim)
  real(8),intent(in)  :: x
  real(8),intent(in)  :: y
  integer,intent(in)  :: atyp
  real(8),intent(out) :: aval

  real(8) :: Fb
  integer :: i, ibs

  aval = 0.d0
  if (atyp==1) then ! psi
     ibs=0    ! 1 -> natm
     do i = 1, natm
        if (awavenum(i)%typ=='A') then
           Fb = sqrt(2.0d0)*cos( awavenum(i)%P * y )
        elseif (awavenum(i)%typ=='K') then
           Fb = 2.0d0*cos( awavenum(i)%M * n0 * x )*sin( awavenum(i)%P * y )
        elseif (awavenum(i)%typ=='L') then
           Fb = 2.0d0*sin( awavenum(i)%H * n0 * x )*sin( awavenum(i)%P * y )
        else
            write(*,*) "[error] s2g_a: unrecognized function type: ", trim(awavenum(i)%typ)
            stop 1
        endif
        aval = aval + s(ibs+i)*Fb
     enddo
  elseif (atyp==2) then         ! theta
     ibs=natm  ! natm+1 -> natm+natm
     do i = 1, natm
        if (awavenum(i)%typ=='A') then
           Fb = sqrt(2.0d0)*cos( awavenum(i)%P * y )
        elseif (awavenum(i)%typ=='K') then
           Fb = 2.0d0*cos( awavenum(i)%M * n0 * x )*sin( awavenum(i)%P * y )
        elseif (awavenum(i)%typ=='L') then
           Fb = 2.0d0*sin( awavenum(i)%H * n0 * x )*sin( awavenum(i)%P * y )
        else
            write(*,*) "[error] s2g_a: unrecognized function type: ", trim(awavenum(i)%typ)
            stop 1
        endif
        aval = aval + s(ibs+i)*Fb
     enddo
     aval = 2*f0*aval/RR
  else
     write(*,*) "[error] s2g_a: unrecognized atyp: " , atyp
     stop 1
  endif

endsubroutine

!
! get the value of ocn physical vars given the location x', y'
!
subroutine s2g_o(ndim,s,x,y,otyp,oval)
  use params,          only : natm, noc, n0 => n
  use inprod_analytic, only : owavenum
  implicit none
 
  integer,intent(in)  :: ndim
  real(8),intent(in)  :: s(0:ndim)
  real(8),intent(in)  :: x
  real(8),intent(in)  :: y
  integer,intent(in)  :: otyp
  real(8),intent(out) :: oval

  real(8) :: Fb, phi
  integer :: j, jbs

  real(8),parameter :: pi = acos(-1.d0)
  
  oval = 0.d0
  if (otyp==1) then ! psi_o
     jbs = 2*natm
     do j = 1, noc
        phi = 2*( (-1.d0)**owavenum(j)%H - 1.d0 )*( (-1.d0)**owavenum(j)%P - 1.d0 )/( owavenum(j)%H * owavenum(j)%P * pi * pi )
        Fb = 2*sin(owavenum(j)%H * n0 * x /2.d0 )*sin(owavenum(j)%P * y )
        !oval = oval + s(jbs+j)*(Fb-phi) 
        oval = oval + s(jbs+j)*Fb ! don't have the 
     enddo
  elseif( otyp==2) then ! deltaT_o
     jbs = 2*natm+noc
     do j = 1, noc
        Fb = 2*sin(owavenum(j)%H * n0 * x /2.d0 )*sin(owavenum(j)%P * y )
        oval = oval + s(jbs+j)*Fb
     enddo
  else
     write(*,*) "[error] s2g_o: unrecognized otyp: " , otyp
     stop 2
  endif

endsubroutine

endmodule m_maooam
