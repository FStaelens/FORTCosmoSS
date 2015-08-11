!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------

module constants
! this module includes :
! - all the physical constants needed
! - the time, length and mass scales used in the code in the int. sys. of units.  
  
  implicit none
  real(kind=8), parameter :: pi = 3.141592653589793238462643383279
  real(kind=8), parameter :: c = 299792458 !(in m/s)
  real(kind=8), parameter :: G = 6.67e-11 !(in m³/kg/s²)
  real(kind=8), parameter :: H0_exp = 2.31e-18 !(in 1/s)
  real(kind=8), parameter :: M_sun = 1.99e30 !(in kg)
  real(kind=8), parameter :: yr2sec = 3.16e7
  real(kind=8), parameter :: Mpc2m = 3.08e22

  real(kind=8) :: t_scale, l_scale, m_scale

end module constants
