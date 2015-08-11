!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module grid
  ! this file includes all the parameters of the numerical stagered grid + the CFL factor
  ! any change to this file must be followed by a recompiling of the complete code with all
  ! dependencies using a command such as "make clean all"
  ! 
  ! also includes the time parameter "t" (real)
  !
  ! the parameter "off" is an offset determining the number of ghost points of negative radius
  ! as well as the number of additional points at a very big radius (where dissipative boundary
  ! conditions are applied)
  !
  implicit none
  real(kind=8), parameter :: dx = 0.05, x_end = 100.0, x_start = 0.0, CFL = 0.5, dt = CFL*dx
  integer, parameter :: nx=int((x_end-x_start)/dx), off = 3
  integer :: grid_i
  real(kind=8), dimension(-off+1:nx+off) :: x_grid = (/(grid_i*dx-dx/2d0, grid_i=-off+1, nx+off)/)
  real(kind=8) :: t = 0
end module grid
