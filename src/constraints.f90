!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
! includes the variables and subroutine used to represent the violation of the BSSN constraint equations
!
module constraints
  use grid
  use BSSN_var
  use matter
  use cosmov
  use boundaries
  use constants
  use sources
  implicit none

  ! Hamiltonian (H) and momentum (M) constraints violation absolute and relative (rel)
  real(kind=8), dimension(-off+1:nx+off) :: Delta_H, Delta_M, Delta_H_rel, Delta_M_rel
  ! components of the H constraint (for debugging purpose)
  real(kind=8), dimension(-off+1:nx+off) :: big_R, H_cons_1, H_cons_2, H_cons_3
  ! L2 norms of the H and M constraints
  real(kind=8) :: L2_H, L2_M
  ! absolute violation of the cosmological Hamiltonian constraint (Friedmann equation)
  real(kind=8) :: Delta_H_cosm
  
  integer :: i

contains
  subroutine H_cons
    !
    ! description : updates the value of the violation of the Hamiltonian constraint and its L2_norm
    !
    implicit none
    call update_derivatives
    call build_matter_sources

    do i = 1,nx+off
       big_R(i) = -big_X(i)**2/big_a(i)/a**2*(0.5d0*d2xa(i)/big_a(i)+d2xb(i)/big_b(i)-big_a(i)*dxDhr(i)-(dxa(i)/big_a(i))**2 & 
            +0.5d0*(dxb(i)/big_b(i))**2+2d0/x_grid(i)/big_b(i)*(3d0-big_a(i)/big_b(i))*dxb(i) & 
            +4d0/x_grid(i)**2*(1d0-big_a(i)/big_b(i)))-8d0/big_a(i)/a**2*(d2xchi_X2(i)+dxchi(i)*dxchi_X2(i)) & 
            +8d0/big_a(i)/a**2*dxchi_X2(i)*(0.5d0*dxa(i)/big_a(i)-dxb(i)/big_b(i)-2d0/x_grid(i))
       Delta_H(i) = big_R(i) + 2d0/3d0*K(i)**2 - (Aa(i)**2+2d0*Ab(i)**2) - 16d0*pi*E(i)
       Delta_H_rel(i) = (big_R(i) + 2d0/3d0*K(i)**2 - (Aa(i)**2+2d0*Ab(i)**2) - 16d0*pi*E(i)) / (16d0*pi*E(i))
       !big_R(i) + 2d0/3d0*K(i)**2 - (Aa(i)**2+2d0*Ab(i)**2) - 16d0*pi*E(i)
    enddo
    H_cons_1 = 2d0/3d0*K**2
    H_cons_2 = -(Aa**2+2d0*Ab**2)
    H_cons_3 = - 16d0*pi*E
    
    call symmetrise(big_R)
    call symmetrise(Delta_H)
    
    call L2_norm(L2_H,Delta_H,1d0,1)
  end subroutine H_cons
  subroutine M_cons
    !
    ! description : updates the value of the violation of the momentum constraint and its L2_norm
    !
    implicit none
    integer :: i

    call update_derivatives
    call build_matter_sources
    
    do i = 1,nx+off
       Delta_M(i) = dxAa(i)-2d0/3d0*dxK(i)+6d0*Aa(i)*dxchi(i)+(Aa(i)-Ab(i))*(2d0/x_grid(i)+dxb(i)/big_b(i))-8d0*pi*jr(i)
       Delta_M_rel(i) = (dxAa(i)-2d0/3d0*dxK(i)+6d0*Aa(i)*dxchi(i)+(Aa(i)-Ab(i))*(2d0/x_grid(i)+dxb(i)/big_b(i))-8d0*pi*jr(i)) &
            / (8d0*pi*jr(i))
    end do

    call symmetrise(Delta_M)  
    call L2_norm(L2_M,Delta_M,0d0,1)

  end subroutine M_cons     
  subroutine L2_norm(L2_func,func,x0,skip)
    !
    ! input  : func (real) array representing the function of which the L2 norm is computed
    !          x0 (real) the lower bound of the integral
    !          skip (real) the integer step for the computation of the integral as a sum of rectangles
    !
    ! output : L2_func (real) the coputed L2_norm
    !
    ! description : computes the L2 norm as the square root of the sum of the square of the function
    !               divided by the number of computational points 
    !
    use grid
    implicit none
    real(kind=8), intent(out) :: L2_func
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: func
    real(kind=8), intent(in) :: x0
    integer, intent(in) :: skip
    integer :: n = 0
    integer :: i    

    L2_func = 0d0

    do i = 1,nx
       if(mod(i,skip).eq.0) then
          if(x_grid(i).gt.x0) then
             L2_func = L2_func + abs(func(i))**2
             n = n+1
          end if
       end if
    end do
    L2_func = sqrt(L2_func/n)
    
  end subroutine L2_norm
  subroutine H_cons_cosm
    !
    ! description : updates the value of the violation of the Friedmann equation
    !
    implicit none
    call update_cosmov
    
    delta_H_cosm = 2d0/3d0*K_cosm**2-16d0*pi*rho_cosm

  end subroutine H_cons_cosm
end module constraints

