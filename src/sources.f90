!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
!
! includes subroutine "build_hydro_sources" used to build the source terms of the hydro equations
!        & subroutine "build_matter_sources" used to combine all the energy components of space-time
!                         into the BSSN source terms E, Sa, Sb and jr (see Rekier et al. 2015)
!
module sources
  use matter
  implicit none
contains
  subroutine build_hydro_sources
    use BSSN_var
    use derivatives_fcn
    use metric
    use boundaries
    implicit none

    real(kind=8), dimension(-off+1:nx+off) :: T_up_00, T_up_0r, T_up_rr, T_up_0_r, T_up_r_r
    real(kind=8), dimension(-off+1:nx+off) :: drgam_rr, dralpha
    real(kind=8), dimension(-off+1:nx+off) :: K_up_r_r

    T_up_00 = rho_m*hentalpy*Wlorentz**2/alpha**2 
    T_up_rr = rho_m*hentalpy*Wlorentz**2*v_up_r**2
    T_up_r_r = rho_m*hentalpy*Wlorentz**2*v_up_r*v_r
    T_up_0r = rho_m*hentalpy*Wlorentz**2*v_up_r
    T_up_0_r = rho_m*hentalpy*Wlorentz**2*v_r

    K_up_r_r = Aa + 1d0/3d0 * K

    drgam_rr(1:nx+off) = dxf(gam_rr)
    call anti_symmetrise(drgam_rr)
    dralpha(1:nx+off) = dxf(alpha)
    call anti_symmetrise(dralpha)

    S(1,:) = 0d0
    S(2,:) = alpha*sqrtdetgam*(-T_up_00*(alpha*dralpha)+0.5d0*T_up_rr*(drgam_rr))
    S(3,:) = alpha*sqrtdetgam*(-T_up_0r*dralpha+T_up_r_r*K_up_r_r)

  end subroutine build_hydro_sources
  subroutine build_matter_sources
    use potential
    use BSSN_var
    use cosmov
    implicit none
    real(kind=8), dimension(-off+1:nx+off) :: V, dV
    integer :: i


    E_m = rho_m * hentalpy * Wlorentz**2 
    Sa_m = rho_m * hentalpy * Wlorentz**2 * v_r * v_up_r 
    Sb_m = 0d0
    jr_m = rho_m * hentalpy * Wlorentz**2 * v_r

    do i = -off+1, nx+off
       call V_phi(phi(i),V(i),dV(i))
    enddo

    E_phi = 0.5d0*(Pi_phi**2+Psi_phi**2*big_X**2/a**2/big_a)+V
    Sa_phi = 0.5d0*(Pi_phi**2+Psi_phi**2*big_X**2/a**2/big_a)-V
    Sb_phi = 0.5d0*(Pi_phi**2-Psi_phi**2*big_X**2/a**2/big_a)-V
    jr_phi = -Pi_phi*Psi_phi

    E = E_m + E_Lamb + E_phi 
    Sa = Sa_m + Sa_Lamb + Sa_phi
    Sb = Sb_m + Sb_Lamb + Sb_phi
    jr = jr_m + jr_Lamb + jr_phi

  end subroutine build_matter_sources
end module sources
