!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module evolution
  !
  ! includes the implementation of the PIRK algorithm (see Cerda-Durran & Cordero-Carrion 2012)
  ! as well as the routines for the building of the differential operators.
  ! any change in the evolution of any dynamical variable should be included in this file
  !
  ! for details on the evolution scheme used and the order in which the variables
  ! are evolved, see Rekier et al. 2015
  !
  ! note : the central value of the density profile is corrected at each step using
  !        a mehtod of interpolation
  !
  use constants
  use grid
  use BSSN_var
  use cosmov
  use matter
  use metric
  use gauge_choice
  use boundaries
  use hydro
  use sources
  use derivatives_fcn
  use interpolation
  use input_param
  implicit none
  real(kind=8), dimension(-off+1:nx+off) :: Pi_phi_star
  real(kind=8) :: err_interp
  integer(kind=8), parameter :: npid = 4 ! half number of points interpolated
  integer(kind=8), parameter :: npi = 5 ! half number of points for interpolations
  real(kind=8), dimension(1:2*npi) :: x_interp, y_interp
  integer :: i

contains
  subroutine PIRK
    implicit none
    
    ! cosmological variables
    !-----------------------
    real(kind=8) :: L1_alpha_cosmn, L1_a_cosmn, L2_K_cosmn, L3_K_cosmn
    real(kind=8) :: L1_alpha_cosm1, L1_a_cosm1, L2_K_cosm1, L3_K_cosm1
    real(kind=8) :: L2_K_cosmnp1
    real(kind=8) :: alpha_cosm_n, a_cosm_n, K_cosm_n, alpha_dot_cosm_n

    real(kind=8) :: L1_rho_hom_mn, L1_phi_homn, L1_Pi_phi_homn
    real(kind=8) :: L1_rho_hom_m1, L1_phi_hom1, L1_Pi_phi_hom1
    real(kind=8) :: rho_hom_m_n, phi_hom_n, Pi_phi_hom_n

    real(kind=8) :: L1_t_syncn
    real(kind=8) :: L1_t_sync1
    real(kind=8) :: t_sync_n

    ! space-time variables
    ! --------------------
    real(kind=8), dimension(-off+1:nx+off) :: L1_alphan, L1_an, L1_bn, L1_chin 
    real(kind=8), dimension(-off+1:nx+off) :: L2_Kn, L2_Aan, L2_Dhrn, L3_Kn, L3_Aan, L3_Dhrn
    real(kind=8), dimension(-off+1:nx+off) :: L1_alpha1, L1_a1, L1_b1, L1_chi1
    real(kind=8), dimension(-off+1:nx+off) :: L2_K1, L2_Aa1, L2_Dhr1, L3_K1, L3_Aa1, L3_Dhr1
    real(kind=8), dimension(-off+1:nx+off) :: L2_Knp1, L2_Aanp1, L2_Dhrnp1 
    real(kind=8), dimension(-off+1:nx+off) :: big_an, big_bn, alphan, big_chin, Kn, Aan, Abn, Dhrn 

    ! hydro variables
    !--------------------
    real(kind=8), dimension(-off+1:nx+off) :: UD, USr 
    real(kind=8), dimension(-off+1:nx+off) :: L1_UDn, L1_USrn
    real(kind=8), dimension(-off+1:nx+off) :: L1_UD1, L1_USr1 
    real(kind=8), dimension(-off+1:nx+off) :: UDn, UD1, USrn, USr1    

    ! scalar field variables
    !-----------------------
    real(kind=8), dimension(-off+1:nx+off) :: L1_phin, L1_Psi_phin, L1_Pi_phin
    real(kind=8), dimension(-off+1:nx+off) :: L2_Pi_phin, L3_Pi_phin, L2_Psi_phin, L3_Psi_phin
    real(kind=8), dimension(-off+1:nx+off) :: L1_phi1, L1_Psi_phi1, L1_Pi_phi1
    real(kind=8), dimension(-off+1:nx+off) :: L2_Pi_phi1, L3_Pi_phi1, L2_Psi_phi1, L3_Psi_phi1
    real(kind=8), dimension(-off+1:nx+off) :: L2_Pi_phinp1, L2_Psi_phinp1
    real(kind=8), dimension(-off+1:nx+off) :: phin, Psi_phin, Pi_phin

    Pi_phi_star = Pi_phi

    call update_metric
    UD = sqrtdetgam*big_D
    USr = sqrtdetgam*big_Sr
    big_X = exp(-2d0*big_chi)

    ! save variables for later

    alpha_cosm_n = alpha_cosm; alpha_dot_cosm_n = alpha_dot_cosm; a_cosm_n = a; K_cosm_n = K_cosm
    phi_hom_n = phi_hom; Pi_phi_hom_n = Pi_phi_hom
    rho_hom_m_n = rho_hom_m
    t_sync_n = t_sync

    UDn = UD; USrn = USr
    phin = phi; Psi_phin = Psi_phi; Pi_phin = Pi_phi
 
    big_an = big_a; big_bn = big_b; alphan = alpha; big_chin = big_chi
    Kn = K; Aan = Aa; Abn = Ab; Dhrn = Dhr

    ! Pirk step one
     
    ! - build operators
    call make_L1_cosmo(L1_alpha_cosmn,L1_a_cosmn,L1_rho_hom_mn,L1_phi_homn,L1_Pi_phi_homn,L1_t_syncn)
    call make_L2_cosmo(L2_K_cosmn); call make_L3_cosmo(L3_K_cosmn)
    call update_derivatives
    call make_L1(L1_alphan,L1_an,L1_bn,L1_chin,L1_UDn,L1_USrn,L1_phin,L1_Psi_phi=L1_Psi_phin)
    call make_L2(L2_Kn,L2_Aan,L2_Pi_phi=L2_Pi_phin); call make_L3(L3_Kn,L3_Aan,L3_Pi_phi=L3_Pi_phin)
    call make_L2bar(L2_Dhrn); call make_L3bar(L3_Dhrn)

    ! - explicit evolution
    alpha_cosm = alpha_cosm + dt * L1_alpha_cosmn
    a = a + dt * L1_a_cosmn
    rho_hom_m = rho_hom_m + dt * L1_rho_hom_mn
    phi_hom = phi_hom + dt * L1_phi_homn
    Pi_phi_hom = Pi_phi_hom + dt * L1_Pi_phi_homn
    t_sync = t_sync_n + dt * L1_t_syncn

    UD = UD + dt * L1_UDn ; call symmetrise(UD)
    USr = USr + dt * L1_USrn ; call symmetrise(USr)

    phi = phi + dt * L1_phin ; call symmetrise(phi)
    Psi_phi = Psi_phi + dt * L1_Psi_phin ; call anti_symmetrise(Psi_phi)
    !Pi_phi_star = Pi_phi + dt * L1_Pi_phin ; call symmetrise(Pi_phi_star) 

    alpha = alpha + dt * L1_alphan ; call symmetrise(alpha)
    big_a = big_a + dt * L1_an ; call symmetrise(big_a)
    big_b = big_b + dt * L1_bn ; call symmetrise(big_b)
    big_chi = big_chi + dt * L1_chin ; call symmetrise(big_chi)

    big_X = exp(-2d0*big_chi)

    call make_L2_cosmo(L2_K_cosm1)

    ! - implicit evolution
    K_cosm = K_cosm + dt * ( 0.5d0*L2_K_cosmn + 0.5d0*L2_K_cosm1 + L3_K_cosmn )

    call update_derivatives 
    call make_L2(L2_K1,L2_Aa1,L2_Pi_phi=L2_Pi_phi1)

    !Psi_phi = Psi_phi + dt * ( 0.5d0*L2_Psi_phin + 0.5d0*L2_Psi_phi1 + L3_Psi_phin ); call anti_symmetrise(Psi_phi)
    Pi_phi = Pi_phi + dt * ( 0.5d0*L2_Pi_phin + 0.5d0*L2_Pi_phi1 + L3_Pi_phin ); call symmetrise(Pi_phi)
    K = K + dt * ( 0.5d0*L2_Kn + 0.5d0*L2_K1 + L3_Kn ); call symmetrise(K)
    Aa = Aa + dt * ( 0.5d0*L2_Aan + 0.5d0*L2_Aa1 + L3_Aan ); call symmetrise(Aa)
    Ab = -0.5d0*Aa

    call update_derivatives
    call make_L2bar(L2_Dhr1)!,L2_Psi_phi=L2_Psi_phi1)

    Dhr = Dhr + dt * ( 0.5d0*L2_Dhrn + 0.5d0*L2_Dhr1 + L3_Dhrn); call anti_symmetrise(Dhr)
    !Psi_phi = Psi_phi + dt * ( 0.5d0*L2_Psi_phin + 0.5d0*L2_Psi_phi1 + L3_Psi_phin); call anti_symmetrise(Psi_phi)

    call cons2prim(UD,USr,rho_m,v_r,v_up_r,Wlorentz)
    call build_matter_sources

    call update_metric
    
    !big_D = UD/sqrtdetgam
    !big_Sr = USr/sqrtdetgam

    !print*, UD(nx:nx+off)
    !read(*,*)

    big_D = rho_m * Wlorentz
    big_Sr = rho_m * hentalpy * Wlorentz**2 * v_r 

    ! Pirk step two

    ! - build operators
    call make_L1_cosmo(L1_alpha_cosm1,L1_a_cosm1,L1_rho_hom_m1,L1_phi_hom1,L1_Pi_phi_hom1,L1_t_sync1)
    call make_L3_cosmo(L3_K_cosm1)
    call update_derivatives
    call make_L1(L1_alpha1,L1_a1,L1_b1,L1_chi1,L1_UD1,L1_USr1,L1_phi=L1_phi1,L1_Psi_phi=L1_Psi_phi1)
    call make_L3(L3_K1,L3_Aa1,L3_Pi_phi=L3_Pi_phi1)
    call make_L3bar(L3_Dhr1,L3_Psi_phi=L3_Psi_phi1)

    ! - explicit evolution
    alpha_cosm = 0.5d0 * (alpha_cosm_n + alpha_cosm + dt * L1_alpha_cosm1)
    a = 0.5d0 * (a_cosm_n + a + dt * L1_a_cosm1)
    rho_hom_m = 0.5d0 * (rho_hom_m_n + rho_hom_m + dt * L1_rho_hom_m1)
    phi_hom = 0.5d0 * (phi_hom_n + phi_hom + dt * L1_phi_hom1)
    Pi_phi_hom = 0.5d0 * (Pi_phi_hom_n + Pi_phi_hom + dt * L1_Pi_phi_hom1)
    t_sync = 0.5d0 * (t_sync_n + t_sync + dt * L1_t_sync1)

    UD = 0.5d0 * (UDn + UD + dt * L1_UD1) ; call symmetrise(UD)
    USr = 0.5d0 * (USrn + USr + dt * L1_USr1) ; call anti_symmetrise(USr)

    phi = 0.5d0 * (phin + phi + dt * L1_phi1) ; call symmetrise(phi)
    Psi_phi = 0.5d0 * (Psi_phin + Psi_phi + dt * L1_Psi_phi1) ; call anti_symmetrise(Psi_phi)
    !Pi_phi_star = 0.5d0 * (Pi_phin + Pi_phi_star + dt * L1_Pi_phi1) ; call symmetrise(Pi_phi_star)

    alpha = 0.5d0 * (alphan + alpha + dt * L1_alpha1) ; call symmetrise(alpha)
    big_a = 0.5d0 * (big_an + big_a + dt * L1_a1) ; call symmetrise(big_a)
    big_b = 0.5d0 * (big_bn + big_b + dt * L1_b1) ; call symmetrise(big_b)
    big_chi = 0.5d0 * (big_chin + big_chi + dt * L1_chi1) ; call symmetrise(big_chi)
    big_X = exp(-2d0*big_chi)

    call make_L2_cosmo(L2_K_cosmnp1)

    ! - implicit evolution
    K_cosm = K_cosm_n + 0.5d0 * dt * (L2_K_cosmn + L2_K_cosmnp1 + L3_K_cosmn + L3_K_cosm1)

    call update_derivatives
    call make_L2(L2_Knp1,L2_Aanp1,L2_Pi_Phi=L2_Pi_phinp1)

    Pi_phi = Pi_phin + 0.5d0 * dt * (L2_Pi_phin + L2_Pi_phinp1 + L3_Pi_phin + L3_Pi_phi1); call symmetrise(Pi_phi)

    K = Kn + 0.5d0 * dt * (L2_Kn + L2_Knp1 + L3_Kn + L3_K1); call symmetrise(K)
    Aa = Aan + 0.5d0 * dt * (L2_Aan + L2_Aanp1 + L3_Aan + L3_Aa1); call symmetrise(Aa)
    Ab = -0.5d0*Aa

    call update_derivatives
    call make_L2bar(L2_Dhrnp1,L2_Psi_phinp1)

    Dhr = Dhrn + 0.5d0 * dt * (L2_Dhrn + L2_Dhrnp1 + L3_Dhrn + L3_Dhr1); call anti_symmetrise(Dhr)  
    !Pi_phi = Pi_phin + 0.5d0 * dt * (L2_Pi_phin + L2_Pi_phinp1 + L3_Pi_phin + L3_Pi_phi1); call symmetrise(Pi_phi)
    !Psi_phi = Psi_phin + 0.5d0 * dt * (L2_Psi_phin + L2_Psi_phinp1 + L3_Psi_phin + L3_Psi_phi1); call anti_symmetrise(Psi_phi)

    big_X = exp(-2d0*big_chi)
    psi = 1/sqrt(big_X)

    call update_metric

    call cons2prim(UD,USr,rho_m,v_r,v_up_r,Wlorentz)
    call build_matter_sources
    call update_cosmov

    x_interp(1:npi) = -x_grid(npid+1:npid+npi)
    x_interp(npi+1:2*npi)= x_grid(npid+1:npid+npi) 

    y_interp(1:npi) = rho_m(npid+1:npid+npi)
    y_interp(npi+1:2*npi) = rho_m(npid+1:npid+npi)

    do i=1,npid
       call polint(x_interp,y_interp,int(2*npi),x_grid(i),rho_m(i),err_interp)
    end do
    call symmetrise(rho_m)

    if(rho_hom_m/=0d0)then
       del_m = rho_m/rho_hom_m-1d0
    else
       del_m = 0d0
    endif

    if(rho_hom_phi/=0d0)then
       del_phi = E_phi/rho_hom_phi-1d0
    else
       del_phi = 0d0
    endif

    call update_metric
    
    !big_D = UD/sqrtdetgam
    !big_Sr = USr/sqrtdetgam

    big_D = rho_m * Wlorentz
    big_Sr = rho_m * hentalpy * Wlorentz**2 * v_r 

  end subroutine PIRK

  subroutine make_L1_cosmo(L1_alpha_cosm,L1_a_cosm,L1_rho_hom_m,L1_phi_hom,L1_Pi_phi_hom,L1_tsync)
    use potential
    implicit none 
    real(kind=8), intent(out) :: L1_alpha_cosm, L1_a_cosm, L1_rho_hom_m, L1_phi_hom, L1_Pi_phi_hom, L1_tsync
    real(kind=8) :: V, dV_dphi

    call update_cosmov
    call V_phi(phi_hom,V,dV_dphi)

    L1_alpha_cosm = alpha_dot_cosm
    L1_a_cosm = adot  
    L1_rho_hom_m = alpha_cosm * rho_hom_m * K_cosm
    L1_tsync = alpha_cosm 
    L1_phi_hom = alpha_cosm * Pi_phi_hom
    L1_Pi_phi_hom = alpha_cosm * Pi_phi_hom * K_cosm - alpha_cosm * dV_dphi
  end subroutine make_L1_cosmo

  subroutine make_L2_cosmo(L2_K_cosm)
    implicit none 
    real(kind=8), intent(out) :: L2_K_cosm
    L2_K_cosm = 0d0
  end subroutine make_L2_cosmo

  subroutine make_L3_cosmo(L3_K_cosm)
    implicit none 
    real(kind=8), intent(out) :: L3_K_cosm
    call update_cosmov
    L3_K_cosm = 1d0/3d0*alpha_cosm*K_cosm**2+4d0*pi*alpha_cosm*(rho_cosm+3*p_cosm)
  end subroutine make_L3_cosmo

  subroutine make_L1(L1_alpha,L1_a,L1_b,L1_chi,L1_UD,L1_USr,L1_phi,L1_Psi_phi,L1_Pi_phi)
    use derivatives_fcn
    use potential
    implicit none
    ! hydro temp variables
    !-----------------------
    real(kind=8), dimension(-off+1:nx+off) :: rho_m_left, rho_m_right
    real(kind=8), dimension(-off+1:nx+off) :: v_up_r_left, v_up_r_right
    real(kind=8), dimension(-off+1:nx+off) :: Flux_D, Flux_Sr
    real(kind=8), dimension(-off+1:nx+off) :: UD, USr 

    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L1_alpha, L1_a, L1_b, L1_chi
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L1_UD, L1_USr, L1_phi
    real(kind=8), dimension(-off+1:nx+off), intent(out), optional :: L1_Psi_phi, L1_Pi_phi
    real(kind=8) :: alpha_left, alpha_right

    real(kind=8), dimension(1:nx+off) :: dxPsi, dxalphaPi
    real(kind=8) :: V, dV_dphi

    dxPsi = dxf(Psi_phi)
    dxalphaPi = dxf(alpha*Pi_phi)

    call update_cosmov
    call update_metric

    UD = sqrtdetgam*big_D
    USr = sqrtdetgam*big_Sr

    ! - hydro rhs
    ! - hydro rhs -- fluxes building
    call build_lr(rho_m,rho_m_left,rho_m_right)
    call build_lr(v_up_r,v_up_r_left,v_up_r_right)
    call hlle(rho_m_left,rho_m_right,v_up_r_left,v_up_r_right,Flux_D,Flux_Sr)
    ! - hydro rhs -- source terms
    call build_hydro_sources

    do i=1,nx 
       alpha_right = 0.5d0*(alpha(i)+alpha(i+1))
       alpha_left = 0.5d0*(alpha(i)+alpha(i-1))
       L1_UD(i) = -(alpha_right*Flux_D(i)-alpha_left*Flux_D(i-1))/dx -eps*Delta4_x_i(UD,i)
       L1_USr(i) = -(alpha_right*Flux_Sr(i)-alpha_left*Flux_Sr(i-1))/dx+S(2,i) -eps*Delta4_x_i(USr,i)
    end do
    do i=1,nx+off
       call V_phi(phi(i),V,dV_dphi) 
       L1_alpha(i) = -alpha(i)**2*f_alpha(alpha(i))*K(i) -eps*Delta4_x_i(alpha,i)
       L1_a(i) = -2d0 * alpha(i) * big_a(i) * Aa(i) -eps*Delta4_x_i(big_a,i)
       L1_b(i) = -2d0 * alpha(i) * big_b(i) * Ab(i) -eps*Delta4_x_i(big_b,i)
       L1_chi(i) = -1d0/6d0 * alpha(i) * K(i) - 0.5d0 * adot/a -eps*Delta4_x_i(big_chi,i)
       L1_phi(i) = alpha(i)*Pi_phi(i) -eps*Delta4_x_i(phi,i)
       if(present(L1_Psi_phi))then
          L1_Psi_phi(i) = dxalphaPi(i) -eps*Delta4_x_i(Psi_phi,i)
       end if
       if(present(L1_Pi_phi))then
          L1_Pi_phi(i) = alpha(i)/big_a(i)*big_X(i)**2/a**2*(dxPsi(i)+Psi_phi(i)*(2d0/x_grid(i)-dxa(i)/2d0/big_a(i) &
            +dxb(i)/big_b(i)+2d0*dxchi(i)+dxalpha(i)/alpha(i)))+alpha(i)*K(i)*Pi_phi(i)-alpha(i)*dV_dphi &
            -eps*Delta4_x_i(Pi_phi,i)
       end if
    end do

    call update_cosmov

    call sommerfeld(alpha,L1_alpha,alpha_cosm,alpha_dot_cosm,sqrt(2d0)) 
    call sommerfeld(big_a,L1_a,1d0,0d0,1d0)
    call sommerfeld(big_b,L1_b,1d0,0d0,1d0)
    call sommerfeld(big_chi,L1_chi,0d0,0d0,1d0)
    call sommerfeld(UD,L1_UD,a**3*rho_hom_m,0d0,0d0) ! v_matter=0 asymptotically
    call sommerfeld(USr,L1_USr,0d0,0d0,1d0)
    call sommerfeld(phi,L1_phi,phi_hom,alpha_cosm*Pi_phi_hom,1d0)
    if(present(L1_Psi_phi))then
       call sommerfeld(Psi_phi,L1_Psi_phi,0d0,0d0,1d0)
    end if
    if(present(L1_Pi_phi))then
       call V_phi(phi_hom,V,dV_dphi)
       call sommerfeld(Pi_phi,L1_Pi_phi,Pi_phi_hom,alpha_cosm*Pi_phi_hom*K_cosm-alpha_cosm*dV_dphi,1d0)
    end if

  end subroutine make_L1

  subroutine make_L2(L2_K,L2_Aa,L2_Psi_phi,L2_Pi_phi)
    use derivatives_fcn
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L2_K, L2_Aa
    real(kind=8), dimension(-off+1:nx+off), intent(out), optional :: L2_Psi_phi, L2_Pi_phi
    real(kind=8), dimension(1:nx+off) :: big_R, Rrr, nabl2_alpha, nabrnabr_alpha
    real(kind=8), dimension(1:nx+off) :: dxalphaPi, dxphi, dxPsi

    dxalphaPi = dxf(alpha*Pi_phi_star)
    dxphi = dxf(phi)
    dxPsi = dxf(Psi_phi)

    do i=1,nx ! dissipation at right boundary in L3 op

       Rrr(i) = -big_X(i)**2/big_a(i)/a**2*(0.5d0*d2xa(i)/big_a(i)-0.75d0*(dxa(i)/big_a(i))**2+& 
         0.5d0*(dxb(i)/big_b(i))**2+dxa(i)/(x_grid(i)*big_b(i))+& 
         (2d0*(1d0-big_a(i)/big_b(i))/(x_grid(i)*x_grid(i)))*(1d0+x_grid(i)*dxb(i)/big_b(i)))-(1d0/big_a(i)/a**2)*&
         (4d0*d2xchi_X2(i)-2d0*dxchi_X2(i)*(dxa(i)/(big_a(i))-dxb(i)/big_b(i)-2d0/x_grid(i)))

       big_R(i) = -big_X(i)**2/big_a(i)/a**2*(0.5d0*d2xa(i)/big_a(i)+d2xb(i)/big_b(i)-& 
         (dxa(i)/big_a(i))*(dxa(i)/big_a(i))+0.5d0*(dxb(i)/big_b(i))*& 
         (dxb(i)/big_b(i))+2d0/(x_grid(i)*big_b(i))*(3d0-big_a(i)/big_b(i))*dxb(i)+& 
         4.d0/(x_grid(i)*x_grid(i))*(1d0-big_a(i)/big_b(i)))-(1d0/big_a(i)/a**2)*8d0*(d2xchi_X2(i)+dxchi(i)*dxchi_X2(i))+& 
         (1d0/big_a(i)/a**2)*8d0*dxchi_X2(i)*(dxa(i)/(2d0*big_a(i))-dxb(i)/(big_b(i))-2d0/x_grid(i))


       nabl2_alpha(i) = big_X(i)**2/big_a(i)/a**2*(d2xalpha(i)-dxalpha(i)*(0.5d0*dxa(i)/big_a(i)- & 
         dxb(i)/big_b(i)-2d0*dxchi(i)-2d0/x_grid(i)))

       nabrnabr_alpha(i) = big_X(i)**2/big_a(i)/a**2*(d2xalpha(i)-dxalpha(i)*(0.5d0*dxa(i)/big_a(i))) & 
         -2d0*dxalpha(i)*dxchi_X2(i)/big_a(i)/a**2

       L2_K(i) = - nabl2_alpha(i)

       L2_Aa(i) = -(nabrnabr_alpha(i)-1d0/3d0*nabl2_alpha(i))+alpha(i)*(Rrr(i)-1d0/3d0*big_R(i))

       ! note: the expressions for 'R_rr' and 'big_R' given above do not include terms with Dhr
       !       those are included in the explicit operators (L_3)

       if(present(L2_Psi_phi))then
          L2_Psi_phi(i) = dxalphaPi(i)
       end if
       if(present(L2_Pi_phi))then
          L2_Pi_phi(i) = alpha(i)/big_a(i)*big_X(i)**2/a**2*(dxphi(i)*(2d0/x_grid(i)-0.5d0*dxa(i)/big_a(i)&
               +dxb(i)/big_b(i)+2d0*dxchi(i)+dxalpha(i)/alpha(i)))
          !alpha(i)/big_a(i)*big_X(i)**2/a**2*(dxPsi(i)+dxphi(i)*(2d0/x_grid(i)-0.5d0*dxa(i)/big_a(i)&
          !+dxb(i)/big_b(i)+2d0*dxchi(i)+dxalpha(i)/alpha(i)))
       End if
    enddo
  end subroutine make_L2

  subroutine make_L3(L3_K,L3_Aa,L3_Psi_phi,L3_Pi_phi)
    use derivatives_fcn
    use potential
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L3_K, L3_Aa
    real(kind=8), dimension(-off+1:nx+off), intent(out), optional :: L3_Psi_phi, L3_Pi_phi
    real(kind=8), dimension(1:nx+off) :: dxPsi
    real(kind=8) :: V, dV_dphi

    dxPsi = dxf(Psi_phi)

    call build_matter_sources
    do i=1,nx+off
       call V_phi(phi(i),V,dV_dphi)
       L3_K(i) = alpha(i)*(Aa(i)**2+2d0*Ab(i)**2+1d0/3d0*K(i)**2)+4d0*pi*alpha(i)*(E(i)+Sa(i)+2d0*Sb(i))&
            -eps*Delta4_x_i(K,i)
       L3_Aa(i) = alpha(i)*big_X(i)**2/big_a(i)/a**2*(0.5d0*Dhr(i)*dxa(i)+2d0/3d0*big_a(i)*dxDhr(i))&
            +alpha(i)*K(i)*Aa(i)-16d0/3d0*pi*alpha(i)*(Sa(i)-Sb(i)) -eps*Delta4_x_i(Aa,i)
       if(present(L3_Psi_phi))then
          L3_Psi_phi(i) = 0d0 -eps*Delta4_x_i(Psi_phi,i)
       end if
       if(present(L3_Pi_phi))then
          L3_Pi_phi(i) = alpha(i)*K(i)*Pi_phi(i)+alpha(i)/big_a(i)*big_X(i)**2/a**2*dxPsi(i)-alpha(i)*dV_dphi&
               -eps*Delta4_x_i(Pi_phi,i)
       end if
    end do

    call update_cosmov
    call V_phi(phi_hom,V,dV_dphi)
    call sommerfeld(K,L3_K,K_cosm,1d0/3d0*alpha_cosm*K_cosm**2+4d0*pi*alpha_cosm*(rho_cosm+3*p_cosm),1d0)
    call sommerfeld(Aa,L3_Aa,0d0,0d0,1d0)
    if(present(L3_Psi_phi))then
       call sommerfeld(Psi_phi,L3_Psi_phi,0d0,0d0,1d0)
    end if
    if(present(L3_Pi_phi))then
       call V_phi(phi_hom,V,dV_dphi)
       call sommerfeld(Pi_phi,L3_Pi_phi,Pi_phi_hom,alpha_cosm*Pi_phi_hom*K_cosm-alpha_cosm*dV_dphi,1d0)
    end if
  end subroutine make_L3

  subroutine make_L2bar(L2_Dhr,L2_Psi_phi,L2_Pi_phi)
    use derivatives_fcn
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L2_Dhr
    real(kind=8), dimension(-off+1:nx+off), intent(out), optional :: L2_Psi_phi, L2_Pi_phi
    real(kind=8), dimension(1:nx+off) :: dxPsi, dxalphaPi

    dxPsi = dxf(Psi_phi)
    dxalphaPi = dxf(alpha*Pi_phi)

    do i=1,nx  ! dissipation at right boundary in L3 op
       L2_Dhr(i) = -2d0/big_a(i)*(Aa(i)*dxalpha(i)+alpha(i)*dxAa(i))-4d0*alpha(i)/x_grid(i)/big_b(i)*(Aa(i)-Ab(i)) &
            +2d0*alpha(i)/big_a(i)*(dxAa(i)-2d0/3d0*dxK(i)+6d0*Aa(i)*dxchi(i)+(Aa(i)-Ab(i))*(2d0/x_grid(i)+dxb(i)/big_b(i)))
       if(present(L2_Psi_phi))then
          L2_Psi_phi(i) = dxalphaPi(i)
       end if
       if(present(L2_Pi_phi))then
          L2_Pi_phi(i) = alpha(i)/big_a(i)*big_X(i)**2/a**2*(dxPsi(i)+Psi_phi(i)*(2d0/x_grid(i)-0.5d0*dxa(i)/big_a(i)&
            +dxb(i)/big_b(i)+2d0*dxchi(i)+dxalpha(i)/alpha(i)))
       end if
       !alpha(i)/big_a(i)*big_X(i)**2/a**2*dxPsi(i) &
       !+alpha(i)*big_X(i)**2/a**2/big_a(i)*Psi_phi(i)*(2d0/x_grid(i)-dxa(i)/big_a(i)+dxalpha(i)/alpha(i)) &
       !+2d0*alpha(i)*dxchi_X2(i)/a**2/big_a(i)*Psi_phi(i)
    end do
  end subroutine make_L2bar

  subroutine make_L3bar(L3_Dhr,L3_Psi_phi,L3_Pi_phi)
    use potential
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: L3_Dhr
    real(kind=8), dimension(-off+1:nx+off), intent(out), optional :: L3_Psi_phi, L3_Pi_phi
    real(kind=8) :: V, dV_dphi

    call build_matter_sources
    do i=1,nx+off
       call V_phi(phi(i),V,dV_dphi)
       L3_Dhr(i) = 2d0*alpha(i)*Aa(i)*Dhr(i)-16d0*pi*jr(i)*alpha(i)/big_a(i) -eps*Delta4_x_i(Dhr,i)
       if(present(L3_Psi_phi))then
          L3_Psi_phi(i) = 0d0 -eps*Delta4_x_i(Psi_phi,i)
       end if
       if(present(L3_Pi_phi))then
          L3_Pi_phi(i) = +alpha(i)*Pi_phi(i)*K(i)-alpha(i)*dV_dphi -eps*Delta4_x_i(Pi_phi,i)
       end if
    end do

    call update_cosmov
    call V_phi(phi_hom,V,dV_dphi)
    call sommerfeld(Dhr,L3_Dhr,0d0,0d0,sqrt(2d0))
    if(present(L3_Psi_phi))then
       call sommerfeld(Psi_phi,L3_Psi_phi,0d0,0d0,1d0)
    end if
    if(present(L3_Pi_phi))then
       call V_phi(phi_hom,V,dV_dphi)
       call sommerfeld(Pi_phi,L3_Pi_phi,Pi_phi_hom,alpha_cosm*Pi_phi_hom*K_cosm-alpha_cosm*dV_dphi,1d0)
    end if
  end subroutine make_L3bar

end module evolution
