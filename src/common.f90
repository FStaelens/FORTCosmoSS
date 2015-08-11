!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
! includes modules that contain all the variables and routines
! shared throughout the code 
!
module gauge_choice
  !
  ! includes the expression for the function defining the Bona-Masso slicing
  !
  ! note : the Bona-Masso slicing is currently the only one supported by the code
  !        this is a good place to start implementing new options
  !        in shuch event, pay attention to change the evolution for the lapse in
  !        file "evolution.f90".
  !
  use input_param
  implicit none
contains
  function f_alpha(alpha)
    real(kind=8) :: alpha, f_alpha

    ! choice of the f(alpha) of Bona-Masso slicing. 

    select case(slicing)
    case('geod.','geodesic','geo')   ! geodesic slicing
       f_alpha = 0d0
    case('harm.','harmonic','h')   ! harmonic slicing
       f_alpha = 1d0
    case('1+log','non-advect')   ! 1+log slicing
       f_alpha = 2d0/alpha 
    case('1thrd','onethird')   ! 1/3 slicing (see Torres et al.) 
       f_alpha = 0.333d0  
    case default    ! default is geodesic slicing
       f_alpha = 0d0  
    end select
  end function f_alpha
end module gauge_choice

module matter
  !
  ! includes the array variables containing all the necessary energy profiles
  !
  use grid
  implicit none

  ! matter source of initial Poisson equation
  real(kind=8), dimension(-off+1:nx+off) :: source_m, source_phi

  ! density and pressure profiles
  real(kind=8), dimension(-off+1:nx+off) :: rho_m
  real(kind=8), dimension(-off+1:nx+off) :: rho_Lamb, p_Lamb
  
  ! density contrasts
  real(kind=8), dimension(-off+1:nx+off) :: del_m
  real(kind=8), dimension(-off+1:nx+off) :: del_phi
    
  ! other dust hydrodynamical primitive variables
  real(kind=8), dimension(-off+1:nx+off) :: v_up_r, v_r, hentalpy, spec_energ, Wlorentz
  
  ! dust conservative variables
  real(kind=8), dimension(-off+1:nx+off) :: big_D, big_Sr
  
  ! scalar field variables
  real(kind=8), dimension(-off+1:nx+off) :: phi, Pi_phi, Psi_phi
  
  ! Einstein eq. source functions
  real(kind=8), dimension(-off+1:nx+off) :: E, Sa, Sb, jr
  real(kind=8), dimension(-off+1:nx+off) :: E_m, Sa_m, Sb_m, jr_m
  real(kind=8), dimension(-off+1:nx+off) :: E_Lamb, Sa_Lamb, Sb_Lamb, jr_Lamb
  real(kind=8), dimension(-off+1:nx+off) :: E_phi, Sa_phi, Sb_phi, jr_phi
  
  ! Hydro sources functions
  real(kind=8), dimension(3,-off+1:nx+off) :: S
end module matter

module BSSN_var
  !
  ! includes the array variables containing all the necessary BSSN var and their spatial derivatives
  !
  ! subroutine update_derivatives does just what the name implies and compute the derivatives
  !                               of the BSSN variables from the updated values of these
  !
  use grid
  implicit none
  real(kind=8), dimension(-off+1:nx+off) :: alpha  ! lapse
  real(kind=8), dimension(-off+1:nx+off) :: big_a, big_b  ! spatial metric
  real(kind=8), dimension(-off+1:nx+off) :: big_chi, big_X, psi  ! conformal factor
  real(kind=8), dimension(-off+1:nx+off) :: K, Aa, Ab  ! extrinsic curvature
  real(kind=8), dimension(-off+1:nx+off) :: Dhr  ! Delta^r aux. variable

  real(kind=8), dimension(-off+1:nx+off) :: psi_i ! initial value of psi
  
  real(kind=8), dimension(1:nx+off) :: dxalpha, d2xalpha
  real(kind=8), dimension(1:nx+off) :: dxa, d2xa
  real(kind=8), dimension(1:nx+off) :: dxb, d2xb
  real(kind=8), dimension(1:nx+off) :: dxX, d2xX
  real(kind=8), dimension(1:nx+off) :: dxchi, d2xchi, dxchi_X2, d2xchi_X2
  real(kind=8), dimension(1:nx+off) :: dxDhr
  real(kind=8), dimension(1:nx+off) :: dxAa, dxAb
  real(kind=8), dimension(1:nx+off) :: dxK
  
contains
  subroutine update_derivatives
    use derivatives_fcn
    implicit none
    dxalpha = dxf(alpha)
    d2xalpha = d2xf(alpha)
    dxa = dxf(big_a)
    d2xa = d2xf(big_a)
    dxb = dxf(big_b)
    d2xb = d2xf(big_b)
    dxX = dxf(big_X)
    d2xX = d2xf(big_X)
    dxDhr = dxf(Dhr)
    dxK = dxf(K)
    dxAa = dxf(Aa)
    dxAb = dxf(Ab)

    dxchi = -0.5d0*dxX/big_X(1:nx+off)
    d2xchi = -0.5d0*(d2xX/big_X(1:nx+off)-dxX**2/big_X(1:nx+off)**2) 

    dxchi_X2 = -0.5d0*dxX*big_X(1:nx+off)
    d2xchi_X2 = -0.5d0*(d2xX*big_X(1:nx+off)-dxX**2)
    
  end subroutine update_derivatives
end module BSSN_var

module cosmov
  !
  ! includes the cosmological variables
  !
  ! subroutine update_cosmov does just what the name implies and compute the updated
  !                          values of the cosmological variables derived from others
  !
  implicit none
  ! space-time functions
  real(kind=8) :: a, K_cosm, alpha_cosm , alpha_dot_cosm 
  real(kind=8) :: t_sync =0  ! synchronous cosmological time
  real(kind=8) :: adot
  
  ! scalar field functions
  real(kind=8) :: phi_hom, Pi_phi_hom
  
  ! energy functions
  real(kind=8) :: rho_hom_m
  real(kind=8) :: rho_hom_Lamb, p_hom_Lamb
  real(kind=8) :: rho_hom_phi, p_hom_phi
  real(kind=8) :: rho_cosm, p_cosm

  real(kind=8) :: Om_m, Om_phi, Om_Lamb

contains 
  subroutine update_cosmov
    use potential
    use gauge_choice
    real(kind=8) :: V, dV
    call V_phi(phi_hom,V,dV)
    rho_hom_phi = 0.5d0*Pi_phi_hom**2+V
    p_hom_phi = 0.5d0*Pi_phi_hom**2-V
    rho_cosm = rho_hom_m + rho_hom_Lamb + rho_hom_phi
    p_cosm = p_hom_Lamb + p_hom_phi
    alpha_dot_cosm = -alpha_cosm**2*f_alpha(alpha_cosm)*K_cosm
    adot = - alpha_cosm/3d0 * K_cosm * a
    Om_m = rho_hom_m/(rho_hom_m+rho_hom_phi)
    Om_phi = rho_hom_phi/(rho_hom_m+rho_hom_phi)
    Om_Lamb = rho_hom_Lamb/(rho_hom_Lamb+rho_hom_Lamb)
  end subroutine update_cosmov
end module cosmov

module metric
  !
  ! includes the metric components and their spatial derivatives
  !
  ! subroutine update_metric does just what the name implies and compute the updated
  !                          values from the updated values of the BSSN and cosmo var.
  !
  use grid
implicit none
real(kind=8), dimension(-off+1:nx+off) :: gam_rr, gam_thth, sqrtdetgam
real(kind=8), dimension(-off+1:nx+off) :: drgam_thth
contains
  subroutine update_metric
    use BSSN_var
    use cosmov
    use derivatives_fcn
    use boundaries
    implicit none
    gam_rr = a**2/big_X**2*big_a
    gam_thth = a**2/big_X**2*big_b*x_grid**2
    sqrtdetgam = sqrt(gam_rr*gam_thth**2)

    drgam_thth(1:nx+off) = dxf(gam_thth)
    call anti_symmetrise(drgam_thth)
  end subroutine update_metric
end module metric

