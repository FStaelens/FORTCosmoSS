!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------! This file was written by J.Rekier
! Last edited : 2015-08-10
!
module hydro
  !
  ! includes the procedures needed for the num. integration of hydrodynamical equations
  !
use grid
implicit none
contains
  subroutine build_lr(var,var_left,var_right)
    !
    ! input  : var (real) array containing the variable to be processed
    ! output : var_left (real) array containing an evaluation of "var" on the left side of the cells
    !          var_right (real) array containing an evaluation of "var" on the right side of the cells
    ! 
    ! desription : computes the values of a variable on the right and left sides of each grid cell
    !              using a slope limiter algorithm
    ! 
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: var
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: var_left, var_right
    real(kind=8) :: slope, ss, ssp, c
    integer :: i

    var_left = 0d0
    var_right = 0d0

    ss = 0d0
    ssp = 0d0
    slope = 0d0
    c = 0d0

    do i = 0, nx
       ss = (var(i)-var(i-1))/dx
       ssp = (var(i+1)-var(i))/dx
       c = 0.5d0*(ss+ssp)

       if((ss*ssp).gt.0d0) then
          slope = min(2d0*abs(ss),2d0*abs(ssp),abs(c))*ss/abs(ss)
       else
          slope = 0d0
       end if
       ! var_left := left value of var at i+1/2
       var_left(i) = var(i) + slope*0.5d0*dx
    enddo

    ss = 0d0
    ssp = 0d0
    slope = 0d0
    c = 0d0

    do i = 0, nx
       ss = (var(i+1)-var(i))/dx
       ssp = (var(i+2)-var(i+1))/dx
       c = 0.5d0*(ss+ssp)

       if((ss*ssp).gt.0d0)then
          slope = min(2d0*abs(ss),2d0*abs(ssp),abs(c))*ss/abs(ss)
       else
          slope = 0d0
       end if
       ! var_right := right value of var at i+1/2
       var_right(i) = var(i+1) - slope*0.5d0*dx
    enddo
    
  end subroutine build_lr

  subroutine hlle(rho_left,rho_right,v_up_r_left,v_up_r_right,Flux_D,Flux_Sr)
    !
    ! input  : v_up_r_left (real) array of left states of the fluid radial velocity (index up)
    !          v_up_r_right (real) array of right states of the fluid radial velocity (index up)
    !          rho_left (real) array of left states of the fluid intrinsic energy density
    !          rho_right (real) array right states of the fluid intrinsic energy density
    !
    ! output : Flux_D (real) array of the conservative energy density flux profile
    !          Flux_Sr (real) array of the conservative momentum density flux profile
    ! 
    ! desription : computes the fluxes of the energy and momentum densities using an HLLE method
    !              (see L. Rezzolla and O. Zanotti, Oxford University Press, 2013)
    !
    use metric
    use BSSN_var
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: v_up_r_left, v_up_r_right
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: rho_left, rho_right
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: Flux_D, Flux_Sr
    
    real(kind=8), dimension(-off+1:nx+off) :: gam_up_rr, gam_up_thth
    real(kind=8) :: gam_rr_mean, gam_thth_mean
    real(kind=8) :: gam_up_rr_mean, gam_up_thth_mean
    real(kind=8) :: sqrtdetgam_mean
    real(kind=8) :: alpha_mean, big_X_mean

    real(kind=8) :: Wlorentz_left, Wlorentz_right
    real(kind=8) :: p_left, p_right
    real(kind=8) :: hentalpy_left, hentalpy_right
    real(kind=8) :: UD_left, UD_right, diff_D
    real(kind=8) :: USr_left, USr_right, diff_Sr

    real(kind=8) :: l_minus_left, l_0_left, l_plus_left  ! eigenvalues left
    real(kind=8) :: l_minus_right, l_0_right, l_plus_right  ! eigenvalues right
    real(kind=8) :: l_min, l_max, l_minmax

    real(kind=8) :: Flux_D_left, Flux_D_right
    real(kind=8) :: Flux_Sr_left, Flux_Sr_right

    integer :: i


    gam_up_rr = 1d0/gam_rr
    gam_up_thth = 1d0/gam_thth

    do i = 0, nx
       gam_up_rr_mean = 0.5d0*(gam_up_rr(i)+gam_up_rr(i+1))
       gam_up_thth_mean = 0.5d0*(gam_up_thth(i)+gam_up_thth(i+1))

       gam_rr_mean = 0.5d0*(gam_rr(i)+gam_rr(i+1))
       gam_thth_mean = 0.5d0*(gam_thth(i)+gam_thth(i+1))

       sqrtdetgam_mean = sqrt(gam_rr_mean*gam_thth_mean**2)
       
       alpha_mean = 0.5d0*(alpha(i)+alpha(i+1))
       big_X_mean = 0.5d0*(big_X(i)+big_X(i+1))

       ! left states

       Wlorentz_left = 1d0/sqrt(1d0-gam_rr_mean*v_up_r_left(i)**2)
       p_left = 0d0
       hentalpy_left = 1d0 

       ! ! primitive to conservative 

       UD_left = sqrtdetgam_mean*rho_left(i)*Wlorentz_left
       USr_left = sqrtdetgam_mean*rho_left(i)*Wlorentz_left**2*hentalpy_left*gam_rr_mean*v_up_r_left(i)

       ! ! compute eigenvalues

       call eigenvalue(rho_left(i),v_up_r_left(i),gam_rr_mean,gam_up_rr_mean,alpha_mean,l_minus_left,l_0_left,l_plus_left)

       ! right states

       Wlorentz_right = 1d0/sqrt(1d0-gam_rr_mean*v_up_r_right(i)**2)
       p_right = 0d0
       hentalpy_right = 1d0 

       ! ! primitive to conservative 

       UD_right = sqrtdetgam_mean*rho_right(i)*Wlorentz_right
       USr_right = sqrtdetgam_mean*rho_right(i)*Wlorentz_right**2*hentalpy_right*gam_rr_mean*v_up_r_right(i)

       ! ! compute eigenvalues

       call eigenvalue(rho_right(i),v_up_r_right(i),gam_rr_mean,gam_up_rr_mean,alpha_mean,l_minus_right,l_0_right,l_plus_right)

       ! extremal eigenvalues

       l_min = min(0d0,l_minus_left,l_0_left,l_plus_left,l_minus_right,l_0_right,l_plus_right)
       l_max = max(0d0,l_minus_left,l_0_left,l_plus_left,l_minus_right,l_0_right,l_plus_right)

       l_minmax = l_max-l_min

       ! Friedrich's Variables (state_right-state_left)

       diff_D = UD_right-UD_left
       diff_Sr = USr_right-USr_left

       Flux_D_left = UD_left * v_up_r_left(i)
       Flux_D_right = UD_right * v_up_r_right(i)

       Flux_Sr_left = USr_left * v_up_r_left(i)     ! dust only (p=0) and beta = 0
       Flux_Sr_right = USr_right * v_up_r_right(i)  ! dust only (p=0) and beta = 0

       if(l_minmax.ne.0d0) then

       Flux_D(i) = (l_max*Flux_D_left - l_min*Flux_D_right + l_max*l_min*diff_D)/l_minmax
       Flux_Sr(i) = (l_max*Flux_Sr_left - l_min*Flux_Sr_right + l_max*l_min*diff_Sr)/l_minmax
       
       else

       Flux_D(i) = 0d0
       Flux_Sr(i) = 0d0

       endif
    end do
    

  end subroutine hlle

  subroutine eigenvalue(rho,v_up_r,gam_up_rr,gam_rr,alpha,l_minus,l_0,l_plus)
    !
    ! input  : rho (real) array of the fluid's intrinsinc density profile
    !          v_up_r (real) array of the fluid's radial velocity profile
    !          gam_rr (real) array of radial metric component (indices down)
    !          game_up_rr (real) array of radial metric component (indicies up)
    !          alpha (real) array of lapse profile
    !
    ! output : l_minus (real) eigenvalue for the left travelling mode
    !          l_0 (real) eignevalue giving the propagation speed of matter
    !          l_plus (real) eigenvalue for the right travelling mode
    ! 
    ! desription : computes the eigenvalues of the fluid relativistic Euler equation 
    !              (see L. Rezzolla and O. Zanotti, Oxford University Press, 2013)
    ! 
    real(kind=8), intent(in) :: rho, v_up_r,gam_rr,gam_up_rr,alpha
    real(kind=8), intent(out) :: l_minus, l_0, l_plus

    real(kind=8) :: v2
    real(kind=8) :: cs

    cs = 0d0
    v2 = gam_rr*v_up_r**2

    l_minus = alpha/(1d0-v2*cs**2)*(v_up_r*(1d0-cs**2)-cs*sqrt((1-v2)*(gam_up_rr*(1-v2*cs**2)-v_up_r**2*(1d0-cs**2))))
    l_0 = alpha*v_up_r
    l_plus = alpha/(1d0-v2*cs**2)*(v_up_r*(1d0-cs**2)+cs*sqrt((1-v2)*(gam_up_rr*(1-v2*cs**2)-v_up_r**2*(1d0-cs**2))))

  end subroutine eigenvalue

  subroutine cons2prim(UD,USr,rho,v_r,v_up_r,Wlorentz)
    !
    ! input  : UD (real) array containing the conservative energy density profile
    !          USr (real) array containing the conservative momentum density profile
    !
    ! output : rho (real) array containing the primitive energy density profile
    !          v_r (real) array containting the primitive fluid velocity (index down)
    !          v_up_r (real) array containting the primitive fluid velocity (index up)
    !          Wlorentz (real) array containing the lorentz factor 
    ! 
    ! desription : computes the primitive density and velocity of the fluid from their conservative form
    !              (see L. Rezzolla and O. Zanotti, Oxford University Press, 2013)
    ! 
    use metric
    use BSSN_var
    use boundaries
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: UD, USr
    real(kind=8), dimension(-off+1:nx+off), intent(out) :: rho, v_r, v_up_r, Wlorentz
    real(kind=8), dimension(-off+1:nx+off) :: big_D, Sr 
    real(kind=8) :: Wlorentztilde

    integer :: i


    call update_metric

    big_D = UD/sqrtdetgam
    Sr = USr/sqrtdetgam

    do i=1,nx+off
       if(big_D(i).eq.0d0) then
          Wlorentztilde = 1d0
          rho(i) = 0d0
          v_r(i) = 0d0
       else
          Wlorentztilde = sqrt(1d0+Sr(i)**2/gam_rr(i)/big_D(i)**2)
          rho(i) = big_D(i)/Wlorentztilde
          v_r(i) = Sr(i)/rho(i)/Wlorentztilde**2
       end if
    end do
    
    call symmetrise(rho)
    call anti_symmetrise(v_r)

    v_up_r = 1d0/gam_rr*v_r
    Wlorentz = 1d0/sqrt(1d0-v_up_r*v_r) 

  end subroutine cons2prim

end module hydro
