!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module initial_profiles
  !
  ! includes all the mathematical expressions for the initial profiles used during initilisation
  ! it is intended to be used in "init.f90" and any new initial profile should be added here.
  ! This file has access to the input parameters so that the expressions can be parameters/user
  ! dependent. 
  !
  use input_param
  use cosmov
  use profiles
  implicit none
  logical :: homogeneous = .false.
contains
  function rho_mix(x) 
    implicit none
    real(kind=8) :: x, rho_mix 
    select case(shape_matter)
    case('bump')
       if(H0==0d0)then
          rho_mix = bump(x,delta_mi,xmax_delta_mi)
       else
          rho_mix = rho_hom_m*(1d0+bump(x,delta_mi,xmax_delta_mi)) 
       endif
    case('step') 
       if(H0==0d0)then
          rho_mix = logistic(x,delta_mi,xmax_delta_mi,steepness_mi)
       else
          rho_mix = rho_hom_m*(1d0+logistic(x,delta_mi,xmax_delta_mi,steepness_mi)) 
       endif
    case default
       if(H0==0d0)then
          rho_mix = bump(x,delta_mi,xmax_delta_mi)
       else
          rho_mix = rho_hom_m*(1d0+bump(x,delta_mi,xmax_delta_mi)) 
       endif
    end select  
  end function rho_mix
  
  function phi_ix(x)
    implicit none
    real(kind=8) :: x, phi_ix
    select case(shape_phi)
    case('gaussian','gauss','sym_gaussian','sym_gauss')
       if(H0==0d0)then
          phi_ix = sym_gaussian(x,delta_phii,x0_delta_phii,sig_delta_phii)
       else
          phi_ix = phi_hom*(1d0+sym_gaussian(x,delta_phii,x0_delta_phii,sig_delta_phii))
       endif
    case('bump')
       if(H0==0d0)then
          phi_ix = bump(x,delta_phii,xmax_delta_phii)
       else
          phi_ix = phi_hom*(1d0+bump(x,delta_phii,xmax_delta_phii))
       endif
    end select
  end function phi_ix

  function Psi_phiix(x)
    implicit none
    real(kind=8) :: x, Psi_phiix 
    select case(shape_phi)
    case('gaussian','gauss','sym_gaussian','sym_gauss')
       if(H0==0d0)then
          Psi_phiix = dxsym_gaussian(x,delta_phii,x0_delta_phii,sig_delta_phii)
       else
          Psi_phiix = phi_hom*dxsym_gaussian(x,delta_phii,x0_delta_phii,sig_delta_phii)
       end if
    case('bump')
       if(H0==0d0)then
          Psi_phiix = dxbump(x,delta_phii,xmax_delta_phii)
       else
          Psi_phiix = phi_hom*dxbump(x,delta_phii,xmax_delta_phii)
       end if
    end select
  end function Psi_phiix

  function V_phi_ix(x)
    use potential
    implicit none
    real(kind=8) :: x, V_phi_ix
    real(kind=8) :: phii, dV_phiix    
    phii = phi_ix(x)
    call V_phi(phii,V_phi_ix,dV_phiix)   
  end function V_phi_ix

end module initial_profiles
