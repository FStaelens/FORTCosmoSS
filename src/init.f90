!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module init
  !
  ! includes the instructions for the initialisation of all variables
  !
  ! any reduction of the Hamiltonian constraint to a BVP can be solved
  ! using a call to routines from "BVP_M" and the routines
  ! "fsub" and "bcsub" at the end of this file
  !
  ! A warning is issued if the solution is in doubt
  !
  ! the results of the initialisation are printed on screen
  !
use potential
use constants
implicit none

! common BVP variables
integer, parameter :: node=2,npar=1,leftbc=1,rightbc=node+npar-leftbc

contains
  subroutine init_data(success)
    use grid
    use profiles
    use cosmov
    use matter
    use BSSN_var
    use metric
    use input_param
    use root_finding
    use gauge_choice
    use bvp_m
    use initial_profiles
    implicit none

    logical, intent(inout) :: success
    character(len=1) :: choice

    real(kind=8) :: rho_c_i
    integer :: i, j

    ! BVP variables
    real(kind=8) :: C_psi(npar)
    type(bvp_sol) :: sol
    real(kind=8) :: Singular(node,node), y(node), errest
    real(kind=8) :: y_guess(node) = (/ 1.0d0, 0.0d0 /)
    real(kind=8) :: C_psi_guess(npar) = (/ 0.5d0 /)

    character(len=20) :: col_1, col_2

    print*, ''
    print*, '#==================================================#'
    print*, '" Spherically Symmetric BSSN code by Jeremy Rekier "'
    print*, '" - contact: jrekier@gmail.com                     "'
    print*, '#==================================================#'
    print*, ''
    print*, 'Initialising...'
    print*, '---------------'
    print*, 'Staggered Grid layout :'
    print*, ''
    print'(f11.4,f4.0,f8.4,2f12.4)',-dx/2,0d0,dx/2,3*dx/2,5*dx/2
    print'(a)',' ... ---x-----|-----x-----------x-----------x--- ...'
    print'(a,f4.2,a)','                     < dx=',dx,' >        '
    print*, ''

    call input_list

    print*, 'Output list:'
    print*, '------------'
    write(col_1,*) 'profiles'; write(col_2,*) 'cosmo. var.';
    print*, col_1, col_2
    do j = 1, max(size(output_names),size(output_cosmo_names))
       if(j>size(output_names))then
          write(col_1,*) ''
       else   
          write(col_1,*) '--> ', trim(adjustl(output_names(j)))
       end if
       if(j>size(output_cosmo_names))then
          write(col_2,*) ''
       else
          write(col_2,*) '--> ', trim(adjustl(output_cosmo_names(j)))
       end if
       print*, col_1, col_2
    enddo
    print*, ''

    if (H0==0d0)then
       print'(a)',' Background is Minkowski'
       print*, ''

       ! set the scales of the code based on Schwarzschild units with the sun as mass scale

       t_scale = G*M_sun/c**3 !(in s)
       l_scale = G*M_sun/c**2 !(in m)
       m_scale = M_sun !(in kg)
       
       call init_potential

       print'(a)', ' unit scales for this run :'
       print'(a,es7.1,a)','   t_scale = ', t_scale,' s' 
       print'(a,es7.1,a)','   l_scale = ', l_scale,' m' 
       print'(a,es7.1,a)','   m_scale = ', m_scale,' kg'  
       print*, ''

       phi_hom = phi_V(0d0)

       rho_c_i = 0d0
       rho_hom_Lamb = 0d0; p_hom_Lamb = 0d0
       rho_hom_m = 0d0
       rho_hom_phi = 0d0; p_hom_phi = 0d0

       a = ai 
       alpha_cosm = 1d0
       K_cosm = 0d0
       alpha_dot_cosm = - alpha_cosm**2 * f_alpha(alpha_cosm) * K_cosm
       adot = 0d0

    else
       print'(a)',' Background is Friedmann, '
       print'(a,es7.1)','   initial scale factor ai = ', ai
       print'(a,es7.1)','   initial expansion factor Hi = ', Hi
       print'(a)','   initial energy density parameters:'
       print'(a,f6.4)','        Om_m_i = ', Om_m_i
       print'(a,f6.4)','        Om_phi_i = ', Om_phi_i
       print'(a,f6.4)','        Om_Lamb_i = ', Om_Lamb_i
       print*, ''   
        
       a = ai 
       alpha_cosm = 1d0
       K_cosm = -3d0*Hi
       alpha_dot_cosm = - alpha_cosm**2 * f_alpha(alpha_cosm) * K_cosm
       adot = Hi*ai

       if(Om_lamb_i<0d0)then
          print*, ''
          print*,'Om_lamb_i = ', Om_lamb_i
          print*,'Error : density parameters must sum up to zero'
          print*,''
          stop
       endif

       Om_m = Om_m_i
       Om_phi = Om_phi_i
       Om_lamb = Om_lamb_i


       ! set the scales of the code based on the value of H0

       t_scale = H0/H0_exp !(in s)
       l_scale = c*t_scale !(in m)
       m_scale = c**3/G*t_scale !(in kg)

       print'(a)', ' unit scales for this run (adjust by changing H0) :'
       print'(a,es7.1,a)','   t_scale = ', t_scale,' s' 
       print'(a,es7.1,a)','   l_scale = ', l_scale,' m' 
       print'(a,es7.1,a)','   m_scale = ', m_scale,' kg' 
       print*, ''

       call init_potential

       rho_c_i = 3d0/8d0/pi*Hi**2

       
       rho_hom_Lamb = Om_lamb_i * rho_c_i
       rho_hom_m = Om_m_i * rho_c_i
       rho_hom_phi = Om_phi_i * rho_c_i

       p_hom_lamb = -rho_hom_lamb

       if((ratio_kin_phi>1d0).or.(ratio_kin_phi<0d0))then
          print*, ''
          print*,'ratio_phi = ', ratio_kin_phi
          print*,'Error : the ratio must have value between 0 and 1'
          print*,''
          stop
       endif
       
       p_hom_phi = (2d0*ratio_kin_phi-1d0)*rho_hom_phi

       phi_hom = phi_V((1d0-ratio_kin_phi)*Om_phi_i*rho_c_i)
       Pi_phi_hom = sqrt(2d0*(ratio_kin_phi*Om_phi_i*rho_c_i))

    end if

    rho_cosm = rho_c_i
    p_cosm = p_hom_phi + p_hom_lamb

    ! initialise gauge

    select case(slicing)
    case('geod.','geodesic','geo')   ! geodesic slicing
       write(*,*) ' *Message : using geodesic slicing'
       write(*,*) ''
    case('harm.','harmonic','h')   ! harmonic slicing
       write(*,*) ' *Message : using harmonic slicing'
       write(*,*) ''
    case('1+log','non-advect')   ! 1+log slicing
       write(*,*) ' *Message : using 1+log slicing'
       print*, ''
    case('1thrd','onethird')   ! 1/3 slicing (see Torres et al.) 
       write(*,*) ' *Message : using Bona-Masso slicing with f=1/3'
       write(*,*) ''
    case default    ! default is geodesic slicing
       write(*,*) '*Message : no slicing specified, working in default geodesic slicing '
       write(*,*) ''
    end select


    ! initialise matter

    select case(shape_matter)
    case('bump')
       if(homogeneous.eqv..false.)then
          write(*,*) ' *Message : using bump matter distribution'
          write(*,'(a,es7.1)') '               - central overdensity, del_m_c = ', delta_mi
          write(*,'(a,es7.1)') '               - support size, x_max          = ', xmax_delta_mi
       end if
    case('step')
       if(homogeneous.eqv..false.)then
          write(*,*) ' *Message : using step-like matter distribution'
          write(*,'(a,es7.1)') '               - central overdensity, del_m_c = ', delta_mi
          write(*,'(a,es7.1)') '               - support size, x_max          = ', xmax_delta_mi
          write(*,'(a,es7.1)') '               - steepness, k                 = ', steepness_mi
       end if
    end select
    
    select case(shape_phi)
    case('gaussian','gauss','sym_gaussian','sym_gauss')
       if(homogeneous.eqv..false.)then
          write(*,*) ' *Message : using sym. gaussian phi distribution'
          write(*,'(a,es7.1)') '               - central overdensity, del_phi_c = ', delta_phii
          write(*,'(a,es7.1)') '               - variance, sig_delta            = ', sig_delta_phii
       endif
    case('bump')
       if(homogeneous.eqv..false.)then
          write(*,*) ' *Message : using bump phi distribution'
          write(*,'(a,es7.1)') '               - central overdensity, del_phi_c = ', delta_phii
          write(*,'(a,es7.1)') '               - support size, x_max            = ', xmax_delta_phii
       end if
    end select

    do i = -off+1, nx+off
       alpha(i) = alpha_cosm+sym_gaussian(x_grid(i),delta_alphai,x0_delta_alphai,sig_delta_alphai)
       phi(i) = phi_ix(x_grid(i)) 
       Pi_phi(i) = Pi_phi_hom !0d0
       Psi_phi(i) = Psi_phiix(x_grid(i))
       rho_m(i) = rho_mix(x_grid(i)) 
       rho_Lamb(i) = rho_hom_Lamb
       p_Lamb(i) = p_hom_Lamb
    enddo

    if(rho_hom_m/=0d0)then
       del_m = rho_m/rho_hom_m-1d0
    else
       del_m = 0d0
    endif

    v_up_r = 0d0
    v_r = gam_rr*v_up_r

    Wlorentz = 1d0/sqrt(1d0-v_up_r*v_r)

    spec_energ = 0d0
    hentalpy = 1d0 + spec_energ

    big_D = rho_m*Wlorentz
    big_Sr = rho_m*hentalpy*Wlorentz**2*v_r

    E_m = rho_m*Wlorentz**2; E_Lamb = rho_Lamb;
    Sa_m = rho_m*Wlorentz**2*v_r*v_up_r; Sa_Lamb = p_Lamb;
    Sb_m = 0d0; Sb_Lamb = p_Lamb;
    jr_m = rho_m*Wlorentz**2*v_r; jr_Lamb = 0d0; 

    if((delta_mi==0.and.delta_phii==0)&
         .or.(H0/=0.and.(Om_phi_i*delta_phii==0.and.delta_mi==0))&
         .or.(H0/=0.and.(Om_m_i*delta_mi==0.and.delta_phii==0)))&
    then
       write(*,'(a)') ''
       write(*,'(a)') ' The energy distribution is homogeneous ==> psi = 1'
       write(*,'(a)') ' If this is not what you expected, please check the values of both'
       write(*,'(a)') ' the cosmological density parameters and the density contrasts'
       print*, ''
       psi = 1d0
       homogeneous= .true.
    else
       ! solving boundary value pblm for conformal factor psi^4
    
       write(*,'(a)') ''
       write(*,'(a)') ' Solving Hamiltonian constraint for conformal factor "psi" ...'
       
       sol = bvp_init(node,leftbc,(/ 0d0, x_grid(nx+off) /),y_guess,C_psi_guess)
       Singular = 0d0
       Singular(2,2) = -2d0
       sol = bvp_solver(sol,fsub,bcsub,singularterm=Singular,yerror=errest,tol=1.5d-14,method=6,trace=0)
       
       print*, '...done!'
       
       call bvp_eval(sol,C_psi)
       
       print'(a,es8.1)', ' The adjustement of the right boundary (see Shibata) gives C_psi = ', C_psi
       print*, ''
       
       do i=-off+1,nx+off
          call bvp_eval(sol,x_grid(i),y)
          psi(i)=y(1)
       enddo
       call bvp_terminate(sol)
       
       if(C_psi(1)<0d0-epsilon(0d0))then
          write(*,'(a)') ' Warning : computation of the conformal factor gives the wrong concavity'
          write(*,'(a)') '           if the distribution of matter is inhomogeneous, the result '
          write(*,'(a)') '           has good chances to be invalid.'
          Write(*,'(a)',advance='no') '           Do you wish to proceed further? (y/n) : '
100       read*, choice
          select case(choice)
          case('n')
             success = .false.
             print*,''
          case('y')
             success = .true.
             print*,''
          case default 
             write(*,'(a)',advance='no') '           please choose y or n : '
             goto 100
          end select
       end if     
    end if

    ! initialise BSSN var

    big_a = 1d0; big_b = 1d0
    big_X = 1d0/psi**2
    big_chi = log(psi)
    K = K_cosm
    Aa = 0d0; Ab = 0d0; Dhr = 0d0

    psi_i = psi
    
    ! initialise metric

    gam_rr = a**2/big_X**2*big_a
    gam_thth = a**2/big_X**2*big_b*x_grid**2
    sqrtdetgam = sqrt(gam_rr*gam_thth**2) 

    do i = -off+1, nx+off
       E_phi(i) = 0.5d0*(Pi_phi_hom**2+Psi_phiix(x_grid(i))**2/ai**2/psi(i)**4)+V_phi_ix(x_grid(i))
       Sa_phi(i) = 0.5d0*(Pi_phi_hom**2+Psi_phiix(x_grid(i))**2/ai**2/psi(i)**4)-V_phi_ix(x_grid(i))
       Sb_phi(i) = 0.5d0*(Pi_phi_hom**2-Psi_phiix(x_grid(i))**2/ai**2/psi(i)**4)-V_phi_ix(x_grid(i))
       jr_phi(i) = -Pi_phi_hom*Psi_phiix(x_grid(i))
       source_phi(i) = pi*(Pi_phi_hom**2)+pi*Psi_phiix(x_grid(i))**2
    enddo

    if(rho_hom_phi/=0d0)then
       del_phi = E_phi/rho_hom_phi-1d0
    else
       del_phi = 0d0
    endif

    E = E_m + E_Lamb + E_phi
    Sa = Sa_m + Sa_Lamb + Sa_phi
    Sb = Sb_m + Sb_Lamb + Sb_phi
    jr = jr_m + jr_Lamb + jr_phi

  end subroutine init_data

  subroutine fsub(x,y,p,f)
    use input_param
    use initial_profiles
    implicit none
    real(kind=8) :: x, y(node), f(node), p(npar), C_psi, rhs
    rhs = (-2*pi*(rho_mix(x)+rho_hom_Lamb+V_phi_ix(x)+0.5d0*Pi_phi_hom**2)+0.75d0*Hi**2)*ai**2*y(1)**5 &
          -(pi*Psi_phiix(x)**2)*y(1)
         ! -(pi*Psi_phiix(x)**2)*ai*y(1)
    f = (/ y(2), rhs /)
  end subroutine fsub

  subroutine bcsub(ya,yb,p,bca,bcb)
    use grid
    implicit none
    real(kind=8) :: ya(node), yb(node), bca(leftbc), bcb(rightbc), p(npar), C_psi
   
    C_psi = p(1)
    
    bca(1) = ya(2)
    bcb(1) = yb(1) -1d0 - C_psi/2d0/x_grid(nx+off)
    bcb(2) = yb(2) + C_psi/2d0/(x_grid(nx+off))**2  
  end subroutine bcsub

end module init

