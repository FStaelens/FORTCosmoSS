!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module potential
  ! includes the routines for the initialisation and the evaluation of a scalar potential function
  ! any new routine using the function parser can be adapted from the structure of the present file
  !
  use constants
  use input_param
  use fparser, only: initf, parsef, evalf, EvalErrType, EvalErrMsg
  implicit none
  real(kind=8), dimension(:), allocatable :: val

contains
  subroutine init_potential
    ! 
    ! desription : process analytical expressions for the potential and its paremeters read
    !              in module "input_param" using module "fparsef"
    !              assign the values of the parameters to array "val" (real) from val(3) onward
    !              val(1) and val(2) are saved to welcome V and phi respectively
    !
    !              once the call to the parser has been issued, the analytical expression for any
    !              parameter or function can be evaluated via a call of the form
    !
    !                        "val(j) = evalf(j, val)"  with j being the index of the function
    !                                 
    !
    integer :: nfunc = 3
    integer :: nvar
    real(kind=8) :: res
    integer :: j
    
    nvar = size(potential_phi%var)

    allocate(val(nvar))
    val = 0d0
    
    call initf(size(potential_phi%param)) 
    do j = 1, size(potential_phi%param)
       call parsef(j,potential_phi%init_param_func(j),potential_phi%var)
    enddo

    do j = 1, size(potential_phi%var)
       if(trim(adjustl(potential_phi%var(j)))=='H0') val(j) = H0
       if(trim(adjustl(potential_phi%var(j)))=='Hi') val(j) = Hi
       if(trim(adjustl(potential_phi%var(j)))=='pi') val(j) = pi
    enddo

    print'(a,a,a)', " scalar potential for this run (from file '",trim(adjustl(potential_filename)),"') :"
    print'(a,a)', '    V(phi) = ', potential_phi%func
    do j = 1, size(potential_phi%param)
       res = evalf(j, val)
       if (EvalErrType > 0)then
          print*,'*** Error: ',EvalErrMsg ()
          stop
       end if
       print'(a,a7,a,es9.3,a,a,a)',&
            '         -- ',trim(adjustl(potential_phi%param(j))),&
            ' = ', res,&
            ' = ( ', trim(adjustl(potential_phi%init_param_func(j))),' )'
       val(2+j) = res
    end do
    print'(a)', ''
    
    call initf(nfunc)
    call parsef(1,potential_phi%func,potential_phi%var)
    call parsef(2,potential_phi%der,potential_phi%var)
    call parsef(3,potential_phi%inv,potential_phi%var)
    
  end subroutine init_potential

  subroutine V_phi(phi,V,dV)
    !
    ! input  : phi (real) the value of the scalar field
    ! output : V (real) the potential evaluated at phi
    !          dV (real) the first derivative of the potential evaluated at phi
    !
    real(kind=8), intent(in) :: phi
    real(kind=8), intent(out) :: V, dV

    val(2) = phi
    val(1) = evalf(1, val)
    V = val(1)
    dV = evalf(2, val)
    
  end subroutine V_phi

  function phi_V(V)
    !
    ! input  : V (real) the value of the potential
    ! output : phi (real) the value of the scalar field obtained after inverting the expression for V
    !
    real(kind=8) :: V
    real(kind=8) :: phi_V
    
    val(1) = V
    val(2) = evalf(3, val)
    phi_V = val(2)
    
  end function phi_V
  
end module potential
