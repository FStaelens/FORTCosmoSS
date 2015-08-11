!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module outputs
  !
  ! includes the array variables containing all the possibly needed output of the main code
  ! any supplementary desired output should be given a name and value and added to the list
  ! of input being allocated during the call to "output_list" (see below).
  ! Output must be of the real type
  !
  ! "output_mat" contains the outputs that are profiles on the numerical grid
  ! "output_cosmo" contains the outputs that are mere real numbers 
  !
  ! for example : if one needs a profile output named "white_rabbit" placed in a new module
  !               named "wonderland" one should give the subroutine "output_list" access to
  !               this module though "use wonderland". The number of elements in output_mat
  !               should then be incremented and two field of the form
  !
  !                      output_mat(n)%name = 'white_rabbit'
  !                      output_mat(n)%val = white_rabbit
  !
  !               where 'n' is the integer index of the new output. Note that the name
  !               of the new element needds to be the name of the variable. This name will
  !               be the one used when listing the names of desired outputs from the
  !               initialisation file (defaut "input.ini")
  !
  implicit none

  type outputtype
     character(len=20) :: name
     real(kind=8), allocatable, dimension(:) :: val
  end type outputtype

  type(outputtype), dimension(44) :: output_mat
  type(outputtype), dimension(35) :: output_cosmo

  logical :: header_exists = .false.

contains
  subroutine output_list
    !
    ! desription : collects all the desired outputs from the accessible modules into
    !              "output_mat" and "output_cosmo" 
    ! 
    use grid
    use BSSN_var
    use matter
    use constraints
    use metric
    use sources
    integer :: j

    do j = 1, size(output_mat)
       if(.not.allocated(output_mat(j)%val)) allocate(output_mat(j)%val(size(x_grid)))
    enddo

    do j = 1, size(output_cosmo)
       if(.not.allocated(output_cosmo(j)%val)) allocate(output_cosmo(j)%val(1))
    enddo
    
    output_mat(1)%name = 'x'
    output_mat(1)%val = x_grid

    output_mat(2)%name = 'rho_m'
    output_mat(2)%val = rho_m
    
    output_mat(3)%name = 'phi'
    output_mat(3)%val = phi

    output_mat(4)%name = 'alpha'
    output_mat(4)%val = alpha

    output_mat(5)%name = 'big_X'
    output_mat(5)%val = big_X
    
    output_mat(6)%name = 'big_chi'
    output_mat(6)%val = big_chi

    output_mat(7)%name = 'psi'
    output_mat(7)%val = psi

    output_mat(8)%name = 'Del_H'
    output_mat(8)%val = Delta_H

    output_mat(9)%name = 'E'
    output_mat(9)%val = E

    output_mat(10)%name = 'del_m'
    output_mat(10)%val = del_m

    output_mat(11)%name = 'K'
    output_mat(11)%val = K

    output_mat(12)%name = 'Sa'
    output_mat(12)%val = Sa

    output_mat(13)%name = 'a'
    output_mat(13)%val = big_a
   
    output_mat(14)%name = 'b'
    output_mat(14)%val = big_b

    output_mat(15)%name = 'Pi_phi'
    output_mat(15)%val = Pi_phi

    output_mat(16)%name = 'Psi_phi'
    output_mat(16)%val = Psi_phi

    output_mat(17)%name = 'Del_M'
    output_mat(17)%val = Delta_M

    output_mat(18)%name = 'E_phi'
    output_mat(18)%val = E_phi

    output_mat(19)%name = 'Sa_phi'
    output_mat(19)%val = Sa_phi
    
    output_mat(20)%name = 'Sb_phi'
    output_mat(20)%val = Sb_phi

    output_mat(21)%name = 'jr_phi'
    output_mat(21)%val = jr_phi

    output_mat(22)%name = 'big_R'
    output_mat(22)%val = big_R

    output_mat(23)%name = 'H_cons_1'
    output_mat(23)%val = H_cons_1

    output_mat(24)%name = 'H_cons_2'
    output_mat(24)%val = H_cons_2
   
    output_mat(25)%name = 'H_cons_3'
    output_mat(25)%val = H_cons_3

    output_mat(26)%name = 'source_m'
    output_mat(26)%val = source_m

    output_mat(27)%name = 'source_phi'
    output_mat(27)%val = source_phi
 
    output_mat(28)%name = 'gam_rr'
    output_mat(28)%val = gam_rr

    output_mat(29)%name = 'gam_thth'
    output_mat(29)%val = gam_thth
    
    output_mat(30)%name = 'drgam_thth'
    output_mat(30)%val = drgam_thth

    output_mat(31)%name = 'big_D'
    output_mat(31)%val = big_D

    output_mat(32)%name = 'v_r'
    output_mat(32)%val = v_r

    output_mat(33)%name = 'Tr(T)'
    output_mat(33)%val = E+Sa+2*Sb

    output_mat(34)%name = 'del_phi'
    output_mat(34)%val = del_phi

    output_mat(35)%name = 'p_anis_phi'
    output_mat(35)%val = Sa_phi-Sb_phi

    output_mat(36)%name = 'Aa'
    output_mat(36)%val = Aa

    output_mat(37)%name = 'x_Mpc'
    output_mat(37)%val = x_grid*l_scale/Mpc2m

    output_mat(38)%name = 'H-H_cosm'
    output_mat(38)%val = -1d0/3d0*K+1d0/3d0*K_cosm

    output_mat(39)%name = 'K-K_cosm'
    output_mat(39)%val = K-K_cosm

    output_mat(40)%name = 'p_phi-p_phi_hom'
    output_mat(40)%val = Sb_phi-p_hom_phi

    output_mat(41)%name = 'del_p_phi'
    output_mat(41)%val = Sb_phi/p_hom_phi-1d0

    output_mat(42)%name = 'phi-phi_hom'
    output_mat(42)%val = phi-phi_hom

    output_mat(43)%name = 'Del_H_rel'
    output_mat(43)%val = Delta_H_rel

    output_mat(44)%name = 'Del_M_rel'
    output_mat(44)%val = Delta_M_rel
    
    ! cosmological output

    output_cosmo(1)%name = 't'
    output_cosmo(1)%val = t

    output_cosmo(2)%name = 'a'
    output_cosmo(2)%val = a
    
    output_cosmo(3)%name = 'adot'
    output_cosmo(3)%val = adot

    output_cosmo(4)%name = 'rho_m'
    output_cosmo(4)%val = rho_hom_m

    output_cosmo(5)%name = 'alpha_cosm'
    output_cosmo(5)%val = alpha_cosm

    output_cosmo(6)%name = 'K_cosm'
    output_cosmo(6)%val = K_cosm

    output_cosmo(7)%name = 'rho'
    output_cosmo(7)%val = rho_cosm

    output_cosmo(8)%name = 'phi_hom'
    output_cosmo(8)%val = phi_hom
    
    output_cosmo(9)%name = 'Pi_phi_hom'
    output_cosmo(9)%val = Pi_phi_hom

    output_cosmo(10)%name = 'alpha_centre'
    output_cosmo(10)%val = alpha(1)

    output_cosmo(11)%name = 'Del_H_cosm'
    output_cosmo(11)%val = Delta_H_cosm

    output_cosmo(12)%name = 'L2_H'
    output_cosmo(12)%val = L2_H

    output_cosmo(13)%name = 'L2_M'
    output_cosmo(13)%val = L2_M

    output_cosmo(14)%name = 'delta_m_c'
    output_cosmo(14)%val = del_m(1)

    output_cosmo(15)%name = 'H'
    output_cosmo(15)%val = adot/a

    output_cosmo(16)%name = 'rdelta_c'
    output_cosmo(16)%val = del_phi(1)/del_m(1)

    output_cosmo(17)%name = 'delta_phi_c'
    output_cosmo(17)%val = del_phi(1)

    output_cosmo(18)%name = 'a_centre'
    output_cosmo(18)%val = a*psi(1)**2/psi_i(1)**2

    output_cosmo(19)%name = 'p'
    output_cosmo(19)%val = p_cosm

    output_cosmo(20)%name = 'rho_phi'
    output_cosmo(20)%val = rho_hom_phi

    output_cosmo(21)%name = 'Om_phi'
    output_cosmo(21)%val = Om_phi

    output_cosmo(22)%name = 'Om_m'
    output_cosmo(22)%val = Om_m

    output_cosmo(23)%name = 'Om_Lamb'
    output_cosmo(23)%val = Om_Lamb

    output_cosmo(24)%name = 'phi_centre'
    output_cosmo(24)%val = phi(1)

    output_cosmo(25)%name = 'Del_H_centre'
    output_cosmo(25)%val = Delta_H(1)

    output_cosmo(26)%name = 't_sync'
    output_cosmo(26)%val = t_sync*t_scale/yr2sec/1e9

    output_cosmo(27)%name = 'del_alpha_c'
    output_cosmo(27)%val = alpha(1)/alpha_cosm-1d0

    output_cosmo(28)%name = 'w_phi_hom'
    output_cosmo(28)%val = p_hom_phi/rho_hom_phi

    output_cosmo(29)%name = 'w_phi_centre'
    output_cosmo(29)%val = 1d0/3d0*(Sa_phi(1)+2d0*Sb_phi(1))/E_phi(1)

    output_cosmo(30)%name = 'rho_phi_centre'
    output_cosmo(30)%val = E_phi(1)

    output_cosmo(31)%name = 'p_phi_centre'
    output_cosmo(31)%val = 1d0/3d0*(Sa_phi(1)+2d0*Sb_phi(1))

    output_cosmo(32)%name = 'p_phi_hom'
    output_cosmo(32)%val = p_hom_phi

    output_cosmo(33)%name = 'Pi_phi_centre'
    output_cosmo(33)%val = Pi_phi(1)

    Output_cosmo(34)%name = 'rho_Lamb'
    output_cosmo(34)%val = rho_hom_Lamb

    Output_cosmo(35)%name = 'a-a_c'
    output_cosmo(35)%val = a-a*psi(1)**2/psi_i(1)**2


  end subroutine output_list

  subroutine append_to_list(list,file_name,header)
    !
    ! input  : list (str) an array containing all the names of the desired output
    !          file_name (str) the name of the output file
    !          header (str) the header for the output file
    ! 
    ! desription : during first call, creates a file ,write the names of the outputs
    !              as column headers under a global header and writes the firs line
    !              of outputs
    !              during later calls, append the new values of the output variables
    !              in their corresponding columns
    !
    !              a warning is issued if one or more output names from "list"
    !              cannot be found within the array variable "output_cosmo(:)%name"
    ! 
    implicit none
    character(len=*), dimension(:), intent(in) :: list
    character(len=*), intent(in) :: file_name, header
    logical :: output_exists = .false.

    integer :: j, k

    open(12,file=file_name,access='append')
    if(header_exists.eqv..false.)then
       write(12,*) header
       do j = 1, size(list)
          output_exists = .false.
          do k = 1, size(output_cosmo)
             if(output_cosmo(k)%name==list(j))then
                if(j==size(list))then
                   write(12,*) ' ', trim(adjustr(list(j)))
                else
                   write(12,'(a,a)',advance='no') ' ', trim(adjustr(list(j)))
                endif
                output_exists = .true.
             end if
          enddo
          if(output_exists.eqv..false.)then
             print*, 'Warning : no matching of "',trim(adjustl(list(j))),'" found in preset output list.'
             print*, '          change in "output.f90".'
             write(12,*) ''
          endif
       enddo
       header_exists = .true.
    endif
    do j = 1, size(list)
       output_exists = .false.
       do k = 1, size(output_cosmo)
          if(output_cosmo(k)%name==list(j))then
             if(j==size(list))then
                write(12,*) ' ', output_cosmo(k)%val(1)
             else
                write(12,'(a,e20.9e3)',advance='no') ' ', output_cosmo(k)%val(1)
             endif
             output_exists = .true.
          endif
       enddo
       if(output_exists.eqv..false.)then
          write(12,*) ''
       endif
    enddo
    close(12)
    
  end subroutine append_to_list

  subroutine print_list(list,file_name,header,warning,mod_output_in)
    !
    ! input  : list (str) an array containing all the names of the desired output
    !          file_name (str) the name of the output file
    !          header (str) the header for the output file
    !          warning (logical, optional) if .true., warnings are enabled
    !          mod_output_in (integer, optional) fixes the number of elements in
    !                                            each output column
    ! 
    ! desription : opens a file and write the desired outputs in columns with
    !              corresponding header labels
    !
    !              a warning is issued if one or more output names from "list"
    !              cannot be found within the array variable "output_mat(:)%name"
    ! 
    use input_param
    implicit none
    character(len=*), dimension(:), intent(in) :: list
    character(len=*), intent(in) :: file_name, header
    logical, optional, intent(in) :: warning
    integer, optional, intent(in) :: mod_output_in
    logical :: output_exists = .false.
    integer :: j, k, l

    integer :: mod_output

    if(present(mod_output_in))then
       mod_output = mod_output_in
    else
       mod_output = 1
    end if
    
    open(11,file=file_name)
    write(11,*) header
    do j = 1, size(list)
       output_exists = .false.
       do k = 1, size(output_mat)
          if(output_mat(k)%name==list(j))then
             if(j==size(list))then
                write(11,*) ' ', trim(adjustr(list(j)))
             else
                write(11,'(a,a)',advance='no') ' ', trim(adjustr(list(j)))
             endif
             output_exists = .true.
          endif
       enddo
       if(output_exists.eqv..false.)then
          if(present(warning))then
             if(warning.eqv..true.)then
                print*, 'Warning : no matching of "',trim(adjustl(list(j))),'" found in preset output list.'
                print*, '          change in "output.f90".'
                write(11,*) ''
             end if
          endif
       endif
    enddo
    
    do l = 1, size(output_mat(1)%val)
       if((mod(l,mod_output).eq.0))then
          do j = 1, size(list)
             output_exists = .false.
             do k = 1, size(output_mat)
                if(output_mat(k)%name==list(j))then
                   if(j==size(list))then
                      write(11,*) ' ', output_mat(k)%val(l)
                   else
                      write(11,'(a,e20.9e3)',advance='no') ' ', output_mat(k)%val(l)
                   endif
                   output_exists = .true.
                endif
             enddo
             if(output_exists.eqv..false.)then
                write(11,*) ''
             endif
          enddo
       endif
    enddo
    close(11)

  end subroutine print_list

  subroutine progress(ratio,percent)
    !
    ! input  : ratio (real) the ratio of the current value of one parameter over
    !                       its final value
    !
    ! output : percent (integer) the percentage corresponding to the value of "ratio"
    !
    ! desription : computes the percentage and print a progress bar in the console
    !              which is printed over the previous one at each call
    ! 
    implicit none
    real(kind=8), intent(in) :: ratio
    integer,intent(out) :: percent
    character(len=48) :: fancy
    character(len=49) :: bar=" \???% |                                        |"
    integer :: j, k

    fancy = "/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\| "
    percent = int(ratio*100)
    if(ratio>1d0) percent=100
    bar(2:2) = fancy(int(percent/2.5)+1:int(percent/2.5)+1)
    write(bar(3:5),fmt="(i3)") percent
    do k = 1, int(percent/2.5)
       bar(8+k:8+k)="="
    enddo
    write(*,'(1a1,a,$)',advance='no') char(13), bar
    if(percent==100)print*,''
  end subroutine progress


end module outputs
