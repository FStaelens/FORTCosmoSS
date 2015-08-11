!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module input_param
  !
  ! this file includes all the input variables required in the main code
  ! all supplementary variable should be added here and assigned using routines from file "inifile.f90"
  !
  ! the names and values of each input variables can be found in a file placed in parent directory
  ! default name is "input.ini" (this can be changed below in first call to "read_input")
  !
  ! an analytic expression for an algebraic definition (e.g. a potential function) can be processed
  ! from a seperate file by adapting the routine "input_potential_list" (see below). The name of
  ! the separate file can be read along with the other variables during the call to "input_list"
  !
  
  implicit none
  real(kind=8) :: delta_mi, xmax_delta_mi, steepness_mi 
  real(kind=8) :: delta_phii, x0_delta_phii, sig_delta_phii, xmax_delta_phii
  real(kind=8) :: delta_alphai, x0_delta_alphai , sig_delta_alphai

  real(kind=8) :: H0, ai, Hi, Om_m_i, Om_phi_i, Om_lamb_i

  real(kind=8) :: w_phi_i, ratio_kin_phi
  
  real(kind=8), allocatable, dimension(:) :: tmp_real_vec
  real(kind=8) :: tmp_real

  real(kind=8) :: t_stop
  real(kind=8) :: eps

  integer :: num_output_x
  integer :: num_output_t
  integer, allocatable, dimension(:) :: tmp_int_vec

  character(len=20), allocatable, dimension(:) :: output_names, output_cosmo_names
  character(len=20) :: tmp_char

  character(len=10) :: slicing = '', shape_phi= ''
  character(len=10) :: shape_matter = ''
  character(len=100), allocatable, dimension(:) :: tmp_char_vec

  character(len=20) :: potential_filename 
  
  type scalar_potential
     character(len=30) :: name
     character(len=100) :: func, inv, der
     character(len=10), allocatable, dimension(:) :: param, var
     character(len=100), allocatable, dimension(:) :: init_param_func
     real(kind=8), allocatable, dimension(:) :: var_val
  end type scalar_potential

  type(scalar_potential) :: potential_phi

contains
  subroutine input_list
    use IniFile
    implicit none
    
    ! assign_ functions only takes array arguments
    ! tmp provided as a temporary 1 element array
    ! tmp_vec provided as a temporaray allocatable array - must be deallocated after each assignement !!!

    call read_input('input.ini')

    call assign_str(tmp_char,output_names,'output')
    call assign_str(tmp_char,output_cosmo_names,'output_cosmo')
    call assign_str(slicing,tmp_char_vec,'slicing'); deallocate(tmp_char_vec)
    call assign_real(eps,tmp_real_vec,'eps'); deallocate(tmp_real_vec)

    call assign_integer(num_output_x,tmp_int_vec,'number_of_x');deallocate(tmp_int_vec)
    call assign_integer(num_output_t,tmp_int_vec,'number_of_t');deallocate(tmp_int_vec)
    
    call assign_real(t_stop,tmp_real_vec,'t_stop');deallocate(tmp_real_vec)

    call assign_real(H0,tmp_real_vec,'H0');deallocate(tmp_real_vec)

    call assign_real(ai,tmp_real_vec,'ai');deallocate(tmp_real_vec)
    call assign_real(Hi,tmp_real_vec,'Hi');deallocate(tmp_real_vec)

    call assign_real(Om_m_i,tmp_real_vec,'Om_mi');deallocate(tmp_real_vec)
    call assign_real(Om_phi_i,tmp_real_vec,'Om_phii');deallocate(tmp_real_vec)
    Om_lamb_i = 1d0 - Om_m_i - Om_phi_i + epsilon(0d0)

    call assign_real(w_phi_i,tmp_real_vec,'w_phi_i');deallocate(tmp_real_vec)
    ratio_kin_phi = 0.5 * ( 1d0 + w_phi_i )
    

    call assign_real(delta_alphai,tmp_real_vec,'delta_alphai');deallocate(tmp_real_vec)
    call assign_real(x0_delta_alphai,tmp_real_vec,'x0_delta_alphai');deallocate(tmp_real_vec)
    call assign_real(sig_delta_alphai,tmp_real_vec,'sig_delta_alphai');deallocate(tmp_real_vec)

    call assign_str(shape_matter,tmp_char_vec,'shape_matter');deallocate(tmp_char_vec)
    call assign_real(delta_mi,tmp_real_vec,'delta_mi');deallocate(tmp_real_vec)
    call assign_real(xmax_delta_mi,tmp_real_vec,'xmax_delta_mi');deallocate(tmp_real_vec)
    call assign_real(steepness_mi,tmp_real_vec,'steepness_mi');deallocate(tmp_real_vec)

    call assign_str(shape_phi,tmp_char_vec,'shape_phi');deallocate(tmp_char_vec)
    call assign_real(delta_phii,tmp_real_vec,'delta_phii');deallocate(tmp_real_vec)
    call assign_real(x0_delta_phii,tmp_real_vec,'x0_delta_phii');deallocate(tmp_real_vec)
    call assign_real(sig_delta_phii,tmp_real_vec,'sig_delta_phii');deallocate(tmp_real_vec)
    call assign_real(xmax_delta_phii,tmp_real_vec,'xmax_delta_phii');deallocate(tmp_real_vec)

    call assign_str(potential_filename,tmp_char_vec,'potential_filename');deallocate(tmp_char_vec)

    deallocate(input)
    
    call input_potential_list
        
  end subroutine input_list

  subroutine input_potential_list
    use Inifile
    implicit none
    integer :: pot_ii
    
    call read_input(potential_filename)

    call assign_str(potential_phi%name,tmp_char_vec,'name');deallocate(tmp_char_vec)
    call assign_str(potential_phi%func,tmp_char_vec,'func');deallocate(tmp_char_vec)
    call assign_str(potential_phi%inv,tmp_char_vec,'inverse');deallocate(tmp_char_vec)
    call assign_str(potential_phi%der,tmp_char_vec,'derivative');deallocate(tmp_char_vec)
    
    call assign_str(tmp_char,potential_phi%param,'parameters')

    allocate(potential_phi%init_param_func(size(potential_phi%param)))

    do pot_ii=1, size(potential_phi%param)
       call assign_str(potential_phi%init_param_func(pot_ii),tmp_char_vec,potential_phi%param(pot_ii));deallocate(tmp_char_vec)
    enddo
    
    if(size(potential_phi%param)/=size(potential_phi%init_param_func))then
       print*, 'Error : number of parameters in potential differs to the number of init. functions provided '
       stop
    endif

    call assign_str(tmp_char,potential_phi%var,'variables')

    deallocate(input)

  end subroutine input_potential_list

end module input_param
