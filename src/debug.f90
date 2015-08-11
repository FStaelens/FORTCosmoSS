!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
program debug
  !
  ! the main program that includes the main loop for the numerical integration in time
  ! as numerical infinities might have physical meaning, the code is designed to keep
  ! on running, should he encounter a NaN and to issue a warning containing the time
  ! at which such event was encountered at the end of the integration.
  !
  ! the process of output are made through a call to "write_output" which calls routines
  ! from "output.f90"
  !
  ! a rudimentary log entry is written at each step using "write_in_log". The content
  ! of the log can be easily eloborated upon by changing this procedure
  !
use grid
use init
use outputs
use constraints
use gauge_choice
use evolution
use sources
implicit none

integer :: iter = 0
integer :: output_iter = 0
logical :: dir_exists = .false., file_exists = .false., init_success = .true.
integer :: precent
integer :: mod_output_t, mod_output_x

type bug
   logical :: exists
   integer :: index
   real(kind=8) :: time  
   character(len=7) :: name
end type bug

type(bug) :: found_NaN 

found_NaN%name = ''
found_NaN%exists = .false.
found_NaN%index = 0
found_NaN%time = 0d0

call init_data(init_success)
call H_cons
call M_cons
call H_cons_cosm

inquire(file='./data/data0000.dat', exist=dir_exists)
if(dir_exists.eqv..true.) call system('rm -rf ./data')
inquire(file='./output_cosmo.dat', exist=file_exists)
if(file_exists.eqv..true.) call system('rm output_cosmo.dat')


mod_output_t = int((t_stop/dt)/num_output_t)+1
mod_output_x = int((x_grid(nx)/dx)/num_output_x)+1

call system('mkdir data')
call write_output

if(init_success.eqv..false.)stop
open(10,file='debug.log')


! time evolution
do while(t<t_stop)
   t=t+dt
   iter=iter+1
   call PIRK
   call build_matter_sources
   call H_cons
   call M_cons
   call H_cons_cosm
   if(isnan(t_sync).eqv..false.)then
      call progress(t/t_stop,precent)
   else
      found_NaN%name = 't_sync'
      found_NaN%exists = .true.
      found_NaN%index = iter
      found_NaN%time = t
      exit
   endif
   if(mod(iter,mod_output_t).eq.0)then
      output_iter=output_iter+1
      call write_in_log
      call write_output
   end if
enddo
if((t_stop/=0).and.(found_NaN%name/='t_sync'))then 
   if(precent/=100)call progress(t/t_stop,precent)
endif
call sleep(1)
print*,''
close(10)

if(found_NaN%exists) then
   print*,'Warning : NaN found in "',trim(adjustl(found_NaN%name)) ,'" at i = ', found_NaN%index, 't = ', found_NaN%time  
else
   print*,'*Message : Programme terminated without error'
endif
print*, '--> Find output in "./data/" and log in "debug.log"'
print*,''

contains
  subroutine write_output
    use input_param
    implicit none
    character(len=40) :: header, file_name
    logical :: warning = .true.
    integer :: i, j

    call output_list
    write(header,'(a,f10.3,a,f10.3)') 't = ', t,' / t_sync = ', t_sync
    write(file_name,'(a,i4.4,a)') './data/data',output_iter,'.dat'
    call print_list(output_names,trim(adjustl(file_name)),header,warning,mod_output_x)
    warning = .false. ! to issue a warning message only on first call to print_list if needed
    write(header,'(a,f10.3)') 'H0 = ', H0
    call append_to_list(output_cosmo_names,'./data/output_cosmo.dat',header)

    if(found_NaN%exists.eqv..false.) then
       do j = 1, size(output_names)
          do i = 1, size(output_mat(j)%val)
             if(isnan(output_mat(j)%val(i))) then
                found_NaN%name = output_names(j)
                found_NaN%exists = .true.
                found_NaN%index = iter
                found_NaN%time = t
             end if
          end do
       end do
    endif
  end subroutine write_output

  subroutine write_in_log
    implicit none
    write(10,'(a,i9,a,f10.3,a,f10.3)')' --> output i =', iter, ' time =', t, ' cosmic_time =', t_sync
  end subroutine write_in_log

end program debug
