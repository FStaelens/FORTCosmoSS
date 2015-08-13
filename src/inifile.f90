!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module IniFile
  ! includes the necessary routines to read the variables from an input file placed in the parent directory
  ! 
  ! read_input    input  : "input_file" (str) the name of the file to be read from the parent directory
  !               output :  assigned directly to variable "input" with custom type "inputtype"
  !
  ! assign_real   input  : "var_array" (real) an allocatable array,
  !                        "name" (str) the name of the variable as it appears in the input file
  !               output : "var_array" (real) if numb. of inputs > 1 : contains the necessary input values
  !                                           if numb. of inputs = 1 : is an array of zero dimension
  !                        "var" (real) if numb. of inputs > 1 : not allocated
  !                                     if numb. of inputs = 1 : contains the necessary input value
  ! 
  ! assign_integer  --- similar to "assign_real" for integer inputs
  !
  ! assign_str      --- similar to "assign_real" for string inputs
  !
  ! isthere      input  : "name" (str) the name of some desired variable
  !              output : "isthere" (logical) is ".true." if input with appropriate name is found
  
  implicit none  
  type inputtype
     character(len=30) :: name
     character(len=100), allocatable, dimension(:) :: val 
  end type inputtype
  type(inputtype), allocatable, dimension(:) :: input

contains 
  subroutine read_input(input_file)
    implicit none
    character(len=*), intent(in) :: input_file
    character(len=2048) :: buffer, label, str, fmt 
    integer :: pos, ios, cnt, i, numb_input

    ios = 0
    cnt = 0
    i = 1
    numb_input = 0
    
    open(10,file=input_file,status='old')

    do while (ios == 0) ! count the number of inputs
       read(10, '(A)', iostat=ios) buffer
       if(ios == 0)then
          
          ! Find the first instance of whitespace. Split label and data.
          pos = scan(buffer, ' ')
          label = buffer(1:pos)
          
          if((label=='') .or. (label( :1)=='#'))cycle  ! skip line if blank or commented with '#'
          numb_input = numb_input + 1 
       endif
    end do
    rewind(10)
    ios = 0  

    allocate(input(numb_input))  
  
    do while (ios == 0) ! read input
       read(10, '(A)', iostat=ios) buffer
       if(ios == 0)then
          
          ! Find the first instance of whitespace. Split label and data.
          pos = scan(buffer, ' ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)
          
          if((label=='') .or. (label( :1)=='#'))cycle  ! skip line if blank or commented with '#'

          str = buffer

          fmt = '(a'
          cnt = 0
          do while((pos/=0).and.(pos/=1))
             pos = scan(trim(adjustl(str)), ' ')
             str = trim(adjustl(str(pos+1:)))
             cnt = cnt + 1
             fmt = trim(adjustl(fmt)) // ',X,a'
          enddo

          fmt = trim(adjustl(fmt)) // ')'
          
          input(i)%name = label
          allocate(input(i)%val(cnt))
          if(cnt==1)then
             read(buffer, '(a)', iostat=ios) input(i)%val
          else
             read(buffer, trim(adjustl(fmt)), iostat=ios) input(i)%val
             read(buffer, *, iostat=ios) input(i)%val
          end if
          !print*,trim(adjustl(fmt))
          !print*,input(i)%val
          
          i=i+1
          
       endif
    end do
    close(10)

  end subroutine read_input

  subroutine assign_real(var,var_array,name)
    implicit none
    real(kind=8), intent(out) :: var
    real(kind=8), allocatable, dimension(:), intent(inout) :: var_array
    character(len=*), intent(in) :: name
    integer :: i, j
    logical :: flag

    flag = .false.
    do i=1, size(input)
       if(input(i)%name==name) then 
          allocate(var_array(size(input(i)%val)))
          do j=1, size(input(i)%val)
             !print*,input(i)%name,trim(adjustl(input(i)%val(j)))
             read(input(i)%val(j),*) var_array(j)
          enddo
          flag = .true.
       end if
    enddo
    if(flag.eqv..false.) then
       allocate(var_array(0))
       print*, 'Warning : no matching of label "', trim(adjustl(name)), '" found in input file'
    endif
    
    var = 0d0
    if(size(var_array)==1) then
       var = var_array(1)
    endif
  end subroutine assign_real

  subroutine assign_integer(var,var_array,name)
    implicit none
    integer, intent(out) :: var
    integer, allocatable, dimension(:), intent(inout) :: var_array
    character(len=*), intent(in) :: name
    integer :: i, j
    logical :: flag

    flag = .false.
    do i=1, size(input)
       if(input(i)%name==name) then
          allocate(var_array(size(input(i)%val)))
          do j=1, size(input(i)%val)
             read(input(i)%val(j),*) var_array(j) 
          enddo
          flag = .true.
       end if
    enddo
    if(flag.eqv..false.)then
       allocate(var_array(0))
       print*, 'Warning : no matching of label "', trim(adjustl(name)), '" found in input file'
    endif

    var = 0
    if(size(var_array)==1) then
       var = var_array(1)
    endif
  end subroutine assign_integer

  subroutine assign_str(var,var_array,name)
    implicit none
    character(len=*), intent(out) :: var
    character(len=*), allocatable, dimension(:), intent(inout) :: var_array
    character(len=*), intent(in) :: name
    integer :: i
    logical :: flag

    flag = .false.
    do i=1, size(input)
       if(input(i)%name==name) then          
          allocate(var_array(size(input(i)%val)))
          var_array = input(i)%val 
          flag = .true.
       end if
    enddo
    if(flag.eqv..false.) then
       allocate(var_array(0))
       print*, 'Warning : no matching of label "', trim(adjustl(name)), '" found in input file'
    endif
    
    var = ''
    if(size(var_array)==1) then
       var = var_array(1)
    endif

  end subroutine assign_str

  function isthere(name)
    implicit none
    logical :: isthere
    character(len=*) :: name
    integer :: i
    
    isthere = .false.
    do i=1, size(input)
       if(input(i)%name==name) then
          isthere = .true.
       end if
    end do
  end function isthere

end module IniFile
