!
! Copyright (c) 2015, Jeremy Rekier. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
! Last edited : 2015-08-10
!
! --------------------------------------------------------------------------
module profiles
  !
  ! includes mathematical function profiles
  ! any new radial distribution profile should be added here
  !
  implicit none
contains

  function logistic(x,d,x0,k)
    implicit none
    real(kind=8) :: x, logistic, d, x0, k

    logistic = d*(1d0-tanh(k*(x-x0)))/2d0 !d/(1+exp(2*k*(x-x0)))

  end function logistic

  function bump(x,d,x0)
    implicit none
    real(kind=8) :: x, bump, d, x0
    bump = 0d0
    
    if (x<x0) then
       bump = d*exp(-x**2/(x0**2-x**2))
    else 
       bump = 0d0
    end if
  end function bump

  function dxbump(x,d,x0)
    implicit none
    real(kind=8) :: x, dxbump, d, x0, sig
    dxbump = 0d0
    
    if (x<x0) then
       dxbump = (-2d0*x*x0**2)/(x**2-x0**2)**2*bump(x,d,x0)
    else
       dxbump = 0d0
    endif
  end function dxbump
  
  function sym_gaussian(x,d,xc,sig)
    implicit none
    real(kind=8) :: x, sym_gaussian, d, xc, sig
    
    sym_gaussian = d*x**2/(1d0+x**2)*(exp(-(x-xc)**2/sig**2)+exp(-(x+xc)**2/sig**2))
    
  end function sym_gaussian

  function dxsym_gaussian(x,d,xc,sig)
    implicit none
    real(kind=8) :: x, dxsym_gaussian, d, xc, sig
    
    dxsym_gaussian = (2*d*x/(1+x**2)**2)*(exp(-(x-xc)**2/sig**2)+exp(-(x+xc)**2/sig**2))&
        -(2*d*x**2)/(1+x**2)*((x-xc)/sig**2*exp(-(x-xc)**2/sig**2)+(x+xc)/sig**2*exp(-(x+xc)**2/sig**2))
    
  end function dxsym_gaussian
end module profiles

module root_finding
  !
  ! includes a root-finding algorithm as the function "root(funcd,x1,x2,xacc)"
  !  - funcd(x,f,df) must be a function taking a real value (x) and returning
  !    the value of a mathematical function evaluated at this point (f) + its derivative (df)
  !
implicit none
contains
  function root(funcd,x1,x2,xacc)   ! root-finding algorithm
    implicit none
    real(kind=8) :: root, x1, x2, xacc
    interface
       subroutine funcd(x,f,df)
         real(kind=8), intent(in) :: x
         real(kind=8), intent(out) :: f, df
       end subroutine funcd
    end interface
    
    real(kind=8), parameter :: EPS = 3.e-8
    integer, parameter :: maxit = 100
    integer :: j
    real(kind=8) df, dx, dxold, f, fh, fl, temp, xh, xl
    
    call funcd(x1,fl,df)
    call funcd(x2,fh,df)

    if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
       print*, '*Warning : the root of the zero-finding function might not be bracketted'
       print*, '           or the function might have multiple roots'
       print*, ''
    end if
    if(fl.eq.0.)then
       root=x1
       return
    elseif(fh.eq.0.)then
       root=x2
       return
    elseif(fl.lt.0.)then    ! orient the search so that f(x1)<0
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    endif
    root=.5*(x1+x2)       ! initialisation
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(root,f,df)
    do j=1,maxit            ! main loop
       if(((root-xh)*df-f)*((root-xl)*df-f).gt.0..or.abs(2.*f).gt.abs(dxold*df))then
          ! bisect if out of range or too slow
          dxold=dx
          dx=0.5*(xh-xl)
          root=xl+dx
          if(xl.eq.root)return  ! change in root negligible. Newton step acceptable. Take it
       else
          dxold=dx
          dx=f/df
          temp=root
          root=root-dx
          if(temp.eq.root)return
       endif
       if(abs(dx).lt.xacc)return
       call funcd(root,f,df)
       if(f.lt.0.)then
          xl=root
       else
          xh=root
       endif
    enddo
    print*, 'routine exceeding maximum iterations'
    return
  end function root

end module root_finding

module derivatives_fcn
  !
  ! includes the functions for the evaluation of derivatives on a staggered grid
  !
use grid
contains
  function dxf(f)
    !
    ! input  : f (real) an array of (nx+2*off-1) values representing a spatial profile
    ! output : dxf (real) an array of (nx+off) values giving the first derivative of f
    ! 
    ! desription : computes the first derivative using a 4th order centred finite diff. scheme
    !              uses virtual mirror points of negative radius to compute the derivative
    !              at points close to r=0.
    !              uses a backward scheme for points close to the right boundary
    ! 
    implicit none
    real(kind=8), dimension(:), intent(in) :: f(-off+1:nx+off)
    real(kind=8), dimension(:) :: dxf(1:nx+off)
    integer :: i
    
    do i = nx+1, nx+off
       dxf(i)=(25d0*f(i)-48d0*f(i-1)+36d0*f(i-2)-16d0*f(i-3)+3d0*f(i-4))/12d0/dx
    end do

    do i = 1, nx
       dxf(i)=(-f(i+2)+8d0*f(i+1)-8d0*f(i-1)+f(i-2))/12d0/dx
    end do   

  end function dxf

  function d2xf(f)
    !
    ! input  : f (real) an array of (nx+2*off-1) values representing a spatial profile
    ! output : d2xf (real) an array of (nx+off) values giving the first derivative of f
    ! 
    ! desription : computes the second derivative using a 4th order centred finite diff. scheme
    !              uses virtual mirror points of negative radius to compute the derivative
    !              at points close to r=0.
    !              uses a backward scheme for points close to the right boundary
    ! 
    implicit none
    real(kind=8), dimension(:), intent(in) :: f(-off+1:nx+off)
    real(kind=8), dimension(:) :: d2xf(1:nx+off)
    integer :: i

    do i = nx+1, nx+off
       d2xf(i)=(45d0*f(i)-154d0*f(i-1)+214d0*f(i-2)-156d0*f(i-3)+61d0*f(i-4)-10d0*f(i-5))/12d0/dx/dx
    end do

    do i = 1, nx
       d2xf(i)=(-f(i+2)+16d0*f(i+1)-30d0*f(i)+16d0*f(i-1)-f(i-2))/12d0/dx/dx
    end do   

  end function d2xf

  function Delta4_x(u)
    !
    ! input  : u (real) an array of (nx+2*off-1) values representing a spatial profile
    ! output : Delta4_x (real) an array of (nx+off) values giving the first derivative of u
    ! 
    ! desription : computes the 4th centred finite difference operator
    !              uses virtual mirror points of negative radius to compute the derivative
    !              at points close to r=0.
    !              values of Delta4_x(nx+1) to Delta4_x(nx+off) are set to zero
    ! 
    real(kind=8), dimension(-off+1:nx+off) :: u
    real(kind=8), dimension(1:nx+off) :: Delta4_x
    integer :: i
    
    Delta4_x = 0d0

    do i=1, nx
       Delta4_x(i) = u(i+2) - 4d0*u(i+1) + 6d0*u(i) - 4d0*u(i-1) + u(i-2)
    end do
  
  end function Delta4_x

  function Delta4_x_i(u,j)
    !
    ! input  : u (real) an array of (nx+2*off-1) values representing a spatial profile
    !          j (integer) the index at which the 
    ! output : Delta4_x_i (real) an array of (nx+off) values giving the first derivative of u
    ! 
    ! desription : computes the 4th centred finite difference operator the point of index j
    !              uses virtual mirror points of negative radius to compute the derivative
    !              at points close to r=0.
    !              values of Delta4_x(nx+1) to Delta4_x(nx+off) are set to zero and a warning is issued
    ! 
    real(kind=8), dimension(-off+1:nx+off) :: u
    integer :: j
    real(kind=8) :: Delta4_x_i

    if((j>nx+1).or.(j<0))then
       Delta4_x_i = 0d0
       print*, "Warning : invalid index reference during call to Delta4_x_i"
    else
       Delta4_x_i = u(j+2) - 4d0*u(j+1) + 6d0*u(j) - 4d0*u(j-1) + u(j-2)
    end if

  end function Delta4_x_i

end module derivatives_fcn

module boundaries
  !
  ! includes the functions for the treatment of (anti)-symmetric boundary condition around r=0 
  ! and dissipative (Sommerfeld) boundary conditions at the right boundary
  !
use grid
use derivatives_fcn
implicit none
contains
  subroutine symmetrise(f)
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(inout) :: f
    integer :: i
    do i = -off+1, 0
       f(i)=f(-i+1)
    end do
  end subroutine symmetrise
  subroutine anti_symmetrise(f)
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(inout) :: f
    integer :: i
    do i = -off+1, 0
       f(i)=-f(-i+1)
    end do
  end subroutine anti_symmetrise
  subroutine symmetrise_centred(f)
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(inout) :: f
    integer :: i
    do i = -off+1, -1
       f(i)=f(-i)
    end do
  end subroutine symmetrise_centred
  subroutine sommerfeld(f,dtf,f0,dtf0,v)
    !
    ! input  : f (real) an array of (nx+2*off-1) values representing a spatial profile
    !          dtf (real) an array of (nx+2*off-1) values representing the first time der. of f
    !          f0 (real) the asymptotic value of f at spatial infinity
    !          dtf0 (real) the asymptotic value of dtf at spatial infinity
    !          drf (real) an array of (nx+off) values representing the first spatial der. of f
    !          v (real) the signal velocity of the field f
    ! 
    ! output : dtf (real) the array provided as input with modified rightmost values 
    ! 
    ! desription : changes the values of the rightmost values of the array dtf to follow a
    !              sommerfeld boundary condition :  dtf = dtf0 - v*drf - v(f-f0)/r
    ! 
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: f
    real(kind=8), dimension(-off+1:nx+off), intent(inout) :: dtf
    real(kind=8), intent(in) :: f0, dtf0
    real(kind=8), intent(in) :: v
    real(kind=8), dimension(1:nx+off) :: drf
    integer :: i
    drf = dxf(f)
    do i = nx+1, nx+off
       dtf(i)=dtf0-v*drf(i)-v*(f(i)-f0)/x_grid(i)
    end do
  end subroutine sommerfeld
end module boundaries

module num_integration
  !
  ! includes a function to compute the numerical integral using the trapeze function on the spatial domain
  !
  use grid
  implicit none
contains
  function integral(f,dx,a,b)
    real(kind=8), dimension(-off+1:nx+off) :: f
    real(kind=8) :: dx, a, b
    real(kind=8) :: integral
    integer :: i, i_start, i_stop
    
    do i = 1, nx+off
       if(x_grid(i)>a)exit
    end do
    i_start = i
    do i = nx+off, 1, -1
       if(x_grid(i)<b)exit
    end do
    i_stop = i

    integral = 0d0
    do i = i_start+1, i_stop
       integral = integral + 0.5d0*dx*(f(i-1)+f(i))
    end do
  end function integral
end module num_integration

module find_value
  !
  ! includes a function to find the value nearest to some real value "x" within an array "vec"
  !
  implicit none
contains
  subroutine find_closest(vec,x,xi,j)
    use grid
    implicit none
    real(kind=8), dimension(-off+1:nx+off), intent(in) :: vec
    real(kind=8), intent(in) :: x
    real(kind=8), intent(out) :: xi
    integer, intent(out) :: j
    integer :: i
    real(kind=8) :: distance, distance_old

    distance_old = huge(1d0)
    xi = 0d0
    do i = 1, nx+off
       distance = abs(vec(i)-x)
       if(distance.lt.distance_old)then
          xi = vec(i)
          j = i
          distance_old = distance
       end if
    end do
    if(j.gt.nx)j=nx

  end subroutine find_closest
end module find_value

module interpolation
contains
  SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
    !*****************************************************
    !*     Polynomial Interpolation or Extrapolation     *
    !*            of a Discreet Function                 *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*   XA:    Table of abcissas  (N)                   *
    !*   YA:    Table of ordinates (N)                   *
    !*    N:    Number of points                         *
    !*    X:    Interpolating abscissa value             *
    !* OUTPUTS:                                          *
    !*    Y:    Returned estimation of function for X    *
    !*   DY:    Estimated error for Y                    *
    !*****************************************************
    PARAMETER(NMAX=25)
    REAL(kind=8) :: XA(N),YA(N), X,Y,DY
    REAL(kind=8) :: C(NMAX),D(NMAX)
    REAL(kind=8) :: DEN,DIF,DIFT,HO,HP,W
    NS=1
    DIF=DABS(X-XA(1))
    DO I=1,N
       DIFT=DABS(X-XA(I))
       IF(DIFT.LT.DIF) THEN
          NS=I                 !index of closest table entry
          DIF=DIFT
       ENDIF
       C(I)=YA(I)             !Initialize the C's and D's
       D(I)=YA(I)
    END DO
    Y=YA(NS)                 !Initial approximation of Y
    NS=NS-1
    DO M=1,N-1
       DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I) 
          DEN=HO-HP
          IF(DEN.EQ.0.) print*, 'Error: two identical abcissas)'
          DEN=W/DEN
          D(I)=HP*DEN          !Update the C's and D's
          C(I)=HO*DEN
       END DO
       IF(2*NS.LT.N-M) THEN   !After each column in the tableau XA
          DY=C(NS+1)          !is completed, we decide which correction,
       ELSE                   !C or D, we want to add to our accumulating 
          DY=D(NS)            !value of Y, i.e. which path to take through 
          NS=NS-1             !the tableau, forking up or down. We do this
       ENDIF                  !in such a way as to take the most "straight 
       Y=Y+DY	              !line" route through the tableau to its apex,
    END DO                    !updating NS accordingly to keep track of 
                              !where we are. This route keeps the partial
    RETURN                    !approximations centered (insofar as possible)
                              !on the target X.The last DY added is thus the
                              !error indication.
    
  END SUBROUTINE POLINT
end module interpolation
