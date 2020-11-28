!--------------------------------------------------
!..Homework 2 Bonus
!..Team 31
!..AE305 - Numerical Methods
!
!..This program calculates the chamber pressure, 
!..burn rate and specific impulse for given rocket
!..using RK4 method with adaptive time step and
!..outputs results to a file.
!--------------------------------------------------
Module rocket_params_vars_bonus
  Implicit none
  integer, parameter :: no_eqs = 2
  real*8, parameter :: p_a = 101325.d0, &
                     rho_p = 1140.d0, &
                         n = 0.305d0, &
                         a = 0.0000555d0, &
                       r_0 = 0.05d0, &
                       r_f = 0.15d0, &
                         L = 1.25d0, &
                       T_c = 2810.d0, &
                     R_cst = 365.d0, &
                      gama = 1.25d0, &
                      grav = 9.81d0, &
                        pi = 4.*atan(1.d0)
  real*8 :: p_c, p_c_dot, &
          r, r_dot, &
          A_star

  contains

  Function m_n_dot(p_c_) 
    ! calculate m_n dot to use in ODE
    real*8 :: m_n_dot, p_c_
    m_n_dot = p_c_*A_star*sqrt(gama/(R_cst*T_c))*((gama+1)/2)**(-(gama+1)/(2*(gama-1)))
    return
  End Function m_n_dot 

  Function rho_c(p_c_)
    ! calculate rho_c to use in ODE
    real*8 :: rho_c, p_c_
    rho_c = p_c_/(R_cst*T_c)
    return
  End Function rho_c

  Function f_cor(rad)
   ! perimeter factor
   Implicit none
    real*8 :: f_cor, eta, rad
     eta = (r_f-rad)/(r_f-r_0)
     f_cor = 1.d0
     if( eta < 0 )then
         f_cor = 0
     else if (eta <= 0.15)then
         f_cor = 1-exp(-7*eta)
     end if
   return
  End Function f_cor

  Function I_sp(p_c_)
   ! specific impulse
   Implicit none
    real*8 :: I_sp, p_c_
    ! Here, absolute of value inside square root is taken because
    ! there is really small error caused by numerical computing makes
    ! value inside of the square root negative for some intervals
    ! which results in runtime error. This negative value is so small
    ! that its practically zero. Therefore, absolute can be used without 
    ! effecting accuracy of result. Using abs solved the runtime error.
    I_sp = (1/grav)*sqrt(abs(((2*gama*R_cst*T_c)/(gama-1))*(1-((p_a/p_c_)**((gama-1)/gama)))))
    return
  End function I_sp   
 
  Function ODEs( p_c_ , r_) result ( k )
   Implicit none
    real*8 :: p_c_, r_
    real*8, dimension( no_eqs ) :: k
     r_dot = a*p_c_**n
     k( 1 ) =  r_dot
     k( 2 ) =  R_cst*T_c*(((f_cor(r_)*2* a * p_c_ ** n)/r_)*(rho_p-rho_c(p_c_))-m_n_dot(p_c_)/(pi*L*r_**2))
   return
  End function

  Subroutine RK4( t, del_t )
    Implicit none
    integer, parameter :: no_stages = 4
    real*8 :: t, del_t
    real*8, dimension( no_eqs )  ::  result(no_eqs)

    ! RK4 has to be calculated such that after the calculation
    ! p_c shouldn't change. Because rk4 is used first while trying
    ! to determine the step size. If p_c changes while calculating 
    ! step size, calculation will be wrong. Therefore, in the bonus
    ! part, RK4 must be calculated in a different function, which
    ! does not change the inputs value. This is why this part of 
    ! code is different than the main code.
    result = RK4_(p_c,r,del_t)

    r = result(1) ! r is updated
    p_c = result(2) ! p_c is updated

    t = t + del_t !time is updated
  return
 End subroutine

  function RK2(p_c_,r_,del_t) result(res)
    implicit none
    integer, parameter :: no_stages = 2
    real*8 ::  p_c_, del_t ,r_
    real*8, dimension( no_eqs )  ::  k( no_stages, no_eqs ) , res(no_eqs)
    k(1,:) = ODEs(p_c_ , r_)

    res(1) = r_ + 0.5d0 * del_t * k(1,1)
    res(2) = p_c_ + 0.5d0 * del_t * k(1,2)
  
    k(2,:) = ODEs(res(2), res(1))
    res(1) = r_ + k(2,1) * del_t
    res(2) = p_c_ + k(2,2) * del_t
    return
  end function

  function RK4_(p_c_, r_, del_t) result(res)
    implicit none
    integer, parameter :: no_stages = 4
    real*8 ::  p_c_, del_t ,r_
    real*8, dimension( no_eqs )  ::  k( no_stages, no_eqs ), phi(no_eqs), res(no_eqs)
    
    k(1,:) = ODEs(p_c_, r_)

    res(2) = p_c_ + 0.5d0 * del_t *k(1,2)
    res(1) = r_ + 0.5d0 * del_t * k(1,1)
    
    k(2,:) = ODEs(res(2), res(1))

    res(2) = p_c_ + 0.5d0 * del_t * k(2,2)
    res(1) = r_ + 0.5d0 * del_t * k(2,1)

    k(3,:) = ODEs(res(2), res(1))

    res(2) = p_c_ + del_t * k(3,2)
    res(1) = r_ +  del_t * k(3,1)

    k(4,:) = ODEs(res(2), res(1))

    phi(:) = (1.0/6)* (k(1, : )+2*k(2, : )+2*k(3, : )+k(4, : )) 
    
    res(2) = p_c_ + phi(2) * del_t
    res(1) = r_ + phi(1) * del_t
    
    return
  end function
End module

Program Rocket_Perf
 Use rocket_params_vars_bonus
 Implicit none
 character*40 :: fname
 real*8 ::  dt, time, isp
 integer :: nstep = 0
 real*8 :: th_radius

 ! Variables for the bonus part
 real*8, parameter :: E_allowed = 0.0001 ! Allowed error for adaptive stepsize
 real*8 :: E_o
 real*8, dimension(no_eqs) :: dummyRK4(no_eqs), dummyRK2(no_eqs)

 write(*,'(a)') 'Enter throat radius in centimeters :'
 write(*,'(a)',advance='no') ':> '
 read*, th_radius
 
A_star = pi * (th_radius/100)**2

 print*,'enter the initial time step size, dt [s] : '
 read*, dt

!.. initial params
  p_c =  p_a    ! chamber pressure
  r = r_0      ! chamber radius
  isp = I_sp(p_c) ! Specific Impulse

  write(*,'(a)',advance='no')' Enter the output file name [rocket.dat]:>'
  read(*,'(a)') fname
  if( fname .eq. ' ') fname = 'rocket.dat'
  open(1 ,file=fname, form='formatted')
  write(1,*)' "t [s]" "p_c [MPa]" "I_sp[s]" "r_dot[m/s]"'

  time = 0.d0
  nstep = 0

  do while( p_c >= p_a )
    nstep = nstep + 1

    ! Dummy variables to hold temporary solutions for step size calculation
    dummyRK4 = RK4_(p_c,r,dt)
    dummyRK2 = RK2(p_c,r,dt)
    
    ! Calculate local truncation error
    E_o = (dummyRK4(2) - dummyRK2(2))/dummyRK4(2)

    ! Step size is updated
    dt = dt * (abs(E_allowed / E_o))**0.20d0    
    
    ! update solution and time
    call RK4(time, dt)

    ! calculate isp
    isp = I_sp(p_c)

    ! print on screen and store soln
    if( nstep == 1 .or. mod(nstep,5)==0) &
      write(1,'(8(4x,e12.6))') time, p_c*1.d-6 ,isp, r_dot
      write(*,'(8(a,e9.3))')' t = ',time,' s,   p_c = ', p_c*1d-6,' MPa, I_sp = ', isp, ' s r_dot = ',r_dot,' m/s'
  enddo
  
  ! Print the interval number for bonus part
  write(*,*) 'Number of time intervals needed = ', nstep
  close(1)

  stop
End


