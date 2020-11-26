Module rocket_params_vars
  Implicit none
  integer, parameter :: no_eqs = 2
  real*8, parameter :: p_a = 101325.d0, &
                     rho_p = 1140.d0, &
                         n = 0.305d0, &
                         a = 0.0000555d0, &
                       r_0 = 5.d0, &
                       r_f = 15.d0, &
                         L = 1.25d0, &
                       T_c = 2810.d0, &
                         R = 365.d0, &
                      gama = 1.25d0, &
                      grav = 9.81d0, &
                        pi = 4.*atan(1.d0)
                        ! A_star tanımlamadım
  real*8 :: p_c, p_c_dot, &
          r, r_dot, &
          p_c_old, r_old, &
          eta_crit, A_star

  contains

  Function m_n_dot(p_c_) 
    real*8 :: m_n_dot, p_c_
    m_n_dot = p_c_*A_star*sqrt(gama/(R*T_c))*((gama+1)/2)**(-(gama+1)/(2*(gama-1)))
    return
  End Function m_n_dot 

  Function rho_c(p_c_)
    real*8 :: rho_c, p_c_
    rho_c = p_c_/(R*T_c)
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
     else if (eta >= 0.or.eta <= 0.15)then
         f_cor = 1-exp(-7*eta)
     end if
   return
  End Function f_cor

  Function I_sp(p_c_)
   ! specific impulse
   Implicit none
     real*8 :: I_sp, p_c_
       I_sp = (1/grav)*sqrt(((2*gama*R*T_c)/(gama-1))*(1-(p_a/p_c_)**((gama-1)/gama)))
     return
  End function      
 
  Function ODEs( x ) result ( k )
   Implicit none
    real*8 :: x
    real*8, dimension( no_eqs ) :: k
     !r_dot=a*p_c**n (bu neden verildi anlamadım çünkü r_dot neden k2 içinde sabit kalsın?)
     k( 1 ) =  a*x**n
     k( 2 ) =  R*T_c(((fcor(r)*2*k(1))/r)*(rho_p-rho_c(x))-m_n_dot(x)/(pi*L*r**2))
   return
  End function

  Subroutine RK4( t, del_t )
    Implicit none
    integer, parameter :: no_stages = 4
    real*8 :: t, del_t
    real*8, dimension( no_eqs )  ::  k( no_stages, no_eqs ), phi(no_eqs)
    p_c_old = p_c
    r_old = r
    k(1, : ) = ODEs( t )
!.. slopes at p_1 = 0.5
    p_c = p_c_old + 0.5d0 * del_t * k(1, 2)
    r = r_old + 0.5d0 * del_t * k(1, 1)
    k(2, : ) = ODEs( p_c )

!.. slopes at p_2 = 0.5
    p_c = p_c_old + 0.5d0 * del_t * k(2, 2)
    r = r_old + 0.5d0 * del_t * k(2, 1)
    k(3, : ) = ODEs( p_c )

!.. slopes at p_3 = 1.0
    p_c = p_c_old + del_t * k(3, 2)
    r = r_old + del_t * k(3, 1)
    k(4, : ) = ODEs( p_c )

    phi(:) = (1.0/6)* (k(1, : )+2*k(2, : )+2*k(3, : )+k(4, : ))         !weighted slope is calculated
    p_c = p_c_old + phi(2) * del_t                                    !p_c is updated
    r = r_old + phi(1) * del_t 

    t = t + del_t                               !time is updated
  return
 End subroutine

End module

Program Rocket_Perf
 Use rocket_params_vars
 Implicit none
 character*40 :: fname
 real*8 ::  V_c, ...
 real*8 ::  dt, time
 integer :: nstep = 0
 real*8 :: th_radius

 write(*,'(a)') 'Enter throat radius in centimeters :'
 write(*,'(a)',advance='no') ':>'
 read*, th_radius
 
A_star = pi * th_radius**2

 print*,'enter time step size, dt [s] : '
 read*, dt

!.. initial params
  p_c =  ...    ! chamber pressure
  r = ...      ! chamber radius
  V_c = ...  ! chamber volume
  V_p = ...  ! propellant volume
  m_p = ...  ! propellant mass
  m_p_init = ...

  write(*,'(a)',advance='no')' Enter the output file name [rocket.dat]:>'
	read(*,'(a)') fname
	if( fname .eq. ' ') fname = 'rocket.dat'
  open(1 ,file=fname, form='formatted')

  write(1,*)' "t [s]" "p_c [MPa]" ... ' ! will be completed later

  time = 0.d0
  nstep = 0

  do while( m_p / m_p_init > 0.05d0 )

    nstep = nstep + 1

    ! update solution and time
    call RK4( ... )

    ! check chamber pressure... abort if p_c < p_a
     ...

    ! compute mass of propellant left
     ...

    ! print on screen and store soln
    if( nstep == 1 .or. mod(nstep,5)==0) &
      write(1,'(8(2x,e12.6))') time, p_c*1.d-6, ...
      write(*,'(8(a,e9.3))')' t = ',time,' s,   p_c = ', p_c*1d-6, ...

  enddo

  close(1)

  stop
End


