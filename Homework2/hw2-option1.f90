Module rocket_params_vars
  Implicit none
  integer, parameter :: no_eqs = 2
  real*8, parameter :: p_a = 101325.d0,  &
                         n = ... ,    &

                       ......

                      grav = 9.81d0,    &
                        pi = 4.*atan(1.d0), &
                    A_star = ...

  real :: p_c, p_c_dot, &
          r, r_dot, &
          p_c_old, r_old, &
          eta_crit

  contains

  Function f_cor()
   ! perimeter factor
   Implicit none
    real*8 :: f_cor, eta
     eta =  ...
     f_cor = 1.d0
     if( eta < ... )then
         f_cor = ...
     endif
   return
  End


  Function I_sp()
   ! specific impulse
   Implicit none
     real*8 :: I_sp
       I_sp = ...
     return
  End function   
 
  Function ODEs( t ) result ( k )
   Implicit none
    real*8 :: t
    real*8, dimension( no_eqs ) :: k
     r_dot = a * p_c**n
     k( 1 ) =  ...
     k( 2 ) =  ...
   return
  End function

  Subroutine RK4( t, del_t )
   Implicit none
    integer, parameter :: no_stages = 4
    real*8 :: t, del_t
    real*8, dimension( no_eqs )  ::  k( no_stages, no_eqs )
!.. soln and slopes at interval beg
     p_c_old = p_c
     ...
     k( 1, : ) = ODEs( t )
!.. slopes at p_1 = 0.5
     p_c = p_c_old + 0.5d0 * ...
     ...
!.. slopes at p_2 = 0.5
     ...
!.. slopes at p_3 = 1.0
    ...
!.. update for end of interval
    ...
!.. update time
    ...
   return
  End subroutine

End module

Program Rocket_Perf
 Use rocket_params_vars
 Implicit none
 real*8 ::  V_c, ...
 real*8 ::  dt, time
 integer :: nstep = 0

 print*,'enter eta (fraction of effective propellant thickness beyond wich burn diminishes) : '
 read*, eta_crit
 
 print*,'enter time step size, dt [s] : '
 read*, dt

!.. initial params
  p_c =  ...    ! chamber pressure
  r = ...      ! chamber radius
  V_c = ...  ! chamber volume
  V_p = ...  ! propellant volume
  m_p = ...  ! propellant mass
  m_p_init = ...

  open( 1, file = 'rocket.dat',form='formatted') 
  write(1,*)' "t [s]" "p_c [MPa]" ... '
  time = 0.d0
  nstep = 0

  do while(  m_p / m_p_init > 0.05d0 )

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


