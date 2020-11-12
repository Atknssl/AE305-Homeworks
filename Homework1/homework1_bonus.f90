!--------------------------------------------------
!..Homework 1 Bonus
!..Team 31
!..AE305 - Numerical Methods
!
! This program uses RK2 method and Trapezoid method
! to calculate ground roll distance of an aircraft
!--------------------------------------------------
Module data
    implicit none
    real,parameter:: WParea=29.24, takeoffweight=88250, TmaxSL=16256, nrofengines=2, CD=0.2150, CL=1.792, friccoeff=0.02
    real,parameter :: airdensitySL=1.225, airdensity1000m=1.112, airdensity2000m=1.007, grav=9.81
    real :: rho, thrust
 
    contains
      Function Liftf(vel) ! Calculate Lift
       real :: vel, Liftf
       Liftf = CL*(1./2)*rho*(vel**2)*WParea
       return
       End Function Liftf
    
      Function Dragf(vel) ! Calculate Drag
       real :: vel, Dragf
       Dragf = CD*(1./2)*rho*(vel**2)*WParea
       return
       End Function Dragf
 
      Function Thrustf(dens) ! Calculate Thrust
       real :: dens, Thrustf
       Thrustf = TmaxSL*(dens/airdensitySL)
       return
       End Function Thrustf
 
       Function ODE(vel) ! Slope function for velocity
          real:: vel
          real :: ODE
          ODE  = (thrust-Dragf(vel)-friccoeff*(takeoffweight-Liftf(vel)))*(grav/takeoffweight)
          return
        End
 End module
 
 PROGRAM HOMEWORK1_BONUS
    use data
 
    implicit none
 
    real :: stepsize, k1, k2, oldvelocity, velocity, velocityx, lift, time
    real :: integration = 0
    real, parameter:: a1 =1./4 , a2=3./4, p1=2./3
    integer :: selectedaltitude
 
 !  Initial values
    time=0.
    velocity=0.
    oldvelocity = 0.
    lift=0.


 !  Select altitude
    write(*,'(/,a)')'Select one of the available altitudes:'
    write(*,'(a)')' [1] 0 m (Sea Level)'
    write(*,'(a)')' [2] 1000 m '
    write(*,'(a)')' [3] 2000 m '
    write(*,'(a)',advance='no')'Input :> '
    read(*,*) selectedaltitude
 
 !  Set proper variables for selected altitude
    Select case(selectedaltitude)
       case(1)
          write(*,*) 'Selected Sea Level'
          ! Set density and thrust for sea level
          rho = airdensitySL
          thrust = TmaxSL*nrofengines
       case(2)
          write(*,*) 'Selected 1000 m'
          ! Set density and thrust for 1000 m
          rho = airdensity1000m
          thrust = Thrustf(rho)*nrofengines
       case(3)
          ! Set density and thrust 2000 m
          write(*,*) 'Selected 2000 m'
          rho = airdensity2000m
          thrust = Thrustf(rho)*nrofengines
       case default
          ! Validation check
          write(*,*) 'Unidentified input, stopping program...'
          stop
    end Select
 
 !  Select time step
    write(*,'(/,(a))',advance='no')'Enter Time Step :> '
    read(*,*) stepsize

 !  Use RK2 Method to find velocity and trapezoidal integration method to find distance 
       do while ( lift .lt. takeoffweight )
          oldvelocity = velocity;                             
          time = time + stepsize                                 ! New time is defined
          k1 = ODE(velocity)                                     ! k1 is calculated which is slope at previous velocity
          velocityx = velocity + k1*p1*stepsize                  ! velocityx is calculated
          k2 = ODE(velocityx)                                    ! k2 is calculated which is slope at velocityx
          velocity = velocity + ( a1*k1 + a2*k2 )* stepsize      ! New velocity is defined

          integration = integration + (1./2)*(oldvelocity + velocity)*stepsize  ! Integral is calculated for every interval
          
          lift = Liftf(velocity)                                 ! Calculate next lift to check if lift off occured
       enddo
       
       write(*,'(a,f12.3)') 'Minumun distance: ', integration
 !..Close the output file
    close(1)
 
    stop
 END PROGRAM HOMEWORK1_BONUS
 