!--------------------------------------------------
!..Homework 1
!..Team 31
!..AE305 - Numerical Methods
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

PROGRAM HOMEWORK1
   use data

   implicit none

   character*40 :: fname
   real :: stepsize, k1, k2, velocity, velocityx, lift, time
   real, parameter:: a1 =1./4 , a2=3./4, p1=2./3
   integer :: method, selectedaltitude

!  Initial values
   time=0.
   velocity=0.
   lift=0.

!  Select numerical Method
   write(*,'(a)')'Select Method:'
   write(*,'(a)')' [1] Eulers Method'
   write(*,'(a)')' [2] RK Method'
   write(*,'(a)',advance='no')'Input :> '
   read(*,*) method

!  Check if selection is valid
   if(.not.(method .eq. 1 .or. method .eq. 2)) then
      write(*,*) 'Unidentified input, stopping program...'
      stop
   endif

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

!  Open the output file
   write(*,'(/,a)',advance='no')'Enter the output file name [velocity.dat] :> '
   read(*,"(a)") fname
   if( fname .eq. " ") fname = "velocity.dat"
   open(1,file=fname,form="formatted")

!  Write inital values to output file
   write(1,"(3f12.3)") time, velocity

!  Use selected method
   Select case(method)
   case(1) ! EULER METHOD 
      write(*,*) 'Using Eulers Method'
      do while ( lift .lt. takeoffweight )                               
         time  = time + stepsize                                ! New time is defined
         velocity = velocity + ODE(velocity)*stepsize           ! New velocity is defined
         lift = Liftf(velocity)                                 ! Calculate next lift to check if lift off occured
         write(1,"(3f12.3)") time, velocity
      enddo
      write(*,'(a,f6.3)') 'Minimum time needed to take off:', time
   case(2) ! RK2 METHOD  
      write(*,*) 'Using RK Method'
      do while ( lift .lt. takeoffweight )                             
         time = time + stepsize                                 ! New time is defined
         k1 = ODE(velocity)                                     ! k1 is calculated which is slope at previous velocity
         velocityx = velocity + k1*p1*stepsize                  ! velocityx is calculated
         k2 = ODE(velocityx)                                    ! k2 is calculated which is slope at velocityx
         velocity = velocity + ( a1*k1 + a2*k2 )* stepsize      ! New velocity is defined
         lift = Liftf(velocity)                                 ! Calculate next lift to check if lift off occured
         write(1,"(3f12.3)") time, velocity
      enddo
      write(*,'(a,f6.3)') 'Minimum time needed to take off:', time
   End select

!..Close the output file
   close(1)

   stop
END PROGRAM HOMEWORK1
