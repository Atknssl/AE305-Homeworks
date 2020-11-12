!--------------------------------------------------
!..Homework 1 
!..AE305 - Numerical Methods
!--------------------------------------------------
Module data
   implicit none
   real,parameter:: WParea=29.24, takeoffweight=88250, TmaxSL=16256, nrofengines=2, CD=0.2150, CL=1.792, friccoeff=0.02, airdensitySL=1.225, airdensity1000m=1.112, airdensity2000m=1.007, grav=9.81
   contains
     Function Lift(vel)
      Lift = CL*(1/2)*rho*(vel**2)*WParea
      return
      End Function Lift
   
     Function Drag(vel)
      Drag = CD*(1/2)*rho*(vel**2)*WParea
      return
      End Function Drag

     Function Thrust(dens)
      Thrust = TmaxSL*(dens/airdensitySL)
      return
      End Function Thrust
End module

PROGRAM HOMEWORK1
   use data

   implicit none

   character*40 :: fname
   real :: stepsize, k1, k2, a1, a2, p1,velocity, velocityx, weight, lift, time, thrust, rho
   integer :: method, selectedaltitude


!  Select numerical Method
   write(*,'(a)')'Select Method:'
   write(*,'(a)')' [1] Eulers Method'
   write(*,'(a)')' [2] RK Method'
   write(*,'(a)',advance='no')'Input :> '
   read(*,*) method

!  Check if selection is valid
!  Commented out for now, i dont know if its necessary

!   if(.not.(method .eq. 1 .or. method .eq. 2)) then
!      write(*,*) 'Unidentified input, stopping program...'
!      stop
!   endif

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
         write(*,*) 'Selected Sea Level' !This will be deleted later
         rho = airdensitySL
         thrust = TmaxSL*nrofengines
         ! Set density and thrust here 
      case(2)
         write(*,*) 'Selected 1000 m' !This will be deleted later
         rho = airdensity1000m
         thrust = Thrust(rho)*nrofengines
         ! Set density and thrust here
      case(3)
         write(*,*) 'Selected 2000 m' !This will be deleted later
         rho = airdensity2000m
         thrust = Thrust(rho)*nrofengines
         ! Set density and thrust here
      case default
         ! Cheap validation check
         write(*,*) 'Unidentified input, stopping program...'
         stop
   end Select

!  Select time step
   write(*,'(/,(a))',advance='no')'Enter Time Step :> '
   read(*,*) stepsize

!..open the output file
   write(*,'(/,a)',advance='no')'Enter the output file name [velocity.dat] :> '
   read(*,"(a)") fname
   if( fname .eq. " ") fname = "velocity.dat"
   open(1,file=fname,form="formatted")

!  Use selected method
   Select case(method)
   case(1)
      do while ( weight .lt. lift )                             ! EULER METHOD   *tanımlanmamış degişkenler*
         time  = time + stepsize                                ! New time is defined
         velocity = velocity + ODE(velocity)*stepsize           ! New velocity is defined by Old Velocity + slope of old velocity*step size  *tanımlanmamış fonksiyon*
         lift = Lift(velocity)
         write(1,"(3f12.3)") time, velocity
      enddo
      write(*,*) 'Using Eulers Method' !This will be deleted later
   case(2)
      do while ( weight .lt. lift )                             ! RK2 METHOD  
         time = time + stepsize                                 ! New time is defined
         k1 = ODE(velocity)                                     ! k1 is calculated which is slope at old velocity *tanımlanmamış fonksiyon*
         velocityx = velocity + k1*p1*stepsize                  ! velocityx is calculated   *tanımlanmamış değişken* 
         k2 = ODE(velocityx)                                    ! k1 is calculated which is slope at velocityx  *tanımlanmamış fonksiyon*
         velocity = velocity + ( a1*k1 + a2*k2 )* stepsize      ! New velocity is defined by Old Velocity + average slope of old velocity and velocityx * step size  *tanımlanmamış değişken a1, a2*
         lift = Lift(velocity)                                   ! *tanımlanmamış değişken ve fonksiyon*     
      enddo
      write(*,*) 'Using RK Method' !This will be deleted later
   End select

!..Close the output file
   close(1)

   stop
END PROGRAM HOMEWORK1

Function ODE(vel)
   use data
   ODE  = (thrust-Drag(vel)-friccoeff*(takeoffweight-Lift(vel)))*(grav/takeoffweight)
   return
 End