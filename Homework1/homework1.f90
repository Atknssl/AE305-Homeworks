!--------------------------------------------------
!..Homework 1 
!..AE305 - Numerical Methods
!--------------------------------------------------
Module data
   implicit none

End module

PROGRAM HOMEWORK1
   use data

   implicit none

   character*40 :: fname
   real :: stepsize
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
         ! Set density and thrust here 
      case(2)
         write(*,*) 'Selected 1000 m' !This will be deleted later
         ! Set density and thrust here
      case(3)
         write(*,*) 'Selected 2000 m' !This will be deleted later
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
      !Eulers Method code here
      write(*,*) 'Using Eulers Method' !This will be deleted later
   case(2)
      !RK Method code here
      write(*,*) 'Using RK Method' !This will be deleted later
   End select

!..Close the output file
   close(1)

   stop
END PROGRAM HOMEWORK1
