!--------------------------------------------------
!..A SIMPLE EULER SOLVER for ODEs
!..AE305 - Numerical Methods
!--------------------------------------------------
Module data
  parameter ( grav = 9.81 )
  parameter ( fric = 12.5 )
  parameter ( xmass = 70. )
End module

Program EULER
  character*40 fname

!..Read the stepsize
   print*," "
   write(*,'(/,(a))',advance='no')'  Enter TimeStep and FinalTime :> '
   read(*,*) stepsize, finaltime

!..open the output file
   write(*,'(a)',advance='no')'  Enter the output file name [velocity.dat]: '
   read(*,"(a)") fname
   if( fname .eq. " ") fname = "velocity.dat"
   open(1,file=fname,form="formatted")

!..Set the Initial Conditions and output them
   time      = 0.
   velocity  = 0.
   write(1,"(3f12.3)") time, velocity


!..Solution loop
   do while ( time .lt. finaltime )
      time  = time + stepsize
      velocity = velocity + ODE(velocity)*stepsize
      write(1,"(3f12.3)") time, velocity
   enddo

!..Close the output file
   close(1)

!..Plot the solution with "xgraph", "xmgr" or "gnuplot" on Linux X-windows
!     call SYSTEM("xg velocity.dat &")
!     call SYSTEM("xmgr velocity.dat &")
!     call SYSTEM("gnuplot plot.gpl")

   stop
End


!--------------------------------------------------
!..Define the ODE as a function
Function ODE(vel)
  use data
  ODE  = grav - fric/xmass * vel
  return
End
