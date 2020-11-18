!------------------------------------------------------------------------
!  RK4 SOLVER for a system of ODEs                                      |
!  Course:  AE305                                                       |
!------------------------------------------------------------------------
module data
	integer, parameter :: neq=2
	real*8, parameter :: rho_p=1140., a=5.55e-5 ! ..., ..., ...
end module
	
program sysRK4
	use data
	character*40 :: fname

!..read the input data
	write(*,'(/(a))',advance='no')' Enter the step size and the end point :>'
	read(*,*) dt,tmax
!..open the output file 
	write(*,'(a)',advance='no')' Enter the output file name [solution.dat]:>'
	read(*,'(a)') fname
	if( fname .eq. ' ') fname = 'solution.dat'
	open(1,file=fname,form='formatted')

!..set the initial conditions
	time = 0.
!~ 	y(1) = ...
!~ 	y(2) = ...
	write(1,'(6E14.5)') time,(y(i),i=1,neq)

!..solution loop
	DO WHILE (time .lt. tmax)      
		call SRK4(time,dt,y)
		time = time + dt 
		write(1,'(6E14.5)') time,(y(i),i=1,neq)
	ENDDO

	close(1)
	stop
end
	
!------------------------------------------------------------------------
subroutine SRK4(time,dt,y)
	real*8 :: y(neq),ytmp(neq),k1(neq),k2(neq),k3(neq),k4(neq)
	dt2 = 0.5*dt
		
	call odes(time,y,k1)
	do i = 1,neq
         ytmp(i)  = y(i) + k1(i)*dt2
	enddo
!~ 	call ODES( ... )
!~ 	...
!~ 	...

!..obtain the solution at t+dt and update y for the next step
	do i = 1,neq
		phi  = (k1(i) + 2*(k2(i)+k3(i)) + k4(i))/6.
		y(i) = y(i) + phi*dt
	enddo

	return
	end

!------------------------------------------------------------------------
subroutine ODEs(time,y,f)
      use data
	real*8 :: f(neq)

!..define the ODE's & return the slopes in the "f" array

!~ 	f(1) = ...
!~ 	f(2) = ...

	return
end
