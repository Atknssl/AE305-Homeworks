!-------------------------------------------------------------------
!..AN EXPLICIT FD SOLVER FOR THE 1-D LINEAR CONVECTION EQUATION  
!-------------------------------------------------------------------
Module vars
  integer, parameter :: imax=201, ntout=1 
  integer :: ntmax
  real, dimension(imax) :: wave_n, wave_np1, f, wave_n0
  real, dimension(imax) :: a, b1, c 
  real :: sigma, d, dx=0.1, x1 , b=0.025, pi= 3.1415926 !ACOS(-1.)
End module

program IMPLICIT_FDE
   use vars

   call INIT             !..Read the input data, and initialize the wave
   DO nt = 1,ntmax       !..Start the solution loop 

      do i = 2, imax-1
         a(i) = -(sigma/2.)-d
         b1(i) = 2.*d + 1. 
         c(i) = (sigma/2.)-d
         if (i.eq. 2) then
            a(i)=0
         elseif (i.eq.imax) then
            c(i)=0
         endif
            
      enddo
      if (nt.eq.1) then
      wave_n(:)=wave_n0(:)
      endif
      print*, wave_n0(3)
      call THOMAS(imax, 3, imax-2, a, b1, c, wave_n)
      wave_n(:)= f(:)
      
      
      
!..Output intermediate solutions
      if( MOD(nt,ntout) .eq. 0 .or. nt .eq. ntmax ) call IO(nt)
   ENDDO                        

   stop
end program IMPLICIT_FDE
    
!------------------------------------------------------------------------
subroutine INIT
  use vars

  write(*,'(/(a))',advance='no')'  Enter sigma, d and ntmax : '
  read(*,*) sigma, d,  ntmax
  x1=-imax*dx/2.
   x = x1
   do i = 1,imax             !..Initialize the wave 
      wave_n(i) = exp(-b * log(2.) * (x/dx)**2)
      x = x+dx
   enddo
   wave_n0(:)=wave_n(:)
  call IO(0)

  return 
end subroutine INIT

!-------------------------------------------------------------------
subroutine  IO(nt)
   use vars
   character :: fname*32,string*6,ext*3
   write(string,'(f5.3)') float(nt)/1000
   read(string,'(2x,a3)') ext
   fname = 'wave-'//ext//'.dat' 
   open(1,file=fname,form='formatted')
   !write(1,*)'ZONE'
   x = x1
   do i=1,imax
     write(1,'(2e14.6)')x, wave_n(i)
     x = x+dx
   enddo
   close(1)
   return 
end
!-------------------------------------------------------------------
subroutine  THOMAS(imax, il, iu, a, b, c, f)
   integer ::  imax, il ,iu, ilp1, iupil
   real, dimension(imax) :: a, b, c, f, x
   real :: z
   x(il) = c(il) / b(il)
   f(il) = f(il) / b(il)
   print*, f(il)
   ilp1 = il+1
   do i = ilp1, iu
      z= 1. / (b(i)- a(i)*x(i-1))
      x(i)= c(i) * z 
      f(i)= ( f(i) - a(i) * f(i-1) ) * z
      
   ENDDO
   iupil = iu + il
   do ii = ilp1, iu
      i = iupil - ii
      f(i) = f(i) - x(i) * f(i+1)
   ENDDO
   
   return
   end