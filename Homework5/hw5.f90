Module vars
  integer, parameter :: imax=201, jmax=121
  integer, parameter :: kres=10, kqout=100
  integer :: kmax,irs,ire,jrs,jre
  real, parameter :: xl=10., yl=6., ck = 0.02, omega=1.2, rl2allow=-5.
  real, dimension(imax,jmax) :: q_k,q_kp1, qx=0., qy=0.
  real :: dx,dy,beta2,dx2i,dy2i, x(imax),y(jmax)
End module

!------------------------------------------------------------------------------|
!..A POINT/LINE ITERATIVE SOLVER FOR ELLIPTIC PDEs                             |
!  Course:  AE305                                                              |
!------------------------------------------------------------------------------|
 program ELLIPTIC
 use vars
!..Read the input data, generate the grid data and initialize the solution
  call INIT()
!..Start the iterative solution loop 
  k = 1
  DO WHILE (k .lt. kmax .and. rl2 .gt. rl2allow )
     k = k+1
!..Apply BCs 
     call BC()
!..Point iterative solutions
     call POINT_ITERATE()
!..or Line iterative solutions
!    call LINE_ITERATE()
!..Evaluate the L2 norm of the residual and update q_k
     rl2 = SQRT(SUM((q_kp1-q_k)**2))
     if(k .eq. 2) rl2_1=rl2
     rl2 = ALOG10(rl2/rl2_1)
     if( MOD(k,kres).eq.0 .or. k .eq. 2) then
       print*, 'Residual @ k =',k,rl2
       write(2,*) k, rl2
     endif
     q_k = q_kp1
!..Output intermediate solutions
     if( MOD(k,kqout) .eq. 0 .and. k .ne. kmax) call QOUT(k)
  ENDDO
  print*, 'Residual @ k =',k,rl2
  write(2,*) k, rl2
  close(2)
  call QOUT(k)

 stop
end

!------------------------------------------------------------------------
subroutine INIT()
 use vars

  write(*,'(a)',advance='no') 'Enter kmax: '
  read(*,*) kmax

  dx = xl/(imax-1)
  dy = yl/(jmax-1)
  beta2 = (dx/dy)**2

!..Grid generation 
  do i=1,imax
     x(i)= dx*(i-1) 
  enddo
  do j=1,jmax
     y(j)= dy*(j-1) 
  enddo
  q_k = 0.       !..Initial guess

!..Set radiator location
  irs = 2/dx + 1
  ire = 8/dx + 1
  jrs = 5.6/dy + 1
  jre = 5.8/dy + 1

 open(2,file='residual.dat',form='formatted')
 return 
end

!-------------------------------------------------------------------
subroutine POINT_ITERATE()
 use vars

  cm = 1.0/(2.0*(1. + beta2)) 
  do j = 2,jmax-1
  do i = 2,imax-1
!..Exclude the radiator from the solution domain
!...
!..Implement, Point Jacobi, Gauss-Seidel and SOR methods
   q_kp1(i,j) = cm*( q_k(i-1,j) + q_k(i+1,j) + beta2*(q_k(i,j-1) + q_k(i,j+1)) )
  enddo
  enddo

 return 
end

!-------------------------------------------------------------------
subroutine BC()
 use vars
  !..Set the bottom/top farfield BC
  q_k(:,1)      = 20.
  q_k(:,jmax)   = 0.
  q_kp1(:,1)    = q_k(:,1)
  q_kp1(:,jmax) = q_k(:,jmax)

  !..Set left/right farfield BC
  q_k  (1,:)    = q_k(2,:)
  q_k  (imax,:) = 20.
  q_kp1(1,:)    = q_k(1,:)
  q_kp1(imax,:) = q_k(imax,:)

!..Set BCs for the radiator
! T_k( irs:ire, jrs:jre )   = 50
! T_kp1( irs:ire, jrs:jre ) = 50

 return 
end

!-------------------------------------------------------------------
subroutine FLUX()
  use vars
!..Compute the heat flux vectors
  do j = 2,jmax-1
  do i = 2,imax-1
       qx(i,j) = ..
       qy(i,j) = ..
  enddo
  enddo
 return 
end

!-------------------------------------------------------------------
subroutine QOUT(k) 
 use vars
   character fname*32,string*9,ext*6
   write(string,'(f9.6)') float(k)/1000000
   read(string,'(3x,a6)') ext
   fname = 'q-'//ext//'.tec'
   open(1,file=fname,form='formatted')
   call FLUX
   write(1,*) ' variables="x","y","T","Tx","Ty"'
   write(1,*) ' zone i=',imax, ', j=',jmax
   do j=1,jmax
   do i=1,imax
      write(1,*) x(i),y(j),q_k(i,j),qx(i,j),qy(i,j)
   enddo
   enddo
   close(1)
 return
end

!-------------------------------------------------------------------
subroutine THOMAS(is,ie, a,b,c,f)
  use vars, only : imax
  real, dimension(imax) ::  a,b,c,f,x
!
!  Solution of a tridiagonal system of n equations of the form
!  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=is,ie
!  the solution X(i) is stored in F(i)
!  A(is-1) and C(ie+1) are not used.
!  Solution is returned in array F
!
   x(is)=c(is)/b(is)
   f(is)=f(is)/b(is)
   isp1 = is+1
   do i=isp1,ie
      z   =1./(b(i)-a(i)*x(i-1))
      x(i)=z*c(i)
      f(i)=(f(i)-a(i)*f(i-1))*z
   enddo
   iepis=ie+is
   do ii=isp1,ie
      i=iepis-ii
      f(i)=f(i)-x(i)*f(i+1)
   enddo
     
 return
end

!-------------------------------------------------------------------
subroutine GAUSS(n,a,b)
   real a(n,n),b(n)
!..Convert to upper triangular form
   do k = 1,N-1
   if (ABS(a(k,k)).gt.1.E-6) THEN
      do i = k+1, n
      x = a(i,k)/a(k,k)
         do j = k+1, n
            a(i,j) = a(i,j) -a(k,j)*x
         enddo
      b(i) = b(i) - b(k)*x
      enddo
   else
      write (6,*) 'zero pivot found in line:', k
      stop
   endif
   enddo
!..Back substitution
   do i = n,1,-1
     sum = b(i)
     if (i.lt.n) then
       do j= i+1,n
         sum = sum - a(i,j)*b(j)
       enddo
     endif
     b(i) = sum/a(i,i)
   enddo
  return
end

