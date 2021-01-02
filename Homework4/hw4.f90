!-------------------------------------------------------------------
!..AN EXPLICIT FD SOLVER FOR THE 1-D LINEAR CONVECTION EQUATION  
!-------------------------------------------------------------------
Module vars
   integer, parameter :: imax=401, ntout=1 
   integer :: ntmax, method
   real, dimension(imax) :: wave_n, wave_np1, wave_nm1
   real :: sigma, dx=0.1, x1 , pi= 3.1415926 !ACOS(-1.)
   real :: d, b=0.025
 End module
 
 program EXPLICIT_FDE
    use vars
 
    call INIT             !..Read the input data, and initialize the wave
    DO nt = 1,ntmax       !..Start the solution loop 
       do i = 2,imax-1
          if(method .eq. 1) then
          wave_np1(i) = wave_n(i) - 0.5*sigma*(wave_n(i+1) - wave_n(i-1)) + d *(wave_n(i+1) - 2 * wave_n(i) + wave_n(i-1))
          elseif(method .eq. 2) then
                if(i.eq.2) then
                   wave_np1(i) = wave_n(i) - 0.5*sigma*(wave_n(i+1) - wave_n(i-1)) + d *(wave_n(i+1) - 2 * wave_n(i) + wave_n(i-1))
                else
                   wave_np1(i) = wave_n(i) - sigma*(wave_n(i) - wave_n(i-1)) + d *(wave_n(i) - 2 * wave_n(i-1) + wave_n(i-2))
                endif
          elseif(method .eq. 3) then
                if(nt .eq. 1) then
                   wave_np1(i) = wave_n(i) - 0.5*sigma*(wave_n(i+1) - wave_n(i-1)) + d *(wave_n(i+1) - 2 * wave_n(i) + wave_n(i-1))
                else
                   wave_np1(i) = 4*wave_n(i) - 3*wave_nm1(i) + sigma*(wave_n(i+1) - wave_n(i-1)) &
                   - 2*d*(wave_n(i+1) - 2*wave_n(i) + wave_n(i-1))
                endif
          endif
       enddo
       wave_nm1(:) = wave_n(:)
       wave_n(:) = wave_np1(:)
 !..Output intermediate solutions
       if( MOD(nt,ntout) .eq. 0 .or. nt .eq. ntmax ) call IO(nt)
    ENDDO                        
 
    stop
 end program EXPLICIT_FDE
     
 !------------------------------------------------------------------------
 subroutine INIT
   use vars
   write(*,'(a)')'Select method'
   write(*,'(a)')' [1] Forward time central spatial'
   write(*,'(a)')' [2] Forward time backward spatial'
   write(*,'(a)')' [3] Question 4'
   read(*,*) method
   if((method .ne. 1) .and. (method .ne. 2) .and. (method .ne. 3)) then
    write(*,'(a)') "Invalid selection"
    stop
   endif
   
   write(*,'(/(a))',advance='no')'  Enter sigma, d and ntmax : '
   read(*,*) sigma, d, ntmax
   x1=-imax*dx/2.
   x = x1
   do i = 1,imax             !..Initialize the wave 
      wave_n(i) = exp(-b * log(2.) * (x/dx)**2)
      x = x+dx
   enddo
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
 
 