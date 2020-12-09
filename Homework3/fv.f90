!------------------------------------------------------------------------------|
!..A 2-D FINITE VOLUME SOLVER FOR THE DIFFUSION EQUATION                       |
!                      ON UNSTRUCTURED GRIDS                                   |
!  Course: AE305                                                               |
!------------------------------------------------------------------------------|
Module vars
  integer :: ncell, nnode
  integer,parameter :: ntmax=100000, ntio=100
  integer,allocatable,dimension(:,:) :: node, neigh
  real,parameter :: dt=1.E-3, delq_allow = 1.E-6, pi = 4*ATAN(1.0)
  real,allocatable,dimension(:,:) :: xy, qgrad
  real,allocatable,dimension(:) :: qcell,qnode,qbc,area,unode,vnode
  real :: uinf, vinf, alpha_deg = 15.0
End module

 program FiniteVolume
  Use vars
!..Read the input data and initialize the solution
    call INIT()
!..Start the solution loop 
    DO nt=1,ntmax
       call GRADIENT      !..Evaluate Q gradients for cells
       delq_max = 0.
       do n = 1,ncell     !..Sweep all the cells and solve for T^n+1
          delq     = -dt/area(n) * CellFLUX(n)
          qcell(n) = qcell(n) + delq
          delq_max = MAX(ABS(delq),delq_max)
       enddo
!..Output the intermediate solutions 
       if(delq_max .lt. delq_allow) exit
       if( MOD(nt,ntio) .eq. 0 .and. nt .ne. ntmax ) then
		  print*, ' nt, Delq_max :',nt,delq_max
	      call TECout(nt)
	   end if
    ENDDO                        
    call TECout(nt-1)       !..Output the final solution
  stop  'DONE'
 end

!------------------------------------------------------------------------
 subroutine INIT()
  Use vars 
  character :: fn*16
  logical   :: ok

!..Read the grid data
   write(*,'(/(a))',advance='no')'  Enter the grid file name [circle.dat]: '
   read(*,'(a)') fn
   if( fn .eq. ' ') fn = 'circle.dat'
   inquire(FILE=fn,EXIST=ok)
   if( .not. ok ) then
       print*, '  ', fn, ' does not exist! \n\n'
       stop
   endif 
   open(1,file=fn,form='formatted')
   read(1,*) ncell,nnode
   allocate( node(3,ncell),neigh(3,ncell),xy(2,nnode),area(ncell), &
             qcell(ncell),qnode(nnode),unode(nnode),vnode(nnode), &
             qbc(4),qgrad(2,ncell) )
   read(1,*) (no,(xy(i,n),i=1,2),n=1,nnode)
   read(1,*) (no,(node(i,n),i=1,3),(neigh(i,n),i=1,3),n=1,ncell) 
   close(1)
   print*, ' # of cells :',ncell
   print*, ' # of nodes :',nnode

!..Compute cell areas
   do n = 1,ncell
      n1 = node(1,n)
      n2 = node(2,n)
      n3 = node(3,n)
      area(n) = 0.5*( (xy(1,n2)-xy(1,n1))*(xy(2,n3)-xy(2,n1)) &
                     -(xy(2,n2)-xy(2,n1))*(xy(1,n3)-xy(1,n1))  )
   enddo

	alpha = alpha_deg/180.*pi
	uinf  = COS(alpha)
	vinf  = SIN(alpha)

!..Initialize the solution with 0
	qcell = 0.

   call TECout(0)

   return 
 end

!-------------------------------------------------------------------
subroutine  GRADIENT()
  Use vars
  qgrad = 0.
  DO n = 1,ncell
     do ns = 1,3   
        n1 = node(ns,n)
        ni2 = MOD(ns,3) + 1
        n2 = node(ni2,n)
        dx = xy(1,n2)-xy(1,n1)
        dy = xy(2,n2)-xy(2,n1)
        ne = neigh(ns,n)
        if(ne .gt. 0) then       !..real neighbor
           qface = 0.5*(qcell(n) + qcell(ne))
           qgrad(1,n) = qgrad(1,n) + (1/area(n)) * qface * dy
           qgrad(2,n) = qgrad(2,n) - (1/area(n)) * qface * dx
        elseif(ne.eq.-1) then
            qgrad(1,n) = uinf
            qgrad(2,n) = vinf
            exit
        elseif(ne.eq.-2) then
         qface = qcell(n)
         qgrad(1,n) = qgrad(1,n) + (1/area(n)) * qface * dy
         qgrad(2,n) = qgrad(2,n) - (1/area(n)) * qface * dx
        endif
     enddo 
  ENDDO

  return 
 end

!------------------------------------------------------------------------
function CellFLUX(n)
  Use vars

  Cellflux = 0.
  do ns = 1,3               !..Add the surface fluxes
     n1 = node(ns,n)
     ni2 = MOD(ns,3) + 1
     n2 = node(ni2,n)
     dx = xy(1,n2)-xy(1,n1)
     dy = xy(2,n2)-xy(2,n1)
!..Apply proper BCs when computing the face fluxes
     ne = neigh(ns,n)
     if( ne .gt. 0 ) then        !..real neighbor... 
         f = ( qgrad(1,n) + qgrad(1,ne) ) / 2.
         g = ( qgrad(2,n) + qgrad(2,ne) ) / 2.
         elseif(ne.eq.-1) then
            f = uinf
            g = vinf
        elseif(ne.eq.-2) then
            f = 0
            g = 0
        endif
     Cellflux = Cellflux - (f*dy - g*dx)
  enddo 
  return
end

!-------------------------------------------------------------------
subroutine  TECout(nstep) !..Output the solution/grid in TECPLOT format
  Use vars 
  character :: fname*32, str*12, ext*9
   call QNODES()
!..Set the output file name
   write(str,'(f12.9)') float(nstep)/1e9
   read(str,'(3x,a9)') ext
   fname = 'q-'//ext//'.tec'
   open(1,file=fname, form='formatted')
   write(1,100) nnode,ncell
   write(1,101) (xy(1,n),xy(2,n),qnode(n),unode(n),vnode(n),n=1,nnode)
   write(1,102) (node(1,n),node(2,n),node(3,n),n=1,ncell)
   close(1)
  100 format (' VARIABLES= "X", "Y", "Q", "U", "V"'/, &
              ' ZONE N=', I6,' E=', I6,' F=FEPOINT ',' ET=triangle'  )
  101 format (5(1x,e12.5))
  102 format (3(1x,i6))
  return 
end
!-------------------------------------------------------------------
subroutine  QNODES()     !..Evaluate averaged q_node values
  Use vars 
  integer, dimension(:) :: npass(nnode) 
   qnode = 0.
   unode = 0.
   vnode = 0.
   npass = 0.
!..Find the contribution of cells to the node Q
   do n=1,ncell
   do nf=1,3
      nn = node(nf,n)
      qnode(nn)=qnode(nn)+qcell(n)
      unode(nn)=unode(nn)+qgrad(1,n)
      vnode(nn)=vnode(nn)+qgrad(2,n)
      npass(nn)=npass(nn)+1
   enddo
   enddo
!..Average the total node Q with # of contributing cells 
   do n=1,nnode
      qnode(n)=qnode(n)/npass(n)
      unode(n)=unode(n)/npass(n)
      vnode(n)=vnode(n)/npass(n)
   enddo
  return 
end
!------------------------------END----------------------------------
