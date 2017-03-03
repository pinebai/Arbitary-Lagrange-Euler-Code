program ale
use object
use eos
use out_paraview
use hydro
use type_boundary
implicit none 

integer(4) i,j,k,l,n 
integer(4) n_bound,n_node,n_cell

integer(4),allocatable :: bound(:)
!Var
real(8) t,t_end,dt
type (physics),allocatable :: phy(:) !All physical value
type (nodes),allocatable :: node(:)	!All nodes value
type (elements),allocatable :: el(:) !Elements (cell {node1,node2,node3,node4}) and contact elements for Advection
type (boundary) bou !Boundary elements
type (numerical) numer

open(1,file = 'mesh.inp')

read(1,*) n_node !number nodes 
allocate(node(n_node))
do i = 1,n_node
	read(1,*) node(i)%z,node(i)%r !read nodes 
enddo
read(1,*) n_bound !number bound
allocate(bou%var(n_bound,6))
do i = 1,n_bound
	read(1,*) bou%var(i,:) !read all elements
enddo

read(1,*) n_cell !numbe all cells
allocate(el(n_cell))  
do i = 1,n_cell
	read(1,*) el(i)%elem(:) !read elements
enddo

read(1,*)
do i = 1,n_cell
	read(1,*) el(i)%cont(:) !read all contacts
enddo
close(1)

allocate(phy(n_cell)) !allocate physical
!Initial conditions
do i = 1,n_cell
	if (el(i)%elem(1) == 18) then !Number mat - see in GMSH
		phy(i)%rho =  1.0d0
		phy(i)%e = 2.5d0
	endif
	if (el(i)%elem(1) == 14) then 
		phy(i)%rho = 0.125d0
		phy(i)%e = 2.0d0
	endif	
	if (el(i)%elem(1) == 16) then 
		phy(i)%rho = 0.125d0
		phy(i)%e = 2.0d0
	endif	
	if (el(i)%elem(1) == 22) then 
		phy(i)%rho = 0.125d0
		phy(i)%e = 2.0d0
	endif	
	if (el(i)%elem(1) == 20) then 
		phy(i)%rho = 0.125d0
		phy(i)%e = 2.0d0
	endif	
enddo	
numer%art = 0.1d0 !Artification viscosity
numer%grid = 1.0d0 !Grid proportional velocity
numer%a0 = 1d0 ! a0 for Euler stability

el(:)%rad = 2d0 !Radial

n = 12
allocate(bou%type_bound(n))
bou%type_bound = 0
bou%type_bound(1) = 1 !Solid wall - fix all boundary velocity
bou%type_bound(2) = 1
bou%type_bound(3) = 1
bou%type_bound(4) = 5 !Reflection


dt = 0.0001 !time interval
t = 0d0 !inintial time
t_end = 0.2d0 !End time calculation

call phase0(el,node,phy,numer,bou) !Initial 
do while(t<t_end)
call phase1(dt,el,node,phy,numer)
call velocity(dt,el,node,phy,numer)
call boundary_flow(bou,node,phy,el)
call energy(dt,el,node,phy,numer)
call grid(dt,node,numer)
call advect(dt,el,node,phy,numer)
call boundary_flow(bou,node,phy,el)
t = t + dt
enddo

call out(t,el,node,phy)
end
