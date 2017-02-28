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
type (artific) art

open(1,file = 'mesh.inp')

read(1,*) n_node !number nodes 
allocate(node(n_node))
do i = 1,n_node
	read(1,*) node(i)%z,node(i)%r !read nodes 
enddo
read(1,*) n_bound !number bound
allocate(bou%var(n_bound,7))
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
	if (el(i)%elem(1) == 9) then !Number 9 - see in GMSH
		phy(i)%rho = 1d0
		phy(i)%e = 2.5d0
	endif
	if (el(i)%elem(1) == 11) then !Number 11 - see in GMSH
		phy(i)%rho = 0.125d0
		phy(i)%e = 2.0d0
	endif	
enddo	
n = 6
allocate(bou%type_bound(n))


bou%type_bound(1) = 3
bou%type_bound(2) = 3
bou%type_bound(3) = 1
bou%type_bound(4) = 3
bou%type_bound(5) = 3
bou%type_bound(6) = 1




art%grid = 0d0 !Set grid velocity 0:1

node(:)%u = 0d0 !
node(:)%v = 0d0 !


dt = 0.0001 !time interval
t = 0d0 !inintial time
t_end = 0.25d0 !End time calculation

call volume0(el,node,phy) !Initial 
call gas(phy) !initial 

do while(t<t_end)
call volume(el,node,phy)
call artvisc(dt,node,phy)
call velocity(dt,el,node,phy)
call boundary_flow(bou,node)
call energy(dt,el,node,phy)
call boundary_flow(bou,node)
call grid(dt,node,art)
call advect(dt,el,node,phy,art)
call boundary_flow(bou,node)
call gas(phy)
t = t + dt
enddo

call out(t,el,node,phy)
end
