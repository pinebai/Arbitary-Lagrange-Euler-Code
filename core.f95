program ale
use object
use eos
use out_paraview
use hydro
use type_boundary
implicit none 

integer(4) i,j,k,l,n 
integer(4) slide,kl
integer(4) n_bound,n_node,n_cell,num
integer(4),allocatable :: bound(:)
character(20) name_input,type_init,type_area
!Var
real(8) metric,t,t_end,dt
real(8) z0,z,r0,r,rad,zc,rc
real(8) ener,dens
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

write(*,*) 'Write name file (see input.inp)'
read(*,*) name_input
write(*,*) name_input
open(1,file = name_input)
read(1,*)
read(1,*)
read(1,*) metric
do i = 1,n_node
    node(i)%z = metric*node(i)%z
    node(i)%r = metric*node(i)%r 
enddo
!Start time, time interval, end time, number out file
read(1,*) 
read(1,*) t,dt,t_end,slide
read(1,*)
!Radial, Grid proportional velocity, Artification viscosity, a0 for Euler stability
read(1,*) rad, numer%grid, numer%art, numer%a0
read(1,*)
read(1,*) n,l !all line in mesh
el(:)%rad = rad
allocate(bou%type_bound(n))
n = l
bou%type_bound = 0
do i = 1,n
    read(1,*) j,l
    bou%type_bound(j) = l
enddo
!Initial
read(1,*) 

read(1,*) type_init
if (type_init == 'gmsh') then
    read(1,*) n !all regions
    do j = 1,n
        read(1,*) num,dens,ener
        do i = 1,n_cell
            if (el(i)%elem(5) == num) then !Number mat - see in GMSH
                phy(i)%rho = dens
                phy(i)%e = ener 
            endif
        enddo
    enddo
endif
if (type_init == 'master') then
    read(1,*) n !all regions
    do j = 1,n
        read(1,*) type_area 
        if (type_area == 'plane') then
            read(1,*) z0,r0,z,r
            read(1,*) dens,ener
            do i = 1,n_cell
                zc = 0.25d0*sum(node(el(i)%elem(:4))%z)
                rc = 0.25d0*sum(node(el(i)%elem(:4))%r)
                if ((zc>=z0).and.(rc>=r0).and.(zc<=z).and.(rc<=r)) then
                    phy(i)%rho = dens
                    phy(i)%e = ener 
                endif    
            enddo
        endif    
        if (type_area == 'circle') then
            read(1,*) z0,r0,Rad
            read(1,*) dens,ener        
            do i = 1,n_cell
                zc = 0.25d0*sum(node(el(i)%elem(:4))%z)
                rc = 0.25d0*sum(node(el(i)%elem(:4))%r)
                if ((zc-z0)**2+(rc-r0)**2<=rad**2) then
                    phy(i)%rho = dens
                    phy(i)%e = ener 
                endif
            enddo
        endif    
    enddo
endif
close(1)
kl = 0
call phase0(el,node,phy,bou,numer) !Initial 

call system('rm -rf result')
call system('mkdir result')

do while(t<=t_end+dt)
    call out(kl,slide,t,el,node,phy)
	call phase1(dt,el,node,phy,numer)
	call velocity(dt,el,node,phy,bou,numer)
	call energy(dt,el,node,phy,numer)
	call grid(dt,node,numer)
	call advect(dt,el,node,phy,bou,numer)
	t = t + dt
	kl = kl+1
	write(*,*) 'cycle =',kl, 'time', t 
enddo
end
