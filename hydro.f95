module hydro
use object
use eos
use type_boundary
implicit none
contains

subroutine phase0(el,node,cell,bou,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type(boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
integer(4) i,j,l(4)
real(8) atr,abl

node(:)%mas_v = 0d0
node(:)%mark = 0

do i = 1,size(cell(:)%vol)

! 3-----2
! |     |
! 4-----1
!
!Node equivalent 
!Sale
!1 2 3 4
!el(i)%elem
!1 2 3 4

	do j = 1,4
        l(j) = el(i)%elem(j) 
    enddo


!Calculation triangles left and right in cell 
	atr = 0.5d0*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
	abl = 0.5d0*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))
 
    if (el(i)%rad == 2) then
        atr = atr*((node(l(1))%r+node(l(2))%r+node(l(3))%r)/3d0)
        abl = abl*((node(l(1))%r+node(l(3))%r+node(l(4))%r)/3d0)
	endif

	
	cell(i)%sig = sign(1d0,abl+atr) !Sign volume (Very important)
	cell(i)%vol = abs(abl)+abs(atr) !Calculaion volume cell
	
	cell(i)%vol_old = cell(i)%vol !Set first old volume 
	cell(i)%mas = cell(i)%rho*cell(i)%vol ! Calculation massa in cell
	
	call gas(i,cell)

!Set old velocity
	node(l(1))%u_l = node(l(1))%u 
	node(l(2))%u_l = node(l(2))%u
	node(l(3))%u_l = node(l(3))%u
	node(l(4))%u_l = node(l(4))%u

	node(l(1))%v_l = node(l(1))%v
	node(l(2))%v_l = node(l(2))%v
	node(l(3))%v_l = node(l(3))%v
	node(l(4))%v_l = node(l(4))%v

!Calculation vertex mass 
	do j = 1,4 !cycle for all nodes in current cell 
		node(l(j))%mas_v = node(l(j))%mas_v + 0.25d0*cell(i)%mas
	enddo
enddo
call boundary_flow(bou,node,cell,el)
end subroutine phase0


subroutine phase1(dt,el,node,cell,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,l(4)
real(8) atr,abl,u(4),v(4),r(4)
real(8) Dvol

do i = 1,size(cell(:)%vol)
	!Calculation EOS
	call gas(i,cell)
	
	!Calculation volume,mass and density
	!first material then shift to 1
	do j = 1,4
        l(j) = el(i)%elem(j) 
        u(j) = node(l(j))%u  
        v(j) = node(l(j))%v
        r(j) = node(l(j))%r        
    enddo
    
	atr = 0.5d0*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
	abl = 0.5d0*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

    if (el(i)%rad == 2) then
        atr = atr*((node(l(1))%r+node(l(2))%r+node(l(3))%r)/3d0)
        abl = abl*((node(l(1))%r+node(l(3))%r+node(l(4))%r)/3d0)
	endif
	
	cell(i)%vol_old = cell(i)%vol !Set old volume 
	cell(i)%vol = abs(abl)+abs(atr) !Calcualtion new volume 
	cell(i)%rho = cell(i)%mas/cell(i)%vol !Calculation density
	cell(i)%mas = cell(i)%rho*cell(i)%vol !Calculation massa
	
	!Calculation artvisc
	Dvol = 2d0*(cell(i)%vol - cell(i)%vol_old)/(cell(i)%vol+cell(i)%vol_old)/dt !Calculation Divergence (Equivalen Hemp3D)
	cell(i)%q = min(0d0,Dvol)*(numer%art*cell(i)%rho*Dvol*cell(i)%Vol**0.666) !Calculation Divergence (Equivalen Hemp3D)	
	
	call difference(i,l,u,cell(i)%uz,cell(i)%ur,el,node)
	call difference(i,l,v,cell(i)%vz,cell(i)%vr,el,node)
	
	
	call divergence(cell(i)%uz,cell(i)%vr,r,v,el(i)%rad,cell(i)%diver)
    call rotor(cell(i)%ur,cell(i)%vz,r,v,el(i)%rad,cell(i)%rot)
	!rotor(dvec_zdr,dvec_rdz,rot)
enddo
end subroutine phase1


subroutine velocity(dt,el,node,cell,bou,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type (boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l,im,ip,ia
real(8) acc_z,acc_r
real(8) an,ksi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Relaxation grid
!This cycle for calculation realxation lagrange grid (See SALE)
!an and ksi - var parameters
!See more comments in call posit(...)
do i = 1,size(cell(:)) !
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		call posit(i,j,el,im,ip,ia) !Function calculation position + (ip) and - (im)
node(l)%u = node(l)%u + 0.25d0*numer%an*(0.5d0*(1d0+numer%ksi)*(node(ip)%u_l+node(im)%u_l)-numer%ksi*node(ia)%u_l-node(l)%u_l)
node(l)%v = node(l)%v + 0.25d0*numer%an*(0.5d0*(1d0+numer%ksi)*(node(ip)%v_l+node(im)%v_l)-numer%ksi*node(ia)%v_l-node(l)%v_l)
	enddo
enddo

call boundary_flow(bou,node,cell,el)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1,size(cell(:)) !
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		call posit(i,j,el,im,ip,ia) !Function calculation position 
		
		acc_z = 0.5d0*cell(i)%sig*(cell(i)%p+cell(i)%q)*(node(ip)%r-node(im)%r)/node(l)%mas_v !Caluclation z-component acceleration 
		acc_r = 0.5d0*cell(i)%sig*(cell(i)%p+cell(i)%q)*(node(im)%z-node(ip)%z)/node(l)%mas_v !Caluclation r-component acceleration
        
        if (el(i)%rad == 2) then
            acc_z = acc_z*(0.5d0*(node(ip)%r+node(im)%r)) !Radial component for z
            acc_r = acc_r*node(l)%r                       !Radial component for r
        endif
		
		node(l)%u = node(l)%u - acc_z*dt !Caluclation z-component velocity 
		node(l)%v = node(l)%v - acc_r*dt !Caluclation r-component velocity 
	enddo	
enddo	

call boundary_flow(bou,node,cell,el)
contains

subroutine posit(i,j,el,im,ip,ia)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: i,j
integer(4),intent(out) :: im,ip,ia
integer(4) l 

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This countur for node (o)
!im,ip - adjacent nodes
!ia - oposite node
!sig - direction rotation
!
! n-----ip----ia 
! |   / | ^   |
! |  /  |  \  |
! | v   |   \ |
! n----(o)----im 
! | \   |   ^ |
! |  \  |  /  |
! |   v | /   |
! n-----n-----n
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

select case(j)
	case(1) !Node 1
! 3-----2 
! |   / |
! |  /  |
! | v   |
! 4----(1)
		im = el(i)%elem(2)    
		ip = el(i)%elem(4) 
		ia = el(i)%elem(3)
	case(2) !Node 2
! 3----(2) 
! | \   |
! |  \  |
! |   v |
! 4-----1 
		im = el(i)%elem(3)
		ip = el(i)%elem(1)
		ia = el(i)%elem(4)
	case(3) !Node 3
!(3)----2 
! |   ^ |
! |  /  |
! | /   |
! 4-----1 
	    im = el(i)%elem(4)
		ip = el(i)%elem(2)
		ia = el(i)%elem(1)
	case(4) !Node 4
! 3-----2 
! | ^   |
! |  \  |
! |   \ |
!(4)----1
		im = el(i)%elem(1)
		ip = el(i)%elem(3)
		ia = el(i)%elem(2)
end select	
end subroutine posit
end subroutine velocity 



subroutine energy(dt,el,node,cell,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l
real(8) atr,abl,z_til(4),r_til(4),u_til(4),v_til(4),vol_til

cell(:)%um = 0d0 !Set zero cell velocity U-component
cell(:)%vm = 0d0 !Set zero cell velocity V-component

do i = 1,size(cell(:))
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		u_til(j) = 0.5d0*(node(l)%u+node(l)%u_l) ! U_til = 1/2(U+U_old) - intermediate velocity (See SALE2D or Samarsky)
		v_til(j) = 0.5d0*(node(l)%v+node(l)%v_l) ! V_til = 1/2(V+V_old) - intermediate velocity (See SALE2D or Samarsky)
		z_til(j) = node(l)%z + u_til(j)*dt !Calculation new intermediate z-position for define intermediate Volume
		r_til(j) = node(l)%r + v_til(j)*dt !Calculation new intermediate r-position for define intermediate Volume
		cell(i)%um = cell(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		cell(i)%vm = cell(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
	enddo

	atr = 0.5d0*((z_til(3)-z_til(2))*(r_til(1)-r_til(2))-(z_til(1)-z_til(2))*(r_til(3)-r_til(2))) !(right triangle) (also as Volume in subroutine Volume) 
	abl = 0.5d0*((z_til(1)-z_til(4))*(r_til(3)-r_til(4))-(z_til(3)-z_til(4))*(r_til(1)-r_til(4))) !(left triangle) (also as Volume in subroutine Volume)

	if (el(i)%rad == 2) then
        atr = atr*((r_til(1)+r_til(2)+r_til(3))/3d0)
        abl = abl*((r_til(1)+r_til(3)+r_til(4))/3d0)
	endif
	
	vol_til = abs(abl)+abs(atr) !intermedia Volume
	cell(i)%e = cell(i)%e-(cell(i)%p+cell(i)%q)*(vol_til-cell(i)%vol)/cell(i)%mas !Calculation inner energy (See SALE2D or Samarsky)
	cell(i)%etot = cell(i)%e + 0.5d0*(cell(i)%um*cell(i)%um+cell(i)%vm*cell(i)%vm)!Calculation Total energy 
enddo
end subroutine energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Need research interecation with boundary!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Conservative energy calculation!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energy_fromm(dt,el,node,cell,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k(4),ii,l
real(8) atr,abl,u_til,v_til,Pij,side

cell(:)%um = 0d0 !Set zero cell velocity U-component
cell(:)%vm = 0d0 !Set zero cell velocity V-component

do i = 1,size(cell(:))
    u_til = 0d0
    v_til = 0d0
    k = 0
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
        k(j) = l
		cell(i)%um = cell(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		cell(i)%vm = cell(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
		u_til = u_til+0.25d0*node(l)%u_l
		v_til = v_til+0.25d0*node(l)%v_l
	enddo
    cell(i)%etot = cell(i)%e+0.5d0*(u_til*u_til+v_til*v_til)
    
    do j = 1,el(i)%cont(9) !Cycle for all contact elements in el(:)%cont
		ii = el(i)%cont(j) !Number contact element 
		pij = ((cell(i)%p+cell(i)%q)*cell(ii)%mas+cell(i)%mas*(cell(ii)%p+cell(ii)%q))/(cell(ii)%mas+cell(i)%mas) !Weight mass-preassure. See YAQUII
		!    ------2------
		!           \
		!      P1,M1 \ P2,M2    Pij = (M1*P2+M2*P1)/(M1+M2)
		!             \
		!        ------1------
		select case(el(i)%cont(j+4)) !Define number side contact
		!      2
		!      |
		!   3----2   
		! 3-|    |-1 
		!   4----1
		!      |
		!      4
			case(1) !Right side
                side = 0.5d0*((node(k(2))%u+node(k(1))%u)*(node(k(2))%r-node(k(1))%r)-&
                (node(k(2))%v+node(k(1))%v)*(node(k(2))%z-node(k(1))%z)) 
                if (el(i)%rad == 2) side = side*0.5d0*(node(k(2))%r+node(k(1))%r)
                    
            case(2) !Top side
                side = 0.5d0*((node(k(3))%u+node(k(2))%u)*(node(k(3))%r-node(k(2))%r)-&
                (node(k(3))%v+node(k(2))%v)*(node(k(3))%z-node(k(2))%z))
                if (el(i)%rad == 2) side = side*0.5d0*(node(k(3))%r+node(k(2))%r)
                
            case(3) !Left sid
                side = 0.5d0*((node(k(4))%u+node(k(3))%u)*(node(k(4))%r-node(k(3))%r)-&
                (node(k(4))%v+node(k(3))%v)*(node(k(4))%z-node(k(3))%z))
                if (el(i)%rad == 2) side = side*0.5d0*(node(k(4))%r+node(k(3))%r)
                    
            case(4) !Bottom side
                side = 0.5d0*((node(k(1))%u+node(k(4))%u)*(node(k(1))%r-node(k(4))%r)-&
                (node(k(1))%v+node(k(4))%v)*(node(k(1))%z-node(k(4))%z))                
                if (el(i)%rad == 2) side = side*0.5d0*(node(k(1))%r+node(k(4))%r)
         end select
         
         cell(i)%etot = cell(i)%etot-dt/cell(i)%mas*pij*(cell(i)%sig*side) !Calculation Total energy
     enddo    
enddo
end subroutine energy_fromm

subroutine grid(dt,node,numer) 
type(nodes),intent(inout) :: node(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i

do i = 1,size(node(:)) !Cycle for all nodes
node(i)%z =  node(i)%z + node(i)%u*dt*numer%grid !Calculation lagrange z-postion node 
node(i)%r =  node(i)%r + node(i)%v*dt*numer%grid !Calculation lagrange r-postion node 
!numer%grid - coff shift vertex (if lagrange  = 1, Euler = 0, ALE = between 0-1)
enddo
end subroutine grid


subroutine advect(dt,el,node,cell,bou,numer) !Subroutine calculation advection
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
type (boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l,ii
real(8) z_p(4),r_p(4),vol,FV,a_V

if (numer%grid<1) then !if numer%grid = 1 then bypassed this block 
	node(:)%mas_til_v = 0d0
	
	node(:)%mas_til_vu = node(:)%mas_v*node(:)%u
	node(:)%mas_til_vv = node(:)%mas_v*node(:)%v
	
	do i = 1,size(cell(:)%vol)
		cell(i)%mas_til = cell(i)%mas !Set intermedia mass equivalent mass
		cell(i)%Me_til = cell(i)%mas*cell(i)%etot !Set intermedia Mass*energy (Total energy)

		! el(i)%cont(j)  (contact_1 contact_2 contact_3 contact_4 side_1 side_2 side_3 side_4) 

		ext:do j = 1,el(i)%cont(9) !Cycle for all contact elements in el(:)%cont
			ii = el(i)%cont(j) !Number contact element  
			select case(el(i)%cont(j+4)) !Define number side contact
				!      2
				!      |
				!   3----2   
				! 3-|    |-1 
				!   4----1
				!      |
				!      4
				case(1) !Right side
				!         z_p(3)    shift node
				!          /         /
				!         v         v
				!   3----2.....z_p(2)      z_p(3)---z_p(2)
				!   |    |     .       ->  |         |      <- Cell shift
				!   4----1.....z_p(1)      z_p(4)---z_p(1)  
				!        ^       ^
				!         \       \  
				!        z_p(4)   shift node
					z_p(3) = node(el(i)%elem(2))%z
					z_p(4) = node(el(i)%elem(1))%z 
					z_p(1) = z_p(4) + node(el(i)%elem(1))%u*dt*(1d0 - numer%grid)   !z(final position) = z(lagrange) - u(velocity shift node in back position) *dt
					z_p(2) = z_p(3) + node(el(i)%elem(2))%u*dt*(1d0 - numer%grid)
				!r -component
					r_p(3) = node(el(i)%elem(2))%r
					r_p(4) = node(el(i)%elem(1))%r 
					r_p(1) = r_p(4) + node(el(i)%elem(1))%v*dt*(1d0 - numer%grid)
					r_p(2) = r_p(3) + node(el(i)%elem(2))%v*dt*(1d0 - numer%grid)
												
				case(2)! Top side
				!z_p(3)..z_p(2)
				!   .    .  
				!   .    . 
				!   .    .  
				!   3----2
				!   |    |   
				!   4----1
					z_p(1) = node(el(i)%elem(2))%z 
					z_p(4) = node(el(i)%elem(3))%z
					z_p(2) = z_p(1) + node(el(i)%elem(2))%u*dt*(1d0 - numer%grid)
					z_p(3) = z_p(4) + node(el(i)%elem(3))%u*dt*(1d0 - numer%grid)
				
					r_p(1) = node(el(i)%elem(2))%r 
					r_p(4) = node(el(i)%elem(3))%r
					r_p(2) = r_p(1) + node(el(i)%elem(2))%v*dt*(1d0 - numer%grid)
					r_p(3) = r_p(4) + node(el(i)%elem(3))%v*dt*(1d0 - numer%grid)
			
				case(3)! Left side
				! 
				!z_p(3).....3----2
				!           |    |   
				!z_p(4).....4----1
					z_p(1) = node(el(i)%elem(4))%z 
					z_p(2) = node(el(i)%elem(3))%z 
					z_p(3) = z_p(2) + node(el(i)%elem(3))%u*dt*(1d0 - numer%grid)
					z_p(4) = z_p(1) + node(el(i)%elem(4))%u*dt*(1d0 - numer%grid)
				
					r_p(1) = node(el(i)%elem(4))%r 
					r_p(2) = node(el(i)%elem(3))%r 
					r_p(3) = r_p(2) + node(el(i)%elem(3))%v*dt*(1d0 - numer%grid)
					r_p(4) = r_p(1) + node(el(i)%elem(4))%v*dt*(1d0 - numer%grid)
				case(4) !Bottom side
				! 
				!    3----2
				!    |    |   
				!    4----1
				!    .    .  
				!    .    . 
				!    .    . 
				!z_p(4)...z_p(1)
					z_p(2) = node(el(i)%elem(1))%z 
					z_p(3) = node(el(i)%elem(4))%z 
					z_p(1) = z_p(2) + node(el(i)%elem(1))%u*dt*(1d0 - numer%grid)
					z_p(4) = z_p(3) + node(el(i)%elem(4))%u*dt*(1d0 - numer%grid)
				
					r_p(2) = node(el(i)%elem(1))%r 
					r_p(3) = node(el(i)%elem(4))%r 
					r_p(1) = r_p(2) + node(el(i)%elem(1))%v*dt*(1d0 - numer%grid)
					r_p(4) = r_p(3) + node(el(i)%elem(4))%v*dt*(1d0 - numer%grid)
				case default !When we have zero index in el(:)%cont
					exit ext
			end select
		
			FV = -cell(i)%sig*volume_adv(z_p,r_p,el(i)%rad) !Calculation volume cell-shift
			a_V = numer%a0*sign(1d0,Fv) !For stability a0 = 1d0 - change variable (0 - unstable 1- more diffuse, using 0 - 1)  
			cell(i)%mas_til = cell(i)%mas_til+ 0.5d0*Fv*((1d0+a_V)*cell(ii)%rho+(1d0-a_V)*cell(i)%rho) !Calculation mass  Fv(0.5[rho_ii+rho_i]) - See Fletcher (vol. 1) or SALE2D (YAQUII)
			cell(i)%Me_til = cell(i)%Me_til + 0.5d0*Fv*((1d0+a_V)*cell(ii)%rho*cell(ii)%etot+(1d0-a_V)*cell(i)%rho*cell(i)%etot) !Mass*Energy 
            do k = 1,4 !Cycle for all nodes in cell		
                l = el(i)%elem(k) !Number node 
			!We should reestablish our velocity from flux momentum
			!   
			!  5-------3-------2
			!  |   dum2|  dum1 |    
			!  |     \ |  /    |
			!  |      v| v     |
			!  6------(4)------1   (M*U)^(n+1) = (M*U)^n + 0.25*dUM 
			!  |      ^| ^     |   (M*U)^(n+1) = (M*U)^n + 0.25*dVM  
			!  |     / |  \    |
			!  |  dum3 | dum4  |
			!  7-------8-------9
			!
node(l)%mas_til_vu = node(l)%mas_til_vu+0.125*Fv*((1d0+a_V)*cell(ii)%rho*cell(ii)%Um+(1d0-a_V)*cell(i)%rho*cell(i)%Um) !Mass_v*U-Velocity
node(l)%mas_til_vv = node(l)%mas_til_vv+0.125*Fv*((1d0+a_V)*cell(ii)%rho*cell(ii)%Vm+(1d0-a_V)*cell(i)%rho*cell(i)%Vm) !Mass_v*V-Velocity
            enddo
		enddo ext
		do j = 1,4 !Cycle for all nodes in cell
			l = el(i)%elem(j) !Number node
			node(l)%mas_til_v = node(l)%mas_til_v+0.25d0*cell(i)%mas_til ! M_v = 1/4(M1+M2+M3+M4) 
		enddo
	enddo	

	
    cell(:)%mas = cell(:)%mas_til
    node(:)%u = node(:)%mas_til_vu/node(:)%mas_til_v
    node(:)%v = node(:)%mas_til_vv/node(:)%mas_til_v
    cell(:)%etot = cell(:)%Me_til/cell(:)%mas_til
    node(:)%mas_v = node(:)%mas_til_v
endif

cell(:)%um = 0d0
cell(:)%vm = 0d0
!Calculation inner energy from new total energy
do i = 1,size(cell(:)%vol)
	do j = 1,4 !Cycle for all nodes in cell
		l = el(i)%elem(j) !Number node
        cell(i)%um = cell(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		cell(i)%vm = cell(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
	enddo	
    cell(i)%e = cell(i)%etot - 0.5d0*(cell(i)%um*cell(i)%um+cell(i)%vm*cell(i)%vm)
enddo
call boundary_flow(bou,node,cell,el)
node(:)%u_l = node(:)%u
node(:)%v_l = node(:)%v

contains
real(8) function volume_adv(z,r,rad)
real(8) z(4),r(4),rad
real(8) atr,abl

atr = 0.5d0*((z(3)-z(2))*(r(1)-r(2))-(z(1)-z(2))*(r(3)-r(2))) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
abl = 0.5d0*((z(1)-z(4))*(r(3)-r(4))-(z(3)-z(4))*(r(1)-r(4))) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

if (rad == 2) then
    atr = atr*((r(1)+r(2)+r(3))/3d0)
    abl = abl*((r(1)+r(3)+r(4))/3d0)
endif

volume_adv = abl+atr
end function volume_adv
end subroutine advect


subroutine difference(i,l,vec,dvecdz,dvecdr,el,node)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
integer(4),intent(in) :: l(4),i
real(8),intent(in) :: vec(4)
real(8),intent(out) :: dvecdz,dvecdr
real(8) area
area = (node(l(2))%z-node(l(4))%z)*(node(l(3))%r-node(l(1))%r)-(node(l(1))%z-node(l(3))%z)*(node(l(4))%r-node(l(2))%r) !Calculation Area
dvecdz = (vec(2)-vec(4))*(node(l(3))%r-node(l(1))%r)-(vec(1)-vec(3))*(node(l(4))%r-node(l(2))%r)/area
dvecdr = (vec(3)-vec(1))*(node(l(2))%r-node(l(4))%r)-(vec(2)-vec(4))*(node(l(3))%r-node(l(1))%r)/area
end subroutine difference


subroutine divergence(dvec_zdz,dvec_rdr,r,vec,rad,diver) !Divergence vector var
real(8),intent(in) :: dvec_zdz,dvec_rdr,r(4),vec(4),rad
real(8),intent(out) :: diver
diver = dvec_zdz+dvec_rdr+(rad-1d0)*(sum(vec))/(sum(r))
end subroutine divergence


subroutine rotor(dvec_zdr,dvec_rdz,r,vec,rad,rot) !Rotor vector var
real(8),intent(in) :: dvec_zdr,dvec_rdz,r(4),vec(4),rad
real(8),intent(out) :: rot
rot = dvec_zdr-dvec_rdz!+(rad-1d0)*(sum(vec))/(sum(r))
end subroutine rotor

end module hydro
