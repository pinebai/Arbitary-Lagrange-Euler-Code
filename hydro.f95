module hydro
use object
use eos
use type_boundary
contains

subroutine phase0(el,node,phy,bou,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
integer(4) i,j,l(4)
real(8) atr,abl

node(:)%mas_v = 0d0
node(:)%num_cont = 0d0
node(:)%mark = 0
!do i = 1,size(bou%var(:,1))
!    l(1:3) = bou%var(i,1:3)
!	if (bou%type_bound(l(1)) == 5) el(l(3))%rad = 1d0
!enddo

do i = 1,size(phy(:)%vol)

! 3-----2
! |     |	
! 4-----1
!
!Node equivalent 
!Sale
!1 2 3 4
!el(i)%elem
!2 3 4 5

	l(1) = el(i)%elem(2) 
	l(2) = el(i)%elem(3)
	l(3) = el(i)%elem(4)
	l(4) = el(i)%elem(5)

!Calculation triangles left and right in cell 
	atr = 0.5d0*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
	abl = 0.5d0*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))
 
	atr = atr*((node(l(1))%r+node(l(2))%r+node(l(3))%r)/3d0)**(el(i)%rad-1d0)
	abl = abl*((node(l(1))%r+node(l(3))%r+node(l(4))%r)/3d0)**(el(i)%rad-1d0)
	
	phy(i)%sig = sign(1d0,abl+atr) !Sign volume (Very important)
	phy(i)%vol = abs(abl)+abs(atr) !Calculaion volume cell
	
	phy(i)%vol_old = phy(i)%vol !Set first old volume 
	phy(i)%mas = phy(i)%rho*phy(i)%vol ! Calculation massa in cell
	
	call gas(i,phy)

!Set old velocity
	node(l(1))%u_l = node(l(1))%u 
	node(l(2))%u_l = node(l(2))%u
	node(l(3))%u_l = node(l(3))%u
	node(l(4))%u_l = node(l(4))%u

	node(l(1))%v_l = node(l(1))%v
	node(l(2))%v_l = node(l(2))%v
	node(l(3))%v_l = node(l(3))%v
	node(l(4))%v_l = node(l(4))%v

!Calculation vertex mass and number cell near node
	do j = 1,4 !cycle for all nodes in current cell 
		node(l(j))%mas_v = node(l(j))%mas_v + 0.25d0*phy(i)%mas
		node(l(j))%num_cont = node(l(j))%num_cont + 1d0
	enddo
enddo
call boundary_flow(bou,node,phy,el)
end subroutine phase0


subroutine phase1(dt,el,node,phy,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,l(4)
real(8) atr, abl
real(8) Dvol

do i = 1,size(phy(:)%vol)
	!Calculation EOS
	call gas(i,phy)
	
	!Calculation volume,mass and density
	!first material then shift to 1
	l(1) = el(i)%elem(2) 
	l(2) = el(i)%elem(3)
	l(3) = el(i)%elem(4)
	l(4) = el(i)%elem(5)

	atr = 0.5d0*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
	abl = 0.5d0*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

	atr = atr*((node(l(1))%r+node(l(2))%r+node(l(3))%r)/3d0)**(el(i)%rad-1d0)
	abl = abl*((node(l(1))%r+node(l(3))%r+node(l(4))%r)/3d0)**(el(i)%rad-1d0)
	
	phy(i)%vol_old = phy(i)%vol !Set old volume 
	phy(i)%vol = abs(abl)+abs(atr) !Calcualtion new volume 
	phy(i)%rho = phy(i)%mas/phy(i)%vol !Calculation density
	phy(i)%mas = phy(i)%rho*phy(i)%vol !Calculation massa
	
	!Calculation artvisc
	Dvol = 2d0*(phy(i)%vol - phy(i)%vol_old)/(phy(i)%vol + phy(i)%vol_old)/dt !Calculation Divergence (Equivalen Hemp3D)
	phy(i)%q = min(0d0,Dvol)*(numer%art*phy(i)%rho*Dvol*phy(i)%Vol**0.666) !Calculation Divergence (Equivalen Hemp3D)	
	
enddo
end subroutine phase1


subroutine velocity(dt,el,node,phy,bou,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type (boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l,im,ip,ia
real(8) acc_z,acc_r,rad_com(2)
real(8) an,ksi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Relaxation grid
!This cycle for calculation realxation lagrange grid (See SALE)
!an and ksi - var parameters
!See more comments in  call posit(...)
an = 0.005d0
ksi = 0.3d0
do i = 1,size(phy(:)) !
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j+1) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		call posit(i,j,phy,el,im,ip,ia) !Function calculation position + (ip) and - (im)			
node(l)%u = node(l)%u + 0.25d0*an*(0.5d0*(1d0+ksi)*(node(ip)%u_l+node(im)%u_l)-ksi*node(ia)%u_l-node(l)%u_l)
node(l)%v = node(l)%v + 0.25d0*an*(0.5d0*(1d0+ksi)*(node(ip)%v_l+node(im)%v_l)-ksi*node(ia)%v_l-node(l)%v_l)
	enddo	
enddo

call boundary_flow(bou,node,phy,el)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1,size(phy(:)) !
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j+1) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		call posit(i,j,phy,el,im,ip,ia) !Function calculation position 
		
		rad_com(1) = (0.5d0*(node(ip)%r+node(im)%r))**(el(i)%rad-1d0) !Radial component for z
		rad_com(2) = node(l)%r**(el(i)%rad-1d0)						  !Radial component for r
		
		acc_z = 0.5d0*phy(i)%sig*rad_com(1)*(phy(i)%p+phy(i)%q)*(node(ip)%r-node(im)%r)/node(l)%mas_v !Caluclation z-component acceleration 
		acc_r = 0.5d0*phy(i)%sig*rad_com(2)*(phy(i)%p+phy(i)%q)*(node(im)%z-node(ip)%z)/node(l)%mas_v !Caluclation r-component acceleration
			
		node(l)%u = node(l)%u - acc_z*dt !Caluclation z-component velocity 
		node(l)%v = node(l)%v - acc_r*dt !Caluclation r-component velocity 
	enddo	
enddo	

call boundary_flow(bou,node,phy,el)	
contains

subroutine posit(i,j,phy,el,im,ip,ia)
type(elements),intent(inout) :: el(:)
type(physics),intent(inout) :: phy(:)
integer(4),intent(in) :: i,j
integer(4),intent(out) :: im,ip,ia
integer(4) l 
l = 1
select case(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This countur for node (o)
!im,ip - adjacent nodes
!ia - oposite node
!sig - direction rotation
!
! n-----ip----ia 
! |   / | ^   |
! |	 /  |  \  |
! |	v   |	\ |
! n----(o)----im 
! | \   |   ^ |
! |	 \  |  /  |
! |	  v | /   |
! n-----n-----n
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	case(1) !Node 1
! 3-----2 
! |   / |
! |	 /  |
! |	v   |
! 4----(1)			 	  	
		im = el(i)%elem(l+2)    
		ip = el(i)%elem(l+4) 
		ia = el(i)%elem(l+3)
	case(2) !Node 2
! 3----(2) 
! | \   |
! |	 \  |
! |	  v |
! 4-----1		 
		im = el(i)%elem(l+3)
		ip = el(i)%elem(l+1)
		ia = el(i)%elem(l+4)
	case(3) !Node 3
!(3)----2 
! |   ^ |
! |	 /  |
! |	/   |
! 4-----1 		
		im = el(i)%elem(l+4)
		ip = el(i)%elem(l+2)
		ia = el(i)%elem(l+1)
	case(4) !Node 4
! 3-----2 
! | ^   |
! |	 \  |
! |	  \ |
!(4)----1				
		im = el(i)%elem(l+1)
		ip = el(i)%elem(l+3)
		ia = el(i)%elem(l+2)
end select	
end subroutine posit
end subroutine velocity 



subroutine energy(dt,el,node,phy,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l
real(8) atr,abl,z_til(4),r_til(4),u_til(4),v_til(4),vol_til

phy(:)%um = 0d0 !Set zero cell velocity U-component
phy(:)%vm = 0d0 !Set zero cell velocity V-component

do i = 1,size(phy(:))
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j+1) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		u_til(j) = 0.5d0*(node(l)%u+node(l)%u_l) ! U_til = 1/2(U+U_old) - intermediate velocity (See SALE2D or Samarsky)
		v_til(j) = 0.5d0*(node(l)%v+node(l)%v_l) ! V_til = 1/2(V+V_old) - intermediate velocity (See SALE2D or Samarsky)
		z_til(j) = node(l)%z + u_til(j)*dt !Calculation new intermediate z-position for define intermediate Volume
		r_til(j) = node(l)%r + v_til(j)*dt !Calculation new intermediate r-position for define intermediate Volume
		phy(i)%um = phy(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		phy(i)%vm = phy(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
	enddo

	atr = 0.5d0*((z_til(3)-z_til(2))*(r_til(1)-r_til(2))-(z_til(1)-z_til(2))*(r_til(3)-r_til(2))) !(right triangle) (also as Volume in subroutine Volume) 
	abl = 0.5d0*((z_til(1)-z_til(4))*(r_til(3)-r_til(4))-(z_til(3)-z_til(4))*(r_til(1)-r_til(4))) !(left triangle) (also as Volume in subroutine Volume)

	atr = atr*((r_til(1)+r_til(2)+r_til(3))/3d0)**(el(i)%rad-1d0)
	abl = abl*((r_til(1)+r_til(3)+r_til(4))/3d0)**(el(i)%rad-1d0)
	
	vol_til = abs(abl)+abs(atr) !intermedia Volume
	phy(i)%e = phy(i)%e-(phy(i)%p+phy(i)%q)*(vol_til-phy(i)%vol)/phy(i)%mas !Calculation inner energy (See SALE2D or Samarsky)
	phy(i)%etot = phy(i)%e + 0.5d0*(phy(i)%um*phy(i)%um+phy(i)%vm*phy(i)%vm)!Calculation Total energy 
enddo
end subroutine energy

subroutine energy_fromm(dt,el,node,phy,numer)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k(4),ii,kl,l
real(8) atr,abl,u_til,v_til,Pij,side

kl = 1

phy(:)%um = 0d0 !Set zero cell velocity U-component
phy(:)%vm = 0d0 !Set zero cell velocity V-component

do i = 1,size(phy(:))
    u_til = 0d0
    v_til = 0d0
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(kl+j) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
        k(j) = l
		phy(i)%um = phy(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		phy(i)%vm = phy(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
		u_til = u_til+0.25d0*node(l)%u_l
		v_til = v_til+0.25d0*node(l)%v_l
	enddo
    phy(i)%etot = phy(i)%e+0.5d0*(u_til*u_til+v_til*v_til)
    
    do j = 1,el(i)%cont(9) !Cycle for all contact elements in el(:)%cont
		ii = el(i)%cont(j) !Number contact element 
		pij = ((phy(i)%p+phy(i)%q)*phy(ii)%mas+phy(i)%mas*(phy(ii)%p+phy(ii)%q))/(phy(ii)%mas+phy(i)%mas) !Weight mass-preassure 
		!    ------2------
		!           \
		!      P1,M1 \ P2,M2    Pij = (M1*P2+M2*P1)/(M1+M2)
		!             \
		!        ------1------
		select case(el(i)%cont(j+4)) !Define number side contact
			!      2
			!      |
			! 	3----2   
			! 3-|    |-1	 
			! 	4----1
			!	   |
			!      4
			case(1) !Right side
                side = 0.5*((node(k(2))%u+node(k(1))%u)*(node(k(2))%r-node(k(1))%r)-&
(node(k(2))%v+node(k(1))%v)*(node(k(2))%z-node(k(1))%z))*(0.5d0*(node(k(2))%r+node(k(1))%r))**((el(i)%rad-1d0)) 
			case(2) !Top side
                side = 0.5*((node(k(3))%u+node(k(2))%u)*(node(k(3))%r-node(k(2))%r)-&
(node(k(3))%v+node(k(2))%v)*(node(k(3))%z-node(k(2))%z))*(0.5d0*(node(k(3))%r+node(k(2))%r))**((el(i)%rad-1d0)) 
			case(3) !Left sid
                side = 0.5*((node(k(4))%u+node(k(3))%u)*(node(k(4))%r-node(k(3))%r)-&
(node(k(4))%v+node(k(3))%v)*(node(k(4))%z-node(k(3))%z))*(0.5d0*(node(k(4))%r+node(k(3))%r))**((el(i)%rad-1d0)) 
			case(4) !Bottom side
                side = 0.5*((node(k(1))%u+node(k(4))%u)*(node(k(1))%r-node(k(4))%r)-&
(node(k(1))%v+node(k(4))%v)*(node(k(1))%z-node(k(4))%z))*(0.5d0*(node(k(1))%r+node(k(4))%r))**((el(i)%rad-1d0))                   
         end select
         phy(i)%etot = phy(i)%etot-dt/phy(i)%mas*pij*(phy(i)%sig*side) !Calculation Total energy
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


subroutine advect(dt,el,node,phy,bou,numer) !Subroutine calculation advection
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type (boundary),intent(inout) :: bou
type(numerical),intent(inout) :: numer
real(8),intent(in) :: dt
integer(4) i,j,k,l,ii,kl
real(8) z_p(4),r_p(4),vol,FV,a_V
kl = 1 !For shift index in el(:)%elem
if (numer%grid<1) then !if numer%grid = 1 then bypassed this block 

	do i = 1,size(phy(:)%vol)
		phy(i)%dum = 0d0 !Set zero flux U -momentum 
		phy(i)%dvm = 0d0 !Set zero flux V -momentum 
		phy(i)%mas_til = phy(i)%mas !Set intermedia mass equivalent mass
		phy(i)%Me_til = phy(i)%mas*phy(i)%etot !Set intermedia Mass*energy (Total energy)

		! el(i)%cont(j)  (contact_1 contact_2 contact_3 contact_4 side_1 side_2 side_3 side_4) 

		ext:do j = 1,el(i)%cont(9) !Cycle for all contact elements in el(:)%cont
			ii = el(i)%cont(j) !Number contact element  
			select case(el(i)%cont(j+4)) !Define number side contact			
				!      2
				!      |
				! 	3----2   
				! 3-|    |-1	 
				! 	4----1
				!	   |
				!      4
				case(1) !Right side
				!		    z_p(3)    shift node
				!          /         /
				!         v         v
				! 	3----2.....z_p(2)      z_p(3)---z_p(2)
				!   |    | 	   .       ->  |         |      <- Cell shift
				! 	4----1.....z_p(1)      z_p(4)---z_p(1)  
				!        ^	      ^
				!         \        \  
				!        z_p(4)    shift node
					z_p(3) = node(el(i)%elem(kl+2))%z
					z_p(4) = node(el(i)%elem(kl+1))%z 
					z_p(1) = z_p(4) + node(el(i)%elem(kl+1))%u*dt*(1d0 - numer%grid)   !z(final position) = z(lagrange) - u(velocity shift node in back position) *dt
					z_p(2) = z_p(3) + node(el(i)%elem(kl+2))%u*dt*(1d0 - numer%grid)
				!r -component
					r_p(3) = node(el(i)%elem(kl+2))%r
					r_p(4) = node(el(i)%elem(kl+1))%r 
					r_p(1) = r_p(4) + node(el(i)%elem(kl+1))%v*dt*(1d0 - numer%grid)
					r_p(2) = r_p(3) + node(el(i)%elem(kl+2))%v*dt*(1d0 - numer%grid)
												
				case(2)! Top side
				!z_p(3)..z_p(2)	
				!	.	 .  
				!   .    . 
				!   .    .  
				! 	3----2
				!   |    | 	   
				! 	4----1
					z_p(1) = node(el(i)%elem(kl+2))%z 
					z_p(4) = node(el(i)%elem(kl+3))%z
					z_p(2) = z_p(1) + node(el(i)%elem(kl+2))%u*dt*(1d0 - numer%grid)
					z_p(3) = z_p(4) + node(el(i)%elem(kl+3))%u*dt*(1d0 - numer%grid)
				
					r_p(1) = node(el(i)%elem(kl+2))%r 
					r_p(4) = node(el(i)%elem(kl+3))%r
					r_p(2) = r_p(1) + node(el(i)%elem(kl+2))%v*dt*(1d0 - numer%grid)
					r_p(3) = r_p(4) + node(el(i)%elem(kl+3))%v*dt*(1d0 - numer%grid)
			
				case(3)! Left side
				! 
				!z_p(3).....3----2
				!   		|    | 	   
				!z_p(4).....4----1
					z_p(1) = node(el(i)%elem(kl+4))%z 
					z_p(2) = node(el(i)%elem(kl+3))%z 
					z_p(3) = z_p(2) + node(el(i)%elem(kl+3))%u*dt*(1d0 - numer%grid)
					z_p(4) = z_p(1) + node(el(i)%elem(kl+4))%u*dt*(1d0 - numer%grid)
				
					r_p(1) = node(el(i)%elem(kl+4))%r 
					r_p(2) = node(el(i)%elem(kl+3))%r 
					r_p(3) = r_p(2) + node(el(i)%elem(kl+3))%v*dt*(1d0 - numer%grid)
					r_p(4) = r_p(1) + node(el(i)%elem(kl+4))%v*dt*(1d0 - numer%grid)
				case(4) !Bottom side
				! 
				!    3----2
				!    |    | 	   
				!    4----1
				!	 .    .  
				!    .    . 
				!	 .    . 
				!z_p(4)...z_p(1)		
					z_p(2) = node(el(i)%elem(kl+1))%z 
					z_p(3) = node(el(i)%elem(kl+4))%z 
					z_p(1) = z_p(2) + node(el(i)%elem(kl+1))%u*dt*(1d0 - numer%grid)
					z_p(4) = z_p(3) + node(el(i)%elem(kl+4))%u*dt*(1d0 - numer%grid)
				
					r_p(2) = node(el(i)%elem(kl+1))%r 
					r_p(3) = node(el(i)%elem(kl+4))%r 
					r_p(1) = r_p(2) + node(el(i)%elem(kl+1))%v*dt*(1d0 - numer%grid)
					r_p(4) = r_p(3) + node(el(i)%elem(kl+4))%v*dt*(1d0 - numer%grid)
				case default !When we have zero index in el(:)%cont
					exit ext	
			end select
		
			FV = -phy(i)%sig*volume_adv(z_p,r_p,el(i)%rad) !Calculation volume cell-shift
			a_V = numer%a0*sign(1d0,Fv) !For stability a0 = 1d0 - change variable (0 - unstable 1- more diffuse, using 0 - 1)  
			phy(i)%mas_til = phy(i)%mas_til+ 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho+(1d0-a_V)*phy(i)%rho) !Calculation mass  Fv(0.5[rho_ii+rho_i]) - See Fletcher (vol. 1) or SALE2D (YAQUII)
			phy(i)%dum = phy(i)%dum + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%Um+(1d0-a_V)*phy(i)%rho*phy(i)%Um) !flux U -momentum 
			phy(i)%dvm = phy(i)%dvm + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%Vm+(1d0-a_V)*phy(i)%rho*phy(i)%Vm) !flux V -momentum 
			phy(i)%Me_til = phy(i)%Me_til + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%etot+(1d0-a_V)*phy(i)%rho*phy(i)%etot) !Mass*Energy 
		
		enddo ext
	enddo	
	node(:)%mas_til_v = 0d0
	do i = 1,size(phy(:))
		do j = 1,4 !Cycle for all nodes in cell
			l = el(i)%elem(j+1) !Number node
			node(l)%mas_til_v = node(l)%mas_til_v+0.25d0*phy(i)%mas_til ! M_v = 1/4(M1+M2+M3+M4) 
		enddo
	enddo	
	!Using u_l and v_l as intermedia velocity for optimization
	node(:)%u_l = node(:)%u !Set U-velocity
	node(:)%v_l = node(:)%v !Set V-velocity
	node(:)%u = 0d0 !Set U-velocity on zero
	node(:)%v = 0d0 !Set V-velocity on zero
	!Find new vertex velocity
	do i = 1,size(phy(:))
		do j = 1,4 !Cycle for all nodes in cell
			l = el(i)%elem(j+1) !Number node 
			!We should reestablish our velocity from flux momentum
			! But when we have incomplete set cells near node (less <4), we should division on count cell (sum cell near node)
			!   
			!  5-------3-------2			
			!  |   dum2|  dum1 |    
			!  |      \|  /    |
			!  |  	  v| v	   |	
			!  6------(4)------1
			!  |      ^| ^     |
			!  |	 / | \     |
			!  |  dum3 | dum4  |
			!  7-------8-------9	
			! (MU)^(n+1) = (MU)^(n)/count + 1/4dUM  	
			! where "count" - is sum cell near node

			node(l)%u = node(l)%u + 0.25d0*phy(i)%dum/node(l)%mas_til_v + node(l)%mas_v*node(l)%u_l/node(l)%num_cont/node(l)%mas_til_v
			node(l)%v = node(l)%v + 0.25d0*phy(i)%dvm/node(l)%mas_til_v + node(l)%mas_v*node(l)%v_l/node(l)%num_cont/node(l)%mas_til_v
		enddo	
	enddo
	node(:)%mas_v = node(:)%mas_til_v
	do i = 1,size(phy(:)%vol)
		phy(i)%mas = phy(i)%mas_til
		phy(i)%etot = phy(i)%Me_til/phy(i)%mas_til
	enddo
	call boundary_flow(bou,node,phy,el)
endif

node(:)%u_l = node(:)%u
node(:)%v_l = node(:)%v


phy(:)%um = 0d0
phy(:)%vm = 0d0
!Calculation inner energy from new total energy
do i = 1,size(phy(:)%vol)
    us = 0
	do j = 1,4 !Cycle for all nodes in cell
		l = el(i)%elem(j+1) !Number node
        phy(i)%um = phy(i)%um + 0.25d0*node(l)%u !Calculation cell-centered velocity U
		phy(i)%vm = phy(i)%vm + 0.25d0*node(l)%v !Calculation cell-centered velocity V
	enddo	
    phy(i)%e = phy(i)%etot - 0.5d0*(phy(i)%um*phy(i)%um+phy(i)%vm*phy(i)%vm)
enddo

contains
real(8) function volume_adv(z,r,rad)
real(8) z(4),r(4),rad
real(8) atr,abl

atr = 0.5d0*((z(3)-z(2))*(r(1)-r(2))-(z(1)-z(2))*(r(3)-r(2))) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
abl = 0.5d0*((z(1)-z(4))*(r(3)-r(4))-(z(3)-z(4))*(r(1)-r(4))) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

atr = atr*((r(1)+r(2)+r(3))/3d0)**(rad-1d0)
abl = abl*((r(1)+r(3)+r(4))/3d0)**(rad-1d0)

volume_adv = abl+atr
end function volume_adv
end subroutine advect


subroutine difference(i,dudz,dudr,dvdz,dvdr,el,node)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
integer(4),intent(in) :: i
real(8),intent(out) :: dudz,dudr,dvdz,dvdr
real(8) area
integer(4) l(4)

l(1) = el(i)%elem(2) 
l(2) = el(i)%elem(3)
l(3) = el(i)%elem(4)
l(4) = el(i)%elem(5)

area = (node(l(2))%z-node(l(4))%z)*(node(l(3))%r-node(l(1))%r)-(node(l(1))%z-node(l(3))%z)*(node(l(4))%r-node(l(2))%r) !Calculation Area

dudz = (node(l(2))%u-node(l(4))%u)*(node(l(3))%r-node(l(1))%r)-(node(l(1))%u-node(l(3))%u)*(node(l(4))%r-node(l(2))%r)/area !Calculation dudz
dudr = (node(l(3))%u-node(l(1))%u)*(node(l(2))%r-node(l(4))%r)-(node(l(2))%u-node(l(4))%u)*(node(l(3))%r-node(l(1))%r)/area!+& !Calculation dudr
!(el(i)%rad-1d0)*(node(l(1))%u+node(l(2))%u+node(l(3))%u+node(l(4))%u)/(node(l(1))%r+node(l(2))%r+node(l(3))%r+node(l(4))%r)

dvdz = (node(l(2))%v-node(l(4))%v)*(node(l(3))%r-node(l(1))%r)-(node(l(1))%v-node(l(3))%v)*(node(l(4))%r-node(l(2))%r)/area !Calculation dvdz
dvdr = (node(l(3))%v-node(l(1))%v)*(node(l(2))%r-node(l(4))%r)-(node(l(2))%v-node(l(4))%v)*(node(l(3))%r-node(l(1))%r)/area!+& !Calculation dvdr
!(el(i)%rad-1d0)*(node(l(1))%u+node(l(2))%u+node(l(3))%u+node(l(4))%u)/(node(l(1))%r+node(l(2))%r+node(l(3))%r+node(l(4))%r)
end subroutine difference

end module hydro
