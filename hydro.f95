module hydro
use object
contains

subroutine volume0(el,node,phy)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
integer(4) i,j,l(4)
real(8) atr,abl

node(:)%mas_v = 0d0
node(:)%num_cont = 0d0

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
	atr = 0.5*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
	abl = 0.5*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

	phy(i)%sig = sign(1d0,abl+atr) !Sign volume (Very important)
	phy(i)%vol = abs(abl)+abs(atr) !Calculaion volume cell
	
	phy(i)%vol_old = phy(i)%vol !Set first old volume 
	phy(i)%mas = phy(i)%rho*phy(i)%vol ! Calculation massa in cell

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
end subroutine volume0


subroutine volume(el,node,phy)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
integer(4) i,l(4)
real(8) atr, abl

do i = 1,size(phy(:)%vol)
!first material then shift to 1
l(1) = el(i)%elem(2) 
l(2) = el(i)%elem(3)
l(3) = el(i)%elem(4)
l(4) = el(i)%elem(5)

atr = 0.5d0*((node(l(3))%z-node(l(2))%z)*(node(l(1))%r-node(l(2))%r)-(node(l(1))%z-node(l(2))%z)*(node(l(3))%r-node(l(2))%r)) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
abl = 0.5d0*((node(l(1))%z-node(l(4))%z)*(node(l(3))%r-node(l(4))%r)-(node(l(3))%z-node(l(4))%z)*(node(l(1))%r-node(l(4))%r)) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))

phy(i)%vol_old = phy(i)%vol !Set old volume 
phy(i)%vol = abs(abl)+abs(atr) !Calcualtion new volume 
phy(i)%rho = phy(i)%mas/phy(i)%vol !Calculation density
phy(i)%mas = phy(i)%rho*phy(i)%vol !Calculation massa
enddo
end subroutine volume


subroutine velocity(dt,el,node,phy)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
real(8),intent(in) :: dt
integer(4) i,j,k,l,im,ip
real(8) sig

do i = 1,size(phy(:)) !
	do j = 1,4 !Cycle for all nodes cell
		l = el(i)%elem(j+1) !index in el(:)%elem next (N_mat node_1 node_2 node_3 node_4)
		call posit(i,j,phy,el,im,ip,sig) !Function calculation position + (ip) and - (im)

		node(l)%u = node(l)%u - 0.5d0*sig*dt*(phy(i)%p+phy(i)%q)*(node(ip)%r-node(im)%r)/node(l)%mas_v
		node(l)%v = node(l)%v - 0.5d0*sig*dt*(phy(i)%p+phy(i)%q)*(node(im)%z-node(ip)%z)/node(l)%mas_v
	enddo	
enddo	

end subroutine velocity 

subroutine posit(i,j,phy,el,im,ip,sig)
type(elements),intent(inout) :: el(:)
type(physics),intent(inout) :: phy(:)
integer(4),intent(in) :: i,j
integer(4),intent(out) :: im,ip
real(8),intent(out) :: sig

select case(j)
	case(1) !Node 1
! 3+1------2+1    ^--->
! |		   |      |   |
! 4+1-----(1+1)	  ^-- v	 	  	
		im = el(i)%elem(5)    
		ip = el(i)%elem(3)
		sig = -phy(i)%sig  !left round

	case(2) !Node 2
! 3+1----(2+1)    ^--->
! |	     |	  |   |   
! 4+1-----1+1	  ^-- v	 
		im = el(i)%elem(2)
		ip = el(i)%elem(4)
		sig = -phy(i)%sig !left round

	case(3) !Node 3
! (3+1)----2+1   <---^
! |	|	 |   |   
! 4+1-----1+1	 v---^  		
		im = el(i)%elem(5)
		ip = el(i)%elem(3)
		sig = phy(i)%sig !right round

	case(4) !Node 4
! 3+1-----2+1    <---^
! |       |      |   |   
!(4+1)----1+1	 v---^  				
		im = el(i)%elem(2)
		ip = el(i)%elem(4)
		sig = phy(i)%sig !right round
end select	
end subroutine posit


subroutine energy(dt,el,node,phy)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
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

	vol_til = abs(abl)+abs(atr) !intermedia Volume
	phy(i)%e = phy(i)%e-(phy(i)%p+phy(i)%q)*(vol_til-phy(i)%vol)/phy(i)%mas !Calculation inner energy (See SALE2D or Samarsky)
enddo
end subroutine energy


subroutine grid(dt,node,art) 
type(nodes),intent(inout) :: node(:)
type(artific),intent(inout) :: art
real(8),intent(in) :: dt
integer(4) i

do i = 1,size(node(:)) !Cycle for all nodes
node(i)%z =  node(i)%z + node(i)%u*dt*art%grid !Calculation lagrange z-postion node 
node(i)%r =  node(i)%r + node(i)%v*dt*art%grid !Calculation lagrange r-postion node 
!art%grid - coff shift vertex (if lagrange  = 1, Euler = 0, ALE = between 0-1)
enddo
end subroutine grid


subroutine advect(dt,el,node,phy,art) !Subroutine calculation advection
 type(elements),intent(inout) :: el(:)
 type(nodes),intent(inout) :: node(:)
 type(physics),intent(inout) :: phy(:)
 type(artific),intent(inout) :: art
 real(8),intent(in) :: dt
 integer(4) i,j,k,l,ii,kl
 real(8) z_p(4),r_p(4),vol,FV,a_V
 kl = 1 !For shift index in el(:)%elem
 if (art%grid<1) then !if art%grid = 1 then bypassed this block 

 	do i = 1,size(phy(:)%vol)
 		phy(i)%dum = 0d0 !Set zero flux U -momentum 
 		phy(i)%dvm = 0d0 !Set zero flux V -momentum 
 		phy(i)%mas_til = phy(i)%mas !Set intermedia mass equivalent mass
 		phy(i)%Me_til = phy(i)%mas*phy(i)%e !Set intermedia Mass*energy (inner energy)

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
 					z_p(1) = z_p(4) + node(el(i)%elem(kl+1))%u*dt*(1d0 - art%grid)   !z(final position) = z(lagrange) - u(velocity shift node in back position) *dt
 					z_p(2) = z_p(3) + node(el(i)%elem(kl+2))%u*dt*(1d0 - art%grid)
 				!r -component
 					r_p(3) = node(el(i)%elem(kl+2))%r
 					r_p(4) = node(el(i)%elem(kl+1))%r 
 					r_p(1) = r_p(4) + node(el(i)%elem(kl+1))%v*dt*(1d0 - art%grid)
 					r_p(2) = r_p(3) + node(el(i)%elem(kl+2))%v*dt*(1d0 - art%grid)
												
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
 					z_p(2) = z_p(1) + node(el(i)%elem(kl+2))%u*dt*(1d0 - art%grid)
 					z_p(3) = z_p(4) + node(el(i)%elem(kl+3))%u*dt*(1d0 - art%grid)
				
 					r_p(1) = node(el(i)%elem(kl+2))%r 
 					r_p(4) = node(el(i)%elem(kl+3))%r
 					r_p(2) = r_p(1) + node(el(i)%elem(kl+2))%v*dt*(1d0 - art%grid)
 					r_p(3) = r_p(4) + node(el(i)%elem(kl+3))%v*dt*(1d0 - art%grid)
			
 				case(3)! Left side
 				! 
 				!z_p(3).....3----2
 				!   		|    | 	   
 				!z_p(4).....4----1
 					z_p(1) = node(el(i)%elem(kl+4))%z 
 					z_p(2) = node(el(i)%elem(kl+3))%z 
 					z_p(3) = z_p(2) + node(el(i)%elem(kl+3))%u*dt*(1d0 - art%grid)
 					z_p(4) = z_p(1) + node(el(i)%elem(kl+4))%u*dt*(1d0 - art%grid)
				
 					r_p(1) = node(el(i)%elem(kl+4))%r 
 					r_p(2) = node(el(i)%elem(kl+3))%r 
 					r_p(3) = r_p(2) + node(el(i)%elem(kl+3))%v*dt*(1d0 - art%grid)
 					r_p(4) = r_p(1) + node(el(i)%elem(kl+4))%v*dt*(1d0 - art%grid)
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
 					z_p(1) = z_p(2) + node(el(i)%elem(kl+1))%u*dt*(1d0 - art%grid)
 					z_p(4) = z_p(3) + node(el(i)%elem(kl+4))%u*dt*(1d0 - art%grid)
				
 					r_p(2) = node(el(i)%elem(kl+1))%r 
 					r_p(3) = node(el(i)%elem(kl+4))%r 
 					r_p(1) = r_p(2) + node(el(i)%elem(kl+1))%v*dt*(1d0 - art%grid)
 					r_p(4) = r_p(3) + node(el(i)%elem(kl+4))%v*dt*(1d0 - art%grid)
 				case default !When we have zero index in el(:)%cont
 					exit ext	
 			end select
		
 			FV = -phy(i)%sig*volume_adv(z_p,r_p) !Calculation volume cell-shift
 			a_V = 1d0*sign(1d0,Fv) !For stability a0 = 1d0 - change variable (0 - unstable 1- more diffuse, using 0 - 1)  
 			phy(i)%mas_til = phy(i)%mas_til + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho+(1d0-a_V)*phy(i)%rho) !Calculation mass  Fv(0.5[rho_ii+rho_i]) - See Fletcher (vol. 1) or SALE2D
 			phy(i)%dum = phy(i)%dum + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%Um+(1d0-a_V)*phy(i)%rho*phy(i)%Um) !flux U -momentum 
 			phy(i)%dvm = phy(i)%dvm + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%Vm+(1d0-a_V)*phy(i)%rho*phy(i)%Vm) !flux V -momentum 
 			phy(i)%Me_til = phy(i)%Me_til + 0.5d0*Fv*((1d0+a_V)*phy(ii)%rho*phy(ii)%e+(1d0-a_V)*phy(i)%rho*phy(i)%e) !Mass*Energy 
		
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
 			!

 			node(l)%u = node(l)%u + 0.25d0*phy(i)%dum/node(l)%mas_til_v + node(l)%mas_v*node(l)%u_l/node(l)%num_cont/node(l)%mas_til_v
 			node(l)%v = node(l)%v + 0.25d0*phy(i)%dvm/node(l)%mas_til_v + node(l)%mas_v*node(l)%v_l/node(l)%num_cont/node(l)%mas_til_v
 		enddo	
 	enddo
 	node(:)%mas_v = node(:)%mas_til_v
 	do i = 1,size(phy(:)%vol)
 		phy(i)%mas = phy(i)%mas_til
 		phy(i)%e = phy(i)%Me_til/phy(i)%mas_til
 	enddo
 endif

 node(:)%u_l = node(:)%u
 node(:)%v_l = node(:)%v
contains

 real(8) function volume_adv(z,r)
 real(8) z(4),r(4)
 real(8) atr,abl

 atr = 0.5d0*((z(3)-z(2))*(r(1)-r(2))-(z(1)-z(2))*(r(3)-r(2))) !0.5*((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))
 abl = 0.5d0*((z(1)-z(4))*(r(3)-r(4))-(z(3)-z(4))*(r(1)-r(4))) !0.5*((x1-x4)*(y3-y4)-(x3-x4)*(y1-y4))
 volume_adv = abl+atr
 end function volume_adv

 end subroutine advect




subroutine artvisc(dt,node,phy)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
real(8),intent(in) :: dt
integer(4) i
real(8) Dvol

do i = 1,size(phy(:)%vol)
	Dvol = 2d0*(phy(i)%vol - phy(i)%vol_old)/(phy(i)%vol + phy(i)%vol_old)/dt !Calculation Divergence (Equivalen Hemp3D)
	phy(i)%q = min(0d0,Dvol)*(0.1d0*phy(i)%rho*Dvol*phy(i)%Vol**0.666) !Calculation Divergence (Equivalen Hemp3D)
enddo
end subroutine artvisc



end module hydro
