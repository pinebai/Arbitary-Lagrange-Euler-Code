module type_boundary
use object	
implicit none
contains
subroutine boundary_flow(bou,node,phy,el)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(elements),intent(inout) :: el(:)
type(physics),intent(inout) :: phy(:)
integer(4) i,l(6),k

do i = 1,size(bou%var(:,1))
    l(:) = bou%var(i,:)
    k = bou%type_bound(l(1))
    if (k>0) then
        select case(k)
            case(1) !Full fix U and V
                call solid_wall(l,bou,node,phy,el,1)
            case(2) !Fix U
                call solid_wall(l,bou,node,phy,el,2)
            case(3) !Fix V
                call solid_wall(l,bou,node,phy,el,3)
            case(4) !Out flow
                call out_flow(l,bou,node,phy,el,4)
            case(5) !Reflection
                call reflection(l,bou,node,phy,el,5)
        end select
    endif
enddo
contains
subroutine solid_wall(l,bou,node,phy,el,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: l(6),typ
integer(4) i1,i2

call posit(l(2),i1,i2)
i1 = el(l(3))%elem(i1+1)                
i2 = el(l(3))%elem(i2+1)                 
if (typ == 1) then !If boundary = 1 fix U and V	
    node(i1)%u = 0d0
    node(i2)%u = 0d0     

    node(i1)%u_l = 0d0
    node(i2)%u_l = 0d0
    
    node(i1)%v = 0d0
    node(i2)%v = 0d0     

    node(i1)%v_l = 0d0
    node(i2)%v_l = 0d0
endif

if (typ == 2) then !If boundary = 2 fix U
    node(i1)%u = 0d0
    node(i2)%u = 0d0     

    node(i1)%u_l = 0d0
    node(i2)%u_l = 0d0
endif

if (typ == 3) then !If boundary = 3 fix V
    node(i1)%v = 0d0
    node(i2)%v = 0d0     

    node(i1)%v_l = 0d0
    node(i2)%v_l = 0d0
endif
end subroutine solid_wall

subroutine out_flow(l,bou,node,phy,el,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: l(6),typ
integer(4) i1(3),i2(3),c1,c2

c1 = l(3)
c2 = l(5)

call posit(l(2),i1(1),i2(1))
i1(1) = el(l(3))%elem(i1(1)+1)                
i2(1) = el(l(3))%elem(i2(1)+1) 

call posit(l(4),i1(2),i2(2))
i1(2) = el(l(3))%elem(i1(2)+1)                
i2(2) = el(l(3))%elem(i2(2)+1) 

call posit(l(6),i1(3),i2(3))
i1(3) = el(l(5))%elem(i1(3)+1)                
i2(3) = el(l(5))%elem(i2(3)+1) 


!node(i1(2))%u = node(i1(3))%u
!node(i2(2))%u = node(i2(3))%u

!node(i1(2))%u_l = node(i1(3))%u_l
!node(i2(2))%u_l = node(i2(3))%u_l


!node(i1(2))%v = node(i1(3))%v
!node(i2(2))%v = node(i2(3))%v

!node(i1(2))%v_l = node(i1(3))%v_l
!node(i2(2))%v_l = node(i2(3))%v_l



node(i1(1))%u = node(i1(2))%u
node(i2(1))%u = node(i2(2))%u

node(i1(1))%u_l = node(i1(2))%u_l
node(i2(1))%u_l = node(i2(2))%u_l


node(i1(1))%v = node(i1(2))%v
node(i2(1))%v = node(i2(2))%v

node(i1(1))%v_l = node(i1(2))%v_l
node(i2(1))%v_l = node(i2(2))%v_l

phy(c1)%p = phy(c2)%p  
phy(c1)%rho = phy(c2)%rho 
phy(c1)%e = phy(c2)%e
end subroutine out_flow

subroutine reflection(l,bou,node,phy,el,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: l(6),typ
integer(4) i1(3),i2(3),c1,c2

c1 = l(3)
c2 = l(5)

call posit(l(2),i1(1),i2(1))
i1(1) = el(l(3))%elem(i1(1)+1)                
i2(1) = el(l(3))%elem(i2(1)+1) 

call posit(l(4),i1(2),i2(2))
i1(2) = el(l(3))%elem(i1(2)+1)                
i2(2) = el(l(3))%elem(i2(2)+1) 

call posit(l(6),i1(3),i2(3))
i1(3) = el(l(5))%elem(i1(3)+1)                
i2(3) = el(l(5))%elem(i2(3)+1) 

node(i1(2))%v = 0d0
node(i2(2))%v = 0d0     

node(i1(2))%v_l = 0d0
node(i2(2))%v_l = 0d0

node(i1(1))%u = node(i1(2))%u
node(i2(1))%u = node(i2(2))%u

node(i1(1))%u_l = node(i1(2))%u_l
node(i2(1))%u_l = node(i2(2))%u_l

node(i1(1))%v = -node(i1(3))%v
node(i2(1))%v = -node(i2(3))%v

node(i1(1))%v_l = -node(i1(3))%v_l
node(i2(1))%v_l = -node(i2(3))%v_l

phy(c1)%p = phy(c2)%p  
!phy(c1)%rho = phy(c2)%rho 
phy(c1)%e = phy(c2)%e
end subroutine reflection

end subroutine boundary_flow

subroutine posit(num,i1,i2)
integer(4),intent(in) :: num
integer(4),intent(out) :: i1,i2

select case(num)
	case(1)
		i1 = 1
		i2 = 2
	case(2)
		i1 = 2
		i2 = 3	
	case(3)
		i1 = 3
		i2 = 4	
	case(4)
		i1 = 1
		i2 = 4
end select		
end subroutine posit
end module type_boundary
