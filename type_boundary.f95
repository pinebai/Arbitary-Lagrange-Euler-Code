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
    ! l = (number boundary, bound_1_cell_1, cell_1, bound_2_cell_1, cell_2, bound_cell_1)
    k = bou%type_bound(l(1))
    if (k.ne.0) then
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
           ! case(6) !In flow
            !    call in_flow(l,bou,node,phy,el,6)
                
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
i1 = el(l(3))%elem(i1)                
i2 = el(l(3))%elem(i2) 

if (typ == 1) then !If boundary = 1 fix U and V	
    node(i1)%u = 0d0
    node(i2)%u = 0d0     

    node(i1)%u_l = 0d0
    node(i2)%u_l = 0d0
    
    node(i1)%v = 0d0
    node(i2)%v = 0d0     

    node(i1)%v_l = 0d0
    node(i2)%v_l = 0d0
    
    node(i1)%mark = 1
    node(i2)%mark = 1
endif

if (typ == 2) then !If boundary = 2 fix U
    node(i1)%u = 0d0
    node(i2)%u = 0d0     

    node(i1)%u_l = 0d0
    node(i2)%u_l = 0d0

    node(i1)%mark = 2
    node(i2)%mark = 2
endif

if (typ == 3) then !If boundary = 3 fix V
    node(i1)%v = 0d0
    node(i2)%v = 0d0     

    node(i1)%v_l = 0d0
    node(i2)%v_l = 0d0

    node(i1)%mark = 3
    node(i2)%mark = 3
endif
end subroutine solid_wall

! subroutine in_flow(l,bou,node,phy,el,typ)
! type(boundary),intent(inout) :: bou
! type(nodes),intent(inout) :: node(:)
! type(physics),intent(inout) :: phy(:)
! type(elements),intent(inout) :: el(:)
! integer(4),intent(in) :: l(6),typ
! integer(4) i1(3),i2(3),el1,el2
! 
! el1 = l(3) !Set number cell_1
! el2 = l(5) !Set number cell_2
! call posit(l(2),i1(1),i2(1)) !Find index node for element
! i1(1) = el(l(3))%elem(i1(1)) !node_1 boundary cell_1                
! i2(1) = el(l(3))%elem(i2(1)) !node_2 boundary cell_1
! 
! call posit(l(4),i1(2),i2(2)) !Find index node for element
! i1(2) = el(l(3))%elem(i1(2)) !node_3 boundary cell_1                 
! i2(2) = el(l(3))%elem(i2(2)) !node_4 boundary cell_1   
! 
! if ((node(i1(1))%mark.ne.1).and.(node(i1(1))%mark.ne.2)) then
!     node(i1(1))%u = 0.75
!     node(i1(1))%u_l = 0.75
! endif
! 
! if ((node(i1(1))%mark.ne.1).and.(node(i1(1))%mark.ne.3)) then
!     node(i1(1))%v = 0.75
!     node(i1(1))%v_l = 0.75
! endif
! 
! if ((node(i2(1))%mark.ne.1).and.(node(i2(1))%mark.ne.2)) then
!     node(i2(1))%u = 0.75
!     node(i2(1))%u_l = 0.75
! endif
! 
! if ((node(i2(1))%mark.ne.1).and.(node(i2(1))%mark.ne.3)) then
!     node(i2(1))%v = 0.75
!     node(i2(1))%v_l = 0.75
! endif
! 
! 
! 
! 
! end subroutine in_flow


subroutine out_flow(l,bou,node,phy,el,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: l(6),typ
integer(4) i1(3),i2(3),el1,el2

el1 = l(3) !Set number cell_1
el2 = l(5) !Set number cell_2

call posit(l(2),i1(1),i2(1)) !Find index node for element
i1(1) = el(l(3))%elem(i1(1)) !node_1 boundary cell_1                
i2(1) = el(l(3))%elem(i2(1)) !node_2 boundary cell_1

call posit(l(4),i1(2),i2(2)) !Find index node for element
i1(2) = el(l(3))%elem(i1(2)) !node_3 boundary cell_1                 
i2(2) = el(l(3))%elem(i2(2)) !node_4 boundary cell_1   

if ((node(i1(1))%mark.ne.1).and.(node(i1(1))%mark.ne.2)) then
    node(i1(1))%u = node(i1(2))%u 
    node(i1(1))%u_l = node(i1(2))%u_l
endif

if ((node(i1(1))%mark.ne.1).and.(node(i1(1))%mark.ne.3)) then
    node(i1(1))%v = node(i1(2))%v
    node(i1(1))%v_l = node(i1(2))%v_l
endif

if ((node(i2(1))%mark.ne.1).and.(node(i2(1))%mark.ne.2)) then
    node(i2(1))%u = node(i2(2))%u
    node(i2(1))%u_l = node(i2(2))%u_l
endif

if ((node(i2(1))%mark.ne.1).and.(node(i2(1))%mark.ne.3)) then
    node(i2(1))%v = node(i2(2))%v
    node(i2(1))%v_l = node(i2(2))%v_l
endif

phy(el1)%p = phy(el2)%p  
phy(el1)%rho = phy(el2)%rho 
phy(el1)%e = phy(el2)%e
end subroutine out_flow

subroutine reflection(l,bou,node,phy,el,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
type(elements),intent(inout) :: el(:)
integer(4),intent(in) :: l(6),typ
integer(4) i1(3),i2(3),el1,el2

el1 = l(3)
el2 = l(5)
call posit(l(2),i1(1),i2(1))
i1(1) = el(l(3))%elem(i1(1))                
i2(1) = el(l(3))%elem(i2(1)) 

call posit(l(4),i1(2),i2(2))
i1(2) = el(l(3))%elem(i1(2))                
i2(2) = el(l(3))%elem(i2(2)) 

call posit(l(6),i1(3),i2(3)) !Find node in cell_2
i1(3) = el(l(5))%elem(i1(3)) !Set number node_1 for cell_2                
i2(3) = el(l(5))%elem(i2(3)) !Set number node_2 for cell_2

node(i1(2))%v = 0d0 !Set V-velocity zero
node(i1(2))%v_l = 0d0

if ((node(i1(1))%mark.ne.1).and.(node(i1(1))%mark.ne.2)) then
    node(i1(1))%u = node(i1(2))%u !Set u-velocity equivalent
    node(i1(1))%u_l = node(i1(2))%u_l
endif

if ((node(i1(1))%mark.ne.1)) then!.or.(node(i1(1))%mark.ne.3)) then    
    node(i1(1))%v = -node(i1(3))%v !Set V-velocity inverse equivalent between node_1 cell_1 and node_1 cell_2
    node(i1(1))%v_l = -node(i1(3))%v_l
endif

node(i2(2))%v = 0d0 !Set V-velocity zero   
node(i2(2))%v_l = 0d0

if ((node(i2(1))%mark.ne.1).and.(node(i2(1))%mark.ne.2)) then
    node(i2(1))%u = node(i2(2))%u
    node(i2(1))%u_l = node(i2(2))%u_l
endif
if ((node(i2(1))%mark.ne.1)) then!  .or.(node(i2(1))%mark.ne.3)) then
    node(i2(1))%v = -node(i2(3))%v
    node(i2(1))%v_l = -node(i2(3))%v_l
endif

phy(el1)%p = phy(el2)%p  
phy(el1)%rho = phy(el2)%rho 
phy(el1)%e = phy(el2)%e
end subroutine reflection

end subroutine boundary_flow

subroutine posit(num,i1,i2)
integer(4),intent(in) :: num
integer(4),intent(out) :: i1,i2
!     2  
!   3---2
! 3 |   | 1
!   4---1
!     4
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

