module type_boundary
use object	
implicit none
contains
subroutine boundary_flow(bou,node)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
integer(4) i,l(7)

do i = 1,size(bou%var(:,1))
    l(:) = bou%var(i,:)

    select case(bou%type_bound(l(1)))
        case(1) !Full fix U and V
            call solid_wall(l,bou,node,1)
        case(2) !Fix U
            call solid_wall(l,bou,node,2)
        case(3) !Fix V
            call solid_wall(l,bou,node,3)
        case(4) !Out flow
        case(5) !Reflection
    end select
    
enddo
contains
subroutine solid_wall(l,bou,node,typ)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)
integer(4),intent(in) :: l(7),typ

                
if (typ == 1) then !If boundary = 1 fix U and V
    node(l(2))%u = 0d0
    node(l(3))%u = 0d0     

    node(l(2))%u_l = 0d0
    node(l(3))%u_l = 0d0
    
    node(l(2))%v = 0d0
    node(l(3))%v = 0d0     

    node(l(2))%v_l = 0d0
    node(l(3))%v_l = 0d0
endif

if (typ == 2) then !If boundary = 2 fix U
    node(l(2))%u = 0d0
    node(l(3))%u = 0d0     

    node(l(2))%u_l = 0d0
    node(l(3))%u_l = 0d0
    
endif

if (typ == 3) then !If boundary = 3 fix V
    node(l(2))%v = 0d0
    node(l(3))%v = 0d0     

    node(l(2))%v_l = 0d0
    node(l(3))%v_l = 0d0
endif
    
end subroutine solid_wall

subroutine out_flow()
end subroutine out_flow

subroutine reflection()
end subroutine reflection

end subroutine boundary_flow

end module type_boundary
