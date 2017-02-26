module type_boundary
use object	
implicit none
contains
subroutine boundary_vel(bou,node)
type(boundary),intent(inout) :: bou
type(nodes),intent(inout) :: node(:)


integer(4) i,j,l(3)

do i = 1,size(bou%var(:,1))
    l(1) = bou%var(i,1)
    l(2) = bou%var(i,2)
    l(3) = bou%var(i,3)
    
    do j = 1,size(bou%num_bound(:)) !cycle for all boundary 
                
        if (l(1) == bou%num_bound(j)) then  !If var(1) equivalent number boundary
                    
            if (bou%type_bound(j) == 1) then !If boundary = 1 fix U and V
                node(l(2))%u = 0d0
                node(l(3))%u = 0d0     

                node(l(2))%u_l = 0d0
                node(l(3))%u_l = 0d0
                
                node(l(2))%v = 0d0
                node(l(3))%v = 0d0     

                node(l(2))%v_l = 0d0
                node(l(3))%v_l = 0d0
            endif
            
            if (bou%type_bound(j) == 2) then !If boundary = 2 fix U
                node(l(2))%u = 0d0
                node(l(3))%u = 0d0     

                node(l(2))%u_l = 0d0
                node(l(3))%u_l = 0d0
                
            endif
        
            if (bou%type_bound(j) == 3) then !If boundary = 3 fix V
                node(l(2))%v = 0d0
                node(l(3))%v = 0d0     

                node(l(2))%v_l = 0d0
                node(l(3))%v_l = 0d0
            endif
            
        endif
    enddo
enddo
end subroutine boundary_vel

end module type_boundary
