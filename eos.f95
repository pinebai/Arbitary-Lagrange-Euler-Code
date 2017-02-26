module eos
use object
contains
subroutine gas(phy)
type(physics),intent(inout) :: phy(:)
integer(4) i
real(8),parameter :: gamma = 1.4  
do i = 1,size(phy(:))
	phy(i)%p = (gamma - 1d0)*phy(i)%rho*phy(i)%e
enddo
end subroutine gas
end module eos 
