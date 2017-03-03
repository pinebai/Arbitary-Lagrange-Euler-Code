module eos
use object
contains
subroutine gas(i,phy)
type(physics),intent(inout) :: phy(:)
integer(4),intent(in) :: i
real(8),parameter :: gamma = 1.4  
phy(i)%p = (gamma - 1d0)*phy(i)%rho*phy(i)%e !Eos ideal gas
end subroutine gas
end module eos 
