module eos
use object
contains
subroutine gas(i,cell)
type(cells),intent(inout) :: cell(:)
integer(4),intent(in) :: i
real(8),parameter :: gamma = 1.4  
cell(i)%p = (gamma - 1d0)*cell(i)%rho*cell(i)%e !Eos ideal gas
if ((isnan(cell(i)%p).eqv..true.).or.(isnan(cell(i)%rho).eqv..true.)) then
    write(*,*) 'Error! Var NaN'
    stop
endif
end subroutine gas
end module eos 
