module out_paraview
use object
contains

subroutine out(t,el,node,phy)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(physics),intent(inout) :: phy(:)
real(8),intent(in) :: t
integer(4) i

open(1,file = 'res.vtu')
write(1,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(1,*) '<UnstructuredGrid>'

write(1,*) '<Piece NumberOfPoints="',size(node(:)),'" NumberOfCells="',size(phy(:)),'">'
write(1,*) '<Points>'
write(1,*) '<DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
do i = 1,size(node(:))
	write(1,*) node(i)%z,node(i)%r,0.0
enddo
write(1,*) '</DataArray>'
write(1,*) '</Points>'
write(1,*) '<CellData Scalars="scalars">'
!Density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '<DataArray type="Float64" Name="mas" format="ascii">'
do i = 1,size(phy(:))
write(1,*) phy(i)%mas
enddo
write(1,*) '</DataArray>'
!Density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '<DataArray type="Float64" Name="rho" format="ascii">'
do i = 1,size(phy(:))
write(1,*) phy(i)%rho
enddo
write(1,*) '</DataArray>'
!Volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '<DataArray type="Float64" Name="vol" format="ascii">'
do i = 1,size(phy(:))
write(1,*) phy(i)%vol
enddo
write(1,*) '</DataArray>'
!Preasure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '<DataArray type="Float64" Name="P" format="ascii">'
do i = 1,size(phy(:))
write(1,*) phy(i)%p
enddo
write(1,*) '</DataArray>'
!Energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '<DataArray type="Float64" Name="e" format="ascii">'
do i = 1,size(phy(:))
write(1,*) phy(i)%e
enddo
write(1,*) '</DataArray>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(1,*) '</CellData>'
write(1,*) '<Cells>'
write(1,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do i = 1,size(phy(:))
write(1,*) el(i)%elem(2)-1,el(i)%elem(3)-1,el(i)%elem(4)-1,el(i)%elem(5)-1
enddo
write(1,*) '</DataArray>'
write(1,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
do i = 1,size(phy(:))
write(1,*) i*4
enddo
write(1,*) '</DataArray>'
write(1,*) '<DataArray type="Int32" Name="types" Format="ascii">'
do i = 1,size(phy(:))
write(1,*) '9'  
enddo
write(1,*) '</DataArray>'
write(1,*) '</Cells>'
write(1,*) '</Piece>'
write(1,*) '</UnstructuredGrid>'
write(1,*) '</VTKFile>'
close(1)
end subroutine out 
end module out_paraview
