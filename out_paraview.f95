module out_paraview
use object
contains

subroutine out(kl,slide,t,el,node,cell)
type(elements),intent(inout) :: el(:)
type(nodes),intent(inout) :: node(:)
type(cells),intent(inout) :: cell(:)
integer(4),intent(in) :: kl,slide
real(8),intent(in) :: t
character(8) cha
integer(4) i

if (mod(kl,slide) == 0) then
    
    write(cha, '(I8)') kl+10000000
    open(1,file = 'result/'//cha(2:)//'.vtu')
    write(1,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
    write(1,*) '<UnstructuredGrid>'

    write(1,*) '<Piece NumberOfPoints="',size(node(:)),'" NumberOfCells="',size(cell(:)),'">'
    write(1,*) '<Points>'
    write(1,*) '<DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
    do i = 1,size(node(:))
        write(1,*) node(i)%z,node(i)%r,0.0
    enddo
    write(1,*) '</DataArray>'
    write(1,*) '</Points>'
    write(1,*) '<CellData Scalars="scalars">'
    
    !rotor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Rotor velocity" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%rot
    enddo
    write(1,*) '</DataArray>'
    
    !divergence!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Divergence velocity" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%diver
    enddo
    write(1,*) '</DataArray>'
    
    !massa!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Massa" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%mas
    enddo
    write(1,*) '</DataArray>'
    !Density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Density" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%rho
    enddo
    write(1,*) '</DataArray>'
    !Volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Volume" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%vol
    enddo
    write(1,*) '</DataArray>'
    !Preasure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Preassure" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%p
    enddo
    write(1,*) '</DataArray>'
    !Energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Inner energy" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%e
    enddo
    write(1,*) '</DataArray>'
    !Total Energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '<DataArray type="Float64" Name="Total energy" format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) cell(i)%e
    enddo
    write(1,*) '</DataArray>'    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(1,*) '</CellData>'
    write(1,*) '<Cells>'
    write(1,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) el(i)%elem(1)-1,el(i)%elem(2)-1,el(i)%elem(3)-1,el(i)%elem(4)-1
    enddo
    write(1,*) '</DataArray>'
    write(1,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) i*4
    enddo
    write(1,*) '</DataArray>'
    write(1,*) '<DataArray type="Int32" Name="types" Format="ascii">'
    do i = 1,size(cell(:))
    write(1,*) '9'  
    enddo
    write(1,*) '</DataArray>'
    write(1,*) '</Cells>'
    write(1,*) '</Piece>'
    write(1,*) '</UnstructuredGrid>'
    write(1,*) '</VTKFile>'
    close(1)
endif    
end subroutine out 
end module out_paraview
