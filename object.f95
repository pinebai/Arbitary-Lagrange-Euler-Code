module object
type elements
integer(4) elem(5) !Elements
integer(4) cont(9) !Contact elements : el1 el2 el3 el4 side1 side2 side3 side4
real(8) rad !radial
end type elements

type physics
real(8) sig !Sign countours 
real(8) vol !volume
real(8) rho !density
real(8) mas !massa
real(8) p !preasure
real(8) e !inner energy
real(8) Etot !Total energy

real(8) q !artvisc
real(8) vol_old !old volume

real(8) mas_til !massa for advect
real(8) Me_til !Mass*Energy for advect
real(8) um !sum velocity in nodes on center cell (U-componet)
real(8) vm !sum velocity in nodes on center cell (V-componet)
real(8) dum !flux impulse U 
real(8) dvm !flux impulse V

real(8) uz,ur,vz,vr,diver,rot
end type physics

type nodes
real(8) z,r !coordinates
real(8) u,v !velocity u and v
real(8) u_l,v_l !old velocity
real(8) mas_v !vertex mass in node
real(8) mas_til_v !vertex mass in node for advect
real(8) num_cont !number cell near node
integer(4) mark !marker nodes for boundary
end type nodes

type numerical 
real(8) art !Artvisc
real(8) grid  !Proportional between lagrange grid and euler grid (1 - full lagrange, 0 - full euler, other - ALE)
real(8) a0 !For Euler a0
real(8) an,ksi !For relax grid
end type numerical

type boundary
integer(4),allocatable :: num_bound(:),type_bound(:) !Number bound and type_bound (1- fix, 2- fix U- component, 3- fix V-component)
integer(4),allocatable :: var(:,:) !var(:,1) - number boundary, var(:,2) - first node, var(:,3) - second node  
end type boundary

end module object
