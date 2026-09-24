// =====================================================================
//  Diverging microchannel — 3D FLUID domain (tetrahedra)
//  Ref.: Alugoju, Dubey & Javed, IJHMT 160 (2020) 120212 — Table 1, Case 1
//        Wi = 0.2 mm, Wo = 0.4 mm, L = 30 mm, d = 0.4 mm (delta = 2, theta = 0.38 deg)
//
//  x: flow direction (0 -> L), y: width (symmetric about y = 0), z: depth (0 = heated base)
//
//  Export for Rodin (MEDIT, Dimension 3):
//     gmsh DivergingMicrochannel3D.geo -3 -format mesh -o DivergingMicrochannel3D.mesh
//
//  Physical tags (-> MEDIT references / Rodin attributes):
//     Volume : 1  (fluid)
//     Faces  : 1 = Inlet      (x = 0)
//              2 = Outlet     (x = L)
//              3 = Bottom     (z = 0, heated wall: q'' here)
//              4 = Top        (z = d, cover)
//              5 = Side walls (y = +/- W(x)/2)
// =====================================================================

SetFactory("Built-in");

// ---------------- Parameters (mm) ----------------
DefineConstant[
  L    = {30.0, Name "Geometry/L (mm)"},
  Wi   = {0.2,  Name "Geometry/Inlet width Wi (mm)"},
  Wo   = {0.4,  Name "Geometry/Outlet width Wo (mm)"},
  d    = {0.4,  Name "Geometry/Depth d (mm)"},
  unit = {1e-3, Name "Geometry/Scale (1e-3 -> m, 1 -> mm)"},
  structured = {1, Choices{0="Unstructured (Delaunay)", 1="Transfinite (hexes split into tets)"},
                Name "Mesh/Type"},
  nW   = {6,    Name "Mesh/Elements across width"},
  nD   = {10,   Name "Mesh/Elements across depth"},
  dx   = {0.115,  Name "Mesh/Streamwise size dx (mm, structured)"},
  hU   = {0.04, Name "Mesh/Size h (mm, unstructured)"}
];

hp = (structured == 1) ? dx : hU;

// ---------------- Points ----------------
Point(1) = {0,       -Wi/2*unit, 0,      hp*unit};
Point(2) = {L*unit,  -Wo/2*unit, 0,      hp*unit};
Point(3) = {L*unit,   Wo/2*unit, 0,      hp*unit};
Point(4) = {0,        Wi/2*unit, 0,      hp*unit};
Point(5) = {0,       -Wi/2*unit, d*unit, hp*unit};
Point(6) = {L*unit,  -Wo/2*unit, d*unit, hp*unit};
Point(7) = {L*unit,   Wo/2*unit, d*unit, hp*unit};
Point(8) = {0,        Wi/2*unit, d*unit, hp*unit};

// ---------------- Lines ----------------
// streamwise
Line(1)  = {1, 2};  Line(2)  = {4, 3};  Line(3)  = {5, 6};  Line(4)  = {8, 7};
// spanwise
Line(5)  = {1, 4};  Line(6)  = {2, 3};  Line(7)  = {5, 8};  Line(8)  = {6, 7};
// vertical
Line(9)  = {1, 5};  Line(10) = {2, 6};  Line(11) = {3, 7};  Line(12) = {4, 8};

// ---------------- Faces ----------------
Curve Loop(1) = {5, 12, -7, -9};    Plane Surface(1) = {1};   // inlet  x=0
Curve Loop(2) = {6, 11, -8, -10};   Plane Surface(2) = {2};   // outlet x=L
Curve Loop(3) = {1, 6, -2, -5};     Plane Surface(3) = {3};   // bottom z=0
Curve Loop(4) = {3, 8, -4, -7};     Plane Surface(4) = {4};   // top    z=d
Curve Loop(5) = {1, 10, -3, -9};    Plane Surface(5) = {5};   // side y<0
Curve Loop(6) = {2, 11, -4, -12};   Plane Surface(6) = {6};   // side y>0

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// ---------------- Mesh ----------------
If (structured == 1)
  Nx = Ceil(L / dx) + 1;
  Transfinite Curve{1, 2, 3, 4}   = Nx;
  Transfinite Curve{5, 6, 7, 8}   = nW + 1;
  Transfinite Curve{9, 10, 11, 12} = nD + 1;
  Transfinite Surface{1, 2, 3, 4, 5, 6};
  Transfinite Volume{1} = {1, 2, 3, 4, 5, 6, 7, 8};
Else
  Mesh.Algorithm   = 6;   // Frontal-Delaunay (2D faces)
  Mesh.Algorithm3D = 1;   // Delaunay
  Mesh.OptimizeNetgen = 1;
EndIf

Mesh.ElementOrder = 1;
Mesh.SaveAll = 0;         // export only physical entities
Mesh.SaveElementTagType = 2;   // MEDIT refs = PHYSICAL tags (default would be elementary)

// ---------------- Physical groups ----------------
Physical Volume("Fluid", 1)       = {1};
Physical Surface("Inlet", 1)      = {1};
Physical Surface("Outlet", 2)     = {2};
Physical Surface("Bottom", 3)     = {3};
Physical Surface("Top", 4)        = {4};
Physical Surface("SideWalls", 5)  = {5, 6};
