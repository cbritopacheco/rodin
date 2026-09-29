// ===========================================================================
//  lesion_2D.geo  --  idealised stenosis / fusiform aneurysm, 2D meridian plane
//  CFDLab USACH.  Gmsh >= 4.8, built-in kernel.
//
//  Coordinates: x = axial (z), y = radial (r).  Lesion centred at x = 0.
//  Wall profile (single family, signed amplitude a):
//      r(z) = R [ 1 + a (1 + cos(pi z / ell)) / 2 ],   |z| <= ell
//      r(z) = R                                          otherwise
//  stenosis : a = sqrt(1 - S) - 1   (S = area reduction at the throat)
//  aneurysm : a = Gam - 1           (Gam = Dmax / D)
//
//  Usage (all lengths scale with D; default D = 1, i.e. non-dimensional):
//    gmsh lesion_2D.geo -2 -setnumber ltype 1 -setnumber S 0.75 -o s75.msh
//    gmsh lesion_2D.geo -2 -setnumber ltype 2 -setnumber Gam 2.0 -o a200.msh
//
//  Physical tags:  1 inlet | 2 outlet | 3 wall (parent vessel) |
//                  4 axis (symmetry, axisymmetric case only)  |
//                  5 lesion wall (|z| <= ell)  |  10 fluid
// ===========================================================================
SetFactory("Built-in");

DefineConstant[
  ltype  = {1, Choices{0="healthy", 1="stenosis", 2="fusiform aneurysm"},
            Name "1Geometry/0Lesion type"},
  S      = {0.50, Min 0, Max 0.95, Name "1Geometry/1Area reduction S (stenosis)"},
  Gam    = {1.50, Min 1, Max 3.0,  Name "1Geometry/2Dilation Dmax/D (aneurysm)"},
  ell    = {0,    Name "1Geometry/3Half-length ell/D (0 = default: 1 sten., 1.5 aneur.)"},
  Lu     = {10,   Name "1Geometry/4Inlet distance Lu/D (from lesion centre)"},
  Ld     = {25,   Name "1Geometry/5Outlet distance Ld/D (from lesion centre)"},
  D      = {1,    Name "1Geometry/6Parent diameter D"},
  planar = {0, Choices{0="axisymmetric half-plane", 1="full planar channel"},
            Name "1Geometry/7Domain"},
  Np     = {101,  Name "1Geometry/8Profile points in lesion"},
  hb     = {0.05, Name "2Mesh/0Background size hb/D"},
  hl     = {0.025,Name "2Mesh/1Lesion + wake size hl/D"},
  hw     = {0.006,Name "2Mesh/2Wall size in lesion + wake hw/D"},
  hwf    = {0.012,Name "2Mesh/3Wall size elsewhere hwf/D (Stokes layer)"},
  Lw     = {12,   Name "2Mesh/4Refined wake length Lw/D"},
  rf     = {1,    Name "2Mesh/5Refinement factor (all h divided by rf)"},
  order  = {1,    Name "2Mesh/6Element order"}
];

// ---- geometry parameters ---------------------------------------------------
R = D/2;
a = 0; ellD = 1;
If (ltype == 1) a = Sqrt(1 - S) - 1; ellD = 1.0; EndIf
If (ltype == 2) a = Gam - 1;         ellD = 1.5; EndIf
If (ell > 0) ellD = ell; EndIf
L = ellD*D;

// sizes scale with the throat radius so that the throat is always resolved
sc = 1; If (a < 0) sc = 1 + a; EndIf
lcb = hb*D/rf;  lcl = hl*sc*D/rf;  lcw = hw*sc*D/rf;
rmax = R; If (a > 0) rmax = R*(1 + a); EndIf

// ---- points ----------------------------------------------------------------
p1 = newp; Point(p1) = {-Lu*D, 0, 0, lcb};
p2 = newp; Point(p2) = { Ld*D, 0, 0, lcb};
p3 = newp; Point(p3) = { Ld*D, R, 0, lcb};
pw[] = {};
For k In {0:Np-1}
  z  = L - 2*L*k/(Np-1);
  r  = R*(1 + a*(1 + Cos(Pi*z/L))/2);
  pp = newp; Point(pp) = {z, r, 0, lcl};
  pw[] += {pp};
EndFor
p4 = newp; Point(p4) = {-Lu*D, R, 0, lcb};

// ---- curves and surface ----------------------------------------------------
cA = newl; Line(cA)   = {p1, p2};
cO = newl; Line(cO)   = {p2, p3};
cD = newl; Line(cD)   = {p3, pw[0]};
cL = newl; Spline(cL) = {pw[]};
cU = newl; Line(cU)   = {pw[Np-1], p4};
cI = newl; Line(cI)   = {p4, p1};
Curve Loop(1) = {cA, cO, cD, cL, cU, cI};
Plane Surface(1) = {1};

If (planar == 1)
  Symmetry {0, 1, 0, 0} { Duplicata { Surface{1}; } }
  Coherence;
EndIf

// ---- classify boundary curves by bounding box ------------------------------
eps = 1e-6*D;
cin[] = {}; cout[] = {}; cwall[] = {}; caxis[] = {}; cles[] = {};
allc[] = Curve "*";
For i In {0:#allc[]-1}
  c = allc[i];
  bb[] = BoundingBox Curve{c};           // xmin ymin zmin xmax ymax zmax
  If (bb[3] - bb[0] < eps)
    If (Fabs(bb[0] + Lu*D) < eps) cin[] += {c}; EndIf
    If (Fabs(bb[0] - Ld*D) < eps) cout[] += {c}; EndIf
  ElseIf (Fabs(bb[1]) < eps && Fabs(bb[4]) < eps)
    caxis[] += {c};                      // internal in the planar case
  ElseIf (bb[0] > -L - eps && bb[3] < L + eps)
    cles[] += {c};
  Else
    cwall[] += {c};
  EndIf
EndFor

Physical Curve("inlet", 1)  = {cin[]};
Physical Curve("outlet", 2) = {cout[]};
Physical Curve("wall", 3)   = {cwall[]};
If (planar == 0) Physical Curve("axis", 4) = {caxis[]}; EndIf
Physical Curve("lesion", 5) = {cles[]};
Physical Surface("fluid", 10) = {Surface "*"};

// ---- mesh size fields ------------------------------------------------------
Field[1] = Distance;  Field[1].CurvesList = {cwall[], cles[]};  Field[1].Sampling = 2000;
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = lcw;  Field[2].SizeMax = lcb;
Field[2].DistMin = 0.01*D;  Field[2].DistMax = 0.20*D;
Field[3] = Box;   // near-wall refinement only around lesion and wake
Field[3].VIn = 0;  Field[3].VOut = 1e22;
Field[3].XMin = -L - 2*D;  Field[3].XMax = L + Lw*D;
Field[3].YMin = -2*rmax;   Field[3].YMax = 2*rmax;
Field[4] = Max;   Field[4].FieldsList = {2, 3};
Field[5] = Box;   // lesion + wake core refinement
Field[5].VIn = lcl;  Field[5].VOut = lcb;
Field[5].XMin = -L - 1*D;  Field[5].XMax = L + Lw*D;
Field[5].YMin = -2*rmax;   Field[5].YMax = 2*rmax;
Field[5].Thickness = 2*D;
Field[7] = Threshold; Field[7].InField = 1;   // whole wall: Stokes layer at high Wo
Field[7].SizeMin = hwf*D/rf;  Field[7].SizeMax = lcb;
Field[7].DistMin = 0.01*D;    Field[7].DistMax = 0.15*D;
Field[6] = Min;   Field[6].FieldsList = {4, 5, 7};
Background Field = 6;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;              // Frontal-Delaunay
Mesh.ElementOrder = order;
Mesh.HighOrderOptimize = (order > 1) ? 2 : 0;
Mesh.MshFileVersion = 2.2;       // ASCII v2.2: read by MFEM/Rodin, FEniCS, OpenFOAM
