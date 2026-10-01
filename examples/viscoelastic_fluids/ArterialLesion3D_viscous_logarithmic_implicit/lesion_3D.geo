// ===========================================================================
//  lesion_3D.geo  --  idealised stenosis / fusiform aneurysm, 3D pipe
//  CFDLab USACH.  Gmsh >= 4.8, OpenCASCADE kernel.
//
//  The meridian profile of lesion_2D.geo revolved about the x axis:
//      r(x) = R [ 1 + a (1 + cos(pi x / ell)) / 2 ],   |x| <= ell
//      r(x) = R                                          otherwise
//  stenosis : a = sqrt(1 - S) - 1   (S = AREA reduction at the throat)
//  aneurysm : a = Gam - 1           (Gam = Dmax / D)
//  In 3D S is a true area reduction, unlike the planar channel of the 2D case.
//
//  Usage (all lengths scale with D; default D = 1, i.e. non-dimensional):
//    gmsh lesion_3D.geo -3 -setnumber ltype 1 -setnumber S 0.75 -o s75.msh
//
//  Physical tags:  1 inlet | 2 outlet | 3 wall (parent vessel) |
//                  5 lesion wall (|x| <= ell)  |  10 fluid
// ===========================================================================
SetFactory("OpenCASCADE");

DefineConstant[
  ltype  = {1, Choices{0="healthy", 1="stenosis", 2="fusiform aneurysm"},
            Name "1Geometry/0Lesion type"},
  S      = {0.50, Min 0, Max 0.95, Name "1Geometry/1Area reduction S (stenosis)"},
  Gam    = {1.50, Min 1, Max 3.0,  Name "1Geometry/2Dilation Dmax/D (aneurysm)"},
  ell    = {0,    Name "1Geometry/3Half-length ell/D (0 = default: 1 sten., 1.5 aneur.)"},
  Lu     = {10,   Name "1Geometry/4Inlet distance Lu/D (from lesion centre)"},
  Ld     = {25,   Name "1Geometry/5Outlet distance Ld/D (from lesion centre)"},
  D      = {1,    Name "1Geometry/6Parent diameter D"},
  Np     = {41,   Name "1Geometry/8Profile points in lesion"},
  hb     = {0.15, Name "2Mesh/0Background size hb/D"},
  hl     = {0.08, Name "2Mesh/1Lesion + wake size hl/D"},
  hw     = {0.04, Name "2Mesh/2Lesion wall size hw/D"},
  hwf    = {0.10, Name "2Mesh/3Wall size elsewhere hwf/D (Stokes layer)"},
  Lw     = {8,    Name "2Mesh/4Refined wake length Lw/D"},
  rf     = {1,    Name "2Mesh/5Refinement factor (all h divided by rf)"}
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
lcb = hb*D/rf;  lcl = hl*sc*D/rf;  lcw = hw*sc*D/rf;  lcwf = hwf*D/rf;
rmax = R; If (a > 0) rmax = R*(1 + a); EndIf

// ---- meridian half-plane (y >= 0, z = 0) -----------------------------------
pia = newp; Point(pia) = {-Lu*D, 0, 0};
piw = newp; Point(piw) = {-Lu*D, R, 0};
pw[] = {};
For i In {0:Np-1}
  x = -L + 2*L*i/(Np-1);
  pw[] += newp; Point(pw[i]) = {x, R*(1 + a*(1 + Cos(Pi*x/L))/2), 0};
EndFor
pow = newp; Point(pow) = {Ld*D, R, 0};
poa = newp; Point(poa) = {Ld*D, 0, 0};

lin  = newl; Line(lin)  = {pia, piw};
lup  = newl; Line(lup)  = {piw, pw[0]};
If (ltype == 0)
  lles = newl; Line(lles) = {pw[0], pw[Np-1]};
Else
  lles = newl; Spline(lles) = {pw[]};
EndIf
ldn  = newl; Line(ldn)  = {pw[Np-1], pow};
lout = newl; Line(lout) = {pow, poa};
lax  = newl; Line(lax)  = {poa, pia};

Curve Loop(1) = {lin, lup, lles, ldn, lout, lax};
Plane Surface(1) = {1};

// ---- revolve about the x axis ----------------------------------------------
vol[] = Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi} { Surface{1}; };

eps = 1e-6*D;
inlet[]  = Surface In BoundingBox{-Lu*D-eps, -R-eps, -R-eps, -Lu*D+eps, R+eps, R+eps};
outlet[] = Surface In BoundingBox{ Ld*D-eps, -R-eps, -R-eps,  Ld*D+eps, R+eps, R+eps};
lesion[] = Surface In BoundingBox{-L-eps, -rmax-eps, -rmax-eps, L+eps, rmax+eps, rmax+eps};
up[]     = Surface In BoundingBox{-Lu*D-eps, -R-eps, -R-eps, -L+eps, R+eps, R+eps};
dn[]     = Surface In BoundingBox{ L-eps, -R-eps, -R-eps, Ld*D+eps, R+eps, R+eps};
up[] -= {inlet[]};
dn[] -= {outlet[]};

Physical Surface("inlet", 1)  = {inlet[]};
Physical Surface("outlet", 2) = {outlet[]};
Physical Surface("wall", 3)   = {up[], dn[]};
Physical Surface("lesion", 5) = {lesion[]};
Physical Volume("fluid", 10)  = Volume{:};

// ---- size fields -------------------------------------------------------------
Field[1] = Distance;  Field[1].SurfacesList = {up[], dn[], lesion[]};
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = lcwf; Field[2].SizeMax = lcb;
Field[2].DistMin = 0;    Field[2].DistMax = 0.2*D;

Field[3] = Distance;  Field[3].SurfacesList = {lesion[]};
Field[4] = Threshold; Field[4].InField = 3;
Field[4].SizeMin = lcw;  Field[4].SizeMax = lcb;
Field[4].DistMin = 0;    Field[4].DistMax = 0.3*D;

// lesion + wake core
Field[5] = Box;
Field[5].VIn = lcl;  Field[5].VOut = lcb;
Field[5].XMin = -L - 0.5*D;  Field[5].XMax = L + Lw*D;
Field[5].YMin = -rmax;       Field[5].YMax = rmax;
Field[5].ZMin = -rmax;       Field[5].ZMax = rmax;

// inlet disk: the imposed profile is interpolated on it, and its P1 flux
// falls short of Q(t) by several percent on a background-size disk
Field[7] = Distance;  Field[7].SurfacesList = {inlet[]};
Field[8] = Threshold; Field[8].InField = 7;
Field[8].SizeMin = 0.4*lcwf; Field[8].SizeMax = lcb;
Field[8].DistMin = 0;        Field[8].DistMax = 0.5*D;

Field[6] = Min; Field[6].FieldsList = {2, 4, 5, 8};
Background Field = 6;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;
