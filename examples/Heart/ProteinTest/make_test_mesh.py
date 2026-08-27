#!/usr/bin/env python3
"""Genera la malla de prueba: tubo recto con un ensanchamiento (segmento
ectasico). Una entrada, una salida, una pared. El ensanchamiento es lo que
produce recirculacion y tiempo de residencia largo; sin el, el modelo de
coagulacion devuelve cero y el test no ensena nada."""
import numpy as np, sys

R0   = 1.5e-3     # radio del conducto (m)
L    = 40.0e-3    # longitud (m)
AMP  = 2.0        # amplitud del ensanchamiento: R_max = R0*(1+AMP)
ZC   = 0.5*L      # centro del ensanchamiento
W    = 4.0e-3     # anchura gaussiana
NT, NR, NZ = 24, 5, 80
ATTR_VOL, ATTR_WALL, ATTR_IN, ATTR_OUT = 1, 2, 4, 7

def radius(z): return R0*(1.0 + AMP*np.exp(-((z-ZC)/W)**2))

zs = np.linspace(0.0, L, NZ+1)
V, sec = [], []
for z in zs:
    Rz = radius(z); idx = [len(V)]; V.append((0.0,0.0,z))
    for j in range(1, NR+1):
        r = Rz*j/NR
        for k in range(NT):
            a = 2*np.pi*k/NT
            idx.append(len(V)); V.append((r*np.cos(a), r*np.sin(a), z))
    sec.append(idx)
V = np.array(V)

def node(s,j,k):           # j=0 centro, j=1..NR anillos
    return sec[s][0] if j==0 else sec[s][1+(j-1)*NT+(k%NT)]

tris=[]                    # triangulacion de la seccion (indices locales j,k)
for k in range(NT): tris.append(((0,0),(1,k),(1,k+1)))
for j in range(1,NR):
    for k in range(NT):
        tris.append(((j,k),(j+1,k),(j+1,k+1)))
        tris.append(((j,k),(j+1,k+1),(j,k+1)))

def vol(a,b,c,d):
    return np.dot(np.cross(V[b]-V[a], V[c]-V[a]), V[d]-V[a])/6.0

tets=[]
for s in range(NZ):
    for (p,q,r) in tris:
        a1,b1,c1 = node(s,*p),   node(s,*q),   node(s,*r)
        a2,b2,c2 = node(s+1,*p), node(s+1,*q), node(s+1,*r)
        for t in [(a1,b1,c1,c2),(a1,b1,c2,b2),(a1,b2,c2,a2)]:
            if vol(*t) < 0: t = (t[1],t[0],t[2],t[3])
            if abs(vol(*t)) > 1e-18: tets.append(t)

faces=[]
for (p,q,r) in tris:                      # entrada (z=0) y salida (z=L)
    faces.append((node(0,*p),  node(0,*r),  node(0,*q),  ATTR_IN))
    faces.append((node(NZ,*p), node(NZ,*q), node(NZ,*r), ATTR_OUT))
for s in range(NZ):                        # pared lateral
    for k in range(NT):
        a1,b1 = node(s,NR,k),   node(s,NR,k+1)
        a2,b2 = node(s+1,NR,k), node(s+1,NR,k+1)
        faces.append((a1,b1,b2,ATTR_WALL)); faces.append((a1,b2,a2,ATTR_WALL))

out = sys.argv[1] if len(sys.argv)>1 else 'ProteinTest.mesh'
with open(out,'w') as f:
    f.write("MeshVersionFormatted 2\nDimension\n3\nVertices\n%d\n"%len(V))
    for v in V: f.write("%.12g %.12g %.12g 0\n"%tuple(v))
    f.write("Tetrahedra\n%d\n"%len(tets))
    for t in tets: f.write("%d %d %d %d %d\n"%(t[0]+1,t[1]+1,t[2]+1,t[3]+1,ATTR_VOL))
    f.write("Triangles\n%d\n"%len(faces))
    for t in faces: f.write("%d %d %d %d\n"%(t[0]+1,t[1]+1,t[2]+1,t[3]))
    f.write("End\n")

Vt = sum(abs(vol(*t)) for t in tets)
Vex = np.trapezoid(np.pi*radius(zs)**2, zs)
print("vertices %d  tetraedros %d  caras %d"%(len(V),len(tets),len(faces)))
print("volumen malla %.4e m3   exacto %.4e m3   error %.2f%%"%(Vt,Vex,100*abs(Vt/Vex-1)))
print("R0=%.2f mm  Rmax=%.2f mm  razon de expansion %.1f"%(R0*1e3,radius(ZC)*1e3,radius(ZC)/R0))
