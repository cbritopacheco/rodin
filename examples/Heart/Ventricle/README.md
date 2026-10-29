# Ventricle models

Ventricular contraction can be reduced to chamber deformation, active muscle
variables, and vascular pressures without resolving a three-dimensional wall.
For state $z$, the reduced dynamics are written as

$$
\mathcal R_{0D}(z,\dot z,t;\theta)=0,
$$

where $\theta$ contains reference geometry, material and circulation parameters.
The algebraic pressure/flow relations and valve states distinguish this system
from a prescribed pressure waveform.

[CCMLC2014_0D](CCMLC2014_0D/README.md) advances spherical chamber motion,
contractile variables, and proximal/distal arterial compliances. It records
pressure, volume, and flow without a spatial mesh. Coronary drivers use the
same model to supply pressure or motion data; the contraction-driven atrium FSI
example currently has no ventricular coupling.
