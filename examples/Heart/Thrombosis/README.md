# Thrombosis models

A prescribed or computed velocity can transport coagulation species while
local reactions and wall fluxes control conversion. With concentration $c_i$,
diffusivity $D_i$, and reaction $R_i$, the one-way model has the form

$$
\partial_t c_i+u\cdot\nabla c_i-D_i\Delta c_i=R_i(c).
$$

The concentrations do not alter fluid stress, permeability, or geometry. This
assumption restricts interpretation to transport and initiation before an
accumulated clot materially changes the flow. Wall-derived production enters
as a boundary flux, not as a volumetric source.

[EctaticTube](EctaticTube/README.md) combines fixed-wall pulsatile flow with
thrombin/fibrinogen/fibrin transport and a washout tracer. The tracer separates
retention of the original bolus from local reaction and production.
