# Local qualification record

1800-second watchdogs. The complete six-case finite matrix is locally verified:
216 process groups select 432 configurations with 864 successful rank-local
reports. Both thread configurations were compiled freshly for native and
real-PETSc execution. Source identity, dependency freshness, exact case
selection, runtime results and resource guards were audited independently.
The retained dependency libraries and finite mesh sequences do not constitute
a complete dependency-library rebuild, hosted-CI certification or a uniform
stability proof.

## Retained historical observations

An isolated OpenMP two-rank hexahedral P6/P5
patch passed with sampled aggregate participant RSS of approximately
1.01 GiB; the earlier combined process exceeded the 3 GiB validation limit
after preceding rate studies. This comparison motivates fresh-process
validation but does not identify the allocation or retention mechanism.
The current-source combined rate process additionally exceeded the same
memory guard at OpenMP MPI rank four on hexahedra. An isolated full P3/P2
hierarchy passed with sampled aggregate participant RSS of approximately
2.57 GiB, retaining all three levels and the same solver. This supports
separating the two rate processes without reducing their mathematical scope;
it does not identify a particular retained allocation. Earlier complete
matrix records are preserved, but current-source verification of the revised
process partition remains a separate gate and is not yet certified. Hosted
CI certification is also separate. Sampled RSS is a resource-guard observation,
not a continuous peak-memory measurement or a performance benchmark.

Fresh focused
quadratic-geometry verification passed all 60 process groups, 120 configurations
and 264 rank-local reports across the six geometries, Local/MPI
rank-one-through-four contexts and sequential/OpenMP builds. The uncorrected
OpenMP tetrahedral control failed the stricter pressure regression in all five
contexts. This is a fixed-mesh forward-accuracy certificate, not renewed
certification of the complete higher-order rate matrix or hosted PETSc 3.19 CI.
