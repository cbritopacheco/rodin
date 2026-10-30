# Benchmark qualification and measurement procedures


Convergence certification belongs to [the convergence module](../convergence/README.md).
Benchmark guards establish that a measurement still evaluates its stated
discrete operation; they do not certify an approximation rate. Unit,
manufactured and convergence regressions remain the correctness gates.
Benchmark execution contributes to source coverage only when it is explicitly
instrumented and included in the coverage collection; registration alone does
not contribute coverage.

Qualification uses applicable independent identities, original-entry-loop
comparisons, manufactured residuals and coefficient-state changes. Shared
reference dependencies are stated: an entrywise oracle sharing production
mapping does not independently validate that map. Reference memory must follow
the intended cell lifecycle rather than changing the timed working set.

For a discrete matrix comparison, a typical normalized discrepancy is

$$
e_A=\frac{\|A_h-A_{\mathrm{ref}}\|_F}
{\max(1,\|A_{\mathrm{ref}}\|_F)}.
$$

This normalization presumes the workload's nondimensional numerical scale.
For dimensional physical quantities, a declared reference scale and an
absolute budget in matching units are required. Tolerances are fixture-specific
numerical policies, not performance thresholds. Where unchanged arithmetic is
the contract, exact entry comparisons are retained. Logical evaluation counts,
identity, ownership and invalidation are checked by indices/state, not inferred
from coordinate tolerances. Negative controls must demonstrate that an oracle
rejects an incorrect result; a repair's regression must fail without the repair.

Mass scopes use

$$
m_\rho(u,v)=\int_\Omega\rho\langle u,v\rangle\,dx,
\qquad
m_1(\Pi_h f,v)=\int_\Omega\langle f,v\rangle\,dx.
$$

The pairing conjugates the test value for complex fields and is Frobenius for
matrix fields. A complex coefficient does not imply a Hermitian operator.
Component-block checks retain structural zeros; arbitrary fields and component
couplings must not be replaced by an assumed nodal or interleaved DOF layout.
Basis-reuse contracts are checked by exact matrices and counted basis
evaluations, including mixed-degree and independent-copy states. Numerical
qualification alone does not establish a performance benefit.

Solver qualification records true residuals and iterations, not just successful
return codes. A performance repair must preserve applicable operators, loads,
constraints and solution behavior. Changes to numerical stopping policy must
not be presented as source-kernel speedups.

## Portable comparison protocol

Before/after comparisons hold the discrete workload, build configuration,
compiler options, dependency versions, backend, thread/rank counts and timed
scope fixed. Repetitions and reversed execution order distinguish candidate
effects from drift. Unchanged controls and actual work counts accompany the
comparison. An isolated kernel result requires confirmation in its consuming
operation before an end-to-end benefit is claimed.

Elapsed time, CPU time, allocation and memory observations answer different
questions. Sampled inclusive stack weights can overlap and are not additive
stage times. A frame filter is not a timed-region boundary: untimed setup or
oracle calls can appear on the same stacks. Instrumented profiles locate work;
uninstrumented runs establish elapsed comparisons. Platform-specific profilers
are optional evidence collectors, not library requirements or portable gates.

Results retain source revisions, case/filter, dimensions and counts, build and
runtime metadata, repetition policy, numerical qualification and profiler scope.
Machine-specific measurements and exploratory logs belong with their experiment
artifacts or review evidence, not in this architectural coverage document.
No portable speedup follows from one machine. Timing thresholds require a
controlled environment and characterized variability; these standalone
benchmarks impose no CTest timing gates.

Mass reference-storage profiling must distinguish scalar basis structure from
component expansion; its size depends on quadrature and local component/basis
counts. A compact representation requires a proven separability contract, not
numerical zero pruning. First-assembly attribution must distinguish setup,
allocation, kernel work and backend completion before assigning a cause to a
timing difference.

## Building and selecting measurements

Configure the ordinary build with `RODIN_BUILD_BENCHMARKS=ON`, then select an
opt-in target and inspect its registrations before choosing a scope:

```sh
cmake --build build --target RodinNamedMassKernelBenchmarks -j1
build/tests/benchmarks/RodinNamedMassKernelBenchmarks --benchmark_list_tests=true
build/tests/benchmarks/RodinNamedMassKernelBenchmarks \
  --benchmark_filter='^NamedMassKernel/H1/P2/Complex/Matrix/Hexahedron/Weighted/FreshKernel/2/3/real_time$' \
  --benchmark_min_time=0.1s --benchmark_repetitions=5 \
  --benchmark_out=mass-kernel.json --benchmark_out_format=json
```

The duration and repetition count above are invocation examples, not sufficient
qualification or a prescribed statistical confidence level. Executable/filter
selection affects initial cache population; aggregate and single-case processes
must not be compared as though they had identical histories. A numerical-only
smoke run establishes neither scaling nor repeatable performance. Only one
numerical workload is run at a time; a distributed workload may use multiple
ranks and an OpenMP workload may use multiple threads.

## Remaining work and completion criteria

The remaining work comprises the gaps in the ledger: complete mapped/embedded
and space/degree coverage, distributed stage attribution, constraints and
boundary/interface kernels,
nonlinear workflows, solver setup/factorization/preconditioning, and matched
multi-size end-to-end profiles. Follow-up priorities follow the assembly call
tree and measured source attribution. Benchmark availability, numerical
qualification, controlled timings, source attribution and consumer validation
are separate deliverables; no percentage is inferred from registration counts.
