# Test job graphs

Open **Actions → Tests → a workflow run** for the native live dependency graph.
Click a job to inspect its logs. Matrix groups expand as their inputs become available.

The existing CMake / CI job summary links the live run and these commit-matched diagrams.
These generated diagrams show declared stage groups within each compatible configuration;
they do not show expanded matrix children, completion status or coverage percentages.

[Architecture and category overview](../../workflows/README.md#shared-builds-and-isolated-execution)

## Configuration diagrams

- [baseline](baseline.md)
- [complex](complex.md)
- [default_coverage](default_coverage.md)
- [default_sanitized](default_sanitized.md)
- [mpi](mpi.md)
- [native](native.md)
- [optional_coverage](optional_coverage.md)
- [optional_sanitized](optional_sanitized.md)
- [petsc](petsc.md)
- [ranges_complex](ranges_complex.md)
- [ranges_real](ranges_real.md)

Regenerate with `python3 dev/ci_schedule.py --write`; plain invocation checks consistency.
