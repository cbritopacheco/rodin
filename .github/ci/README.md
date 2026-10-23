# Dedicated test lane sources

[targeted-lanes.yml](targeted-lanes.yml) contains the full dedicated execution
matrices, numerical commands, environments, runners and timeout policies.
These definitions are source data, not a separately triggered workflow.

After editing this file or changing registration/grouping, run from the
repository root:

```sh
python3 dev/ci_schedule.py --write
python3 dev/ci_schedule.py
python3 -m unittest discover -s dev -p 'test_ci_*.py' -v
```

The generator writes eleven `TestConfiguration_*.yml` backend/instrumentation
variants and ten single-family `Test<Family>.yml` wrappers under
`.github/workflows/`. Each compatible configuration gets independent hierarchy
producers and parallel execution partitions. An execution group depends only
on its required producers. Configured CMake manifests must validate those
dependency bounds; registration references do not decide which targets exist.
Only applicable execution groups are declared in each variant. Sanitizer
variants omit publishers, and dedicated wrappers have no inactive family jobs.
This keeps unrelated branches out of the check list instead of generating
thousands of skipped checks.

`ConvergenceBaseline` obtains its directory matrix from the native planner;
the other dedicated families declare their full `include` matrices here.
Keep tests, MPI ranks, refinement sequences, environments and numerical
acceptance unchanged when reorganizing execution. Compilation recipes belong
in `dev/ci_builds.py`; ordinary registration belongs in the relevant test
directory's CMake file. The original numerical command/configuration snapshots
remain guarded even when scheduling changes.

See the [extension guide](../../dev/ci-testing.md) for specific changes and the
[workflow architecture](../workflows/README.md) for artifacts, coverage and
failure diagnosis.
