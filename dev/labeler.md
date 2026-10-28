# Pull-request labels

[The label rules](../.github/labeler.yml) classify changed paths. Every pattern
under `any-glob-to-any-file` is an alternative: a single matching file is
sufficient. Do not replace those alternatives with an `all` condition.

[Labeler.yml](../.github/workflows/Labeler.yml) uses a pinned Node.js 24 release
of `actions/labeler`. Dot-directory matching is enabled, so `.github` changes
are included. Label synchronization is disabled: existing manual labels,
including labels also managed by these rules, are preserved. The action can
create newly introduced labels using its scoped issue-write permission.

The write-enabled workflow uses `pull_request_target` and reads rules from the
trusted base commit. It never checks out or executes pull-request code. Changes
to its action version and rules take effect after merging into the target
branch. The separate read-only
[Labeler Rules workflow](../.github/workflows/LabelerRules.yml) tests proposed
changes before that merge.

## Label catalog and colours

[`.github/labels.yml`](../.github/labels.yml) records the name, colour and
description of every repository label, including labels used manually for issue
triage. The labeler rules use `convergence`, `installation` and `nafems` for those
three test categories; unit and manufactured tests retain `unit-tests` and
`manufactured-tests`.

Every label has its own shade; related labels remain in the same colour family.

| Family | Shade range | Labels |
| :----- | :---------- | :----- |
| Core and numerical modules | Dark blues | Core, geometry, algebra, forms, assembly and solvers |
| Models and adaptation | Dark greens | Adaptation, transport, distance, cardiac and solid models |
| External backends | Dark ambers | MFEM, MMG, MPI, PETSc and Scotch |
| Documentation and data interfaces | Light slate blues | Documentation, examples, IO and serialization |
| Tests and validation | Dark teals | Unit, manufactured, convergence, installation, NAFEMS and benchmarks |
| Development infrastructure | Dark purples | CI, CMake, dependencies, tooling, Python, resources and Codex |
| Defects | Dark reds | Bugs and hotfixes |
| Contributions | Light greens | Good first issues and help wanted |
| Disposition | Light greys | Duplicate, invalid and wontfix |

Enhancements use a dark green, questions a light lavender, and work in progress
a light amber. Names and descriptions carry the meaning; colour is a secondary
grouping cue.

GitHub controls the foreground and theme-specific rendering; its label API
accepts a background colour, not a separate text colour. The palette avoids
mid-tones: dark fills are checked against light text and light fills against
dark text, with at least 4.5:1 contrast. Regression tests require distinct
colours, related hues and readable contrast. Labels may be rendered differently
by GitHub's light and dark themes.

When adding a rule, add its metadata to the catalog and create or update the
GitHub label with those values. Rename existing labels in place to preserve
issue and PR assignments. The catalog is checked by the Labeler Rules workflow;
it does not automatically write repository label metadata. The action can still
create an absent label, but its default metadata should then be reconciled with
the catalog. Neither reconciliation nor the labeler should delete manual labels.

## Mapping policy

- Module labels follow every directory under `src/Rodin/` and its umbrella
  header. Existing label names such as `alerts`, `form-language` and `qf` remain.
- The same module labels apply to unit, manufactured and benchmark paths under
  `tests/**/Rodin/`, including nested backend and `Models` directories. For
  example, an MPI Geometry test receives `unit-tests`, `mpi` and `geometry`;
  a PETSc IO test receives `unit-tests`, `petsc` and `io`.
- Core root headers and root-level utility tests receive `core`; Threads tests
  belong to that same core vocabulary.
- Finite element test directories such as P0, P1 and H1, and convergence studies,
  receive `variational`. Distributed convergence paths also receive their MPI
  and PETSc labels. Geometry approximation studies additionally receive
  `geometry`. Linear elasticity paths retain `linear-elasticity` and identify
  their owning `solid` module.
- Category labels distinguish unit, manufactured, convergence, installation,
  NAFEMS and benchmark tests. Documentation includes nested Markdown and
  citation files. CI, developer tooling, Python bindings, resources and
  third-party dependencies have their own path mappings.

Labels describe these declared ownership and category mappings, not every
library a test links against. A test may have several relevant labels; there is
no exclusive first-match rule.

## Adding or changing a mapping

Edit `.github/labeler.yml`, keeping alternative paths together under
`any-glob-to-any-file`. Add a module's production directory, umbrella header and
test hierarchy together. For a test directory named after a finite element or
physics concept instead of a module, add its ownership alias explicitly.

Add an independent single-file regression in [test_labeler.py](test_labeler.py).
The module inventory assertion requires a mapping for every source module.
Include negative cases for similarly named unrelated directories. When changing
the action release, update both workflow pins and the regression's expected
package version together.

## Running the regressions

Use Node.js 24 or later and PyYAML 6.0.2. Check out the exact action revision
referenced in the workflows into a temporary directory:

```sh
git clone https://github.com/actions/labeler.git /tmp/rodin-labeler-action
git -C /tmp/rodin-labeler-action checkout --detach b8dd2d9be0f68b860e7dae5dae7d772984eacd6d
python3 -m pip install PyYAML==6.0.2
LABELER_ACTION_PATH=/tmp/rodin-labeler-action \
  python3 -m unittest discover -s dev -p test_labeler.py -v
```

`NODE_BINARY` can select an explicit Node.js executable. Missing prerequisites
fail the suite rather than silently skipping it.

The tests run the pinned action's bundled entry point against a local HTTP mock
of the GitHub API. They exercise the actual configuration parser, minimatch
semantics, base-commit configuration retrieval, paginated changed-file lists,
label updates and preservation of manual labels. No real token or PR mutation
is involved. Test cases execute serially.
