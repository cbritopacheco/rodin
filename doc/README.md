# Documentation Workspace

The `doc/` directory is Rodin's documentation workspace. It contains the
Doxygen configuration, narrative guide pages, examples index pages, generated
site styling, bibliography, image assets, and the agent knowledge base used by
AI coding assistants.

This directory is for source documentation, not generated output. Build
artifacts are written under the CMake build tree.

## Directory Map

- `Doxygen.in` configures the standard Doxygen build.
- `Doxygen.mcss.in` configures the XML-producing Doxygen pass used by m.css.
- `conf.py.in` configures m.css presentation, navigation, and search.
- `Mainpage.dox`, `Namespaces.dox`, `Notation.dox`, and `Credits.dox` provide
  top-level reference pages.
- `Guides/` contains narrative user guides and conceptual documentation.
- `Examples/` and `ExamplesIndex.dox` organize example-program pages.
- `Images/`, `favicon.ico`, and `rodin.css` provide visual assets and styling.
- `References.bib` stores bibliography entries used by Doxygen citations.
- `agents/` contains the AI-agent knowledge base. It is documentation about
  how Rodin is designed and maintained, not user-facing API reference.
- `README/` stores assets used by the repository-level README.

## Building the Documentation

Rodin builds documentation through the `RodinDoxygen` CMake target.

From a configured build tree:

```sh
cmake --build build --target RodinDoxygen
```

For a fresh build tree:

```sh
cmake -S . -B build -DRODIN_BUILD_DOC=ON
cmake --build build --target RodinDoxygen
```

The generated documentation is placed under `build/doc/`.

## Build Modes

The documentation build has two modes:

| Option | Purpose |
|--------|---------|
| `RODIN_BUILD_DOC=ON` | Enable the `RodinDoxygen` target. |
| `RODIN_USE_MCSS=ON` | Build the styled m.css site instead of plain Doxygen HTML. |

Plain Doxygen:

```sh
cmake -S . -B build -DRODIN_BUILD_DOC=ON -DRODIN_USE_MCSS=OFF
cmake --build build --target RodinDoxygen
```

Styled m.css documentation:

```sh
cmake -S . -B build -DRODIN_BUILD_DOC=ON -DRODIN_USE_MCSS=ON
cmake --build build --target RodinDoxygen
```

Requirements for plain Doxygen:

- Doxygen 1.8.15+

Additional requirements for m.css:

- Python 3.6+
- Jinja2
- Pygments
- LaTeX with the packages required by m.css

The Python packages for the m.css path are listed in `requirements.txt`.

## Warning Policy

Doxygen warnings are tracked by the ratchet in `dev/check_doxygen_warnings.py`.
The accepted baseline lives in `dev/doxygen_warnings.baseline`.

Run the check with:

```sh
python3 dev/check_doxygen_warnings.py
```

If a change fixes old warnings, shrink the baseline with:

```sh
python3 dev/check_doxygen_warnings.py --update-baseline
```

Do not grow the baseline to hide new warnings. Fix the documentation instead.

The m.css renderer can emit Python logging diagnostics such as missing
programlisting metadata or section nesting limits after Doxygen has generated
XML. Those diagnostics should be cleaned up over time, but they are separate
from the Doxygen warning ratchet unless the build target fails.

## Writing Documentation

Prefer documenting the mathematical concept represented by the code. API
comments should state the formula, space, domain, and ownership/lifetime
expectations when those matter. For tests, document the behavior being checked:
read the test body and describe the invariant, regression, or numerical
property the test protects.

Useful conventions:

- Use Doxygen commands such as `@brief`, `@param`, `@tparam`, `@return`,
  `@ingroup`, and `@defgroup` where they map cleanly to the code.
- Use `@f$ ... @f$` for inline mathematics and `@f[ ... @f]` for displayed
  formulas.
- Use links only when the target is documented and resolvable.
- Prefer plain code text over `@ref` for names that are intentionally not part
  of the generated public reference.
- Keep generated files and local build output out of `doc/`.

## Agent Knowledge Base

`doc/agents/` is the entry point for AI coding agents and maintainers who need
architectural context. Start with `doc/agents/README.md`, then read the
philosophy and conventions pages before touching code. Domain-specific pages
cover geometry, variational forms, solvers, PETSc, physics, and testing.
