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
- `Configuration.dox` explains build selections and links to the generated-header reference.
- `templates/` customizes the Modules index and public module headings in m.css.
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
- When a class or class template page describes a family with more than one
  supported specialization, include a table in the class description with the
  exact columns `Specialization` and `Description`. The specialization entry
  must link to that specialization, and the table must list the complete
  supported set for the family.
- In `@see` blocks and prose lists that name public classes, templates,
  specialization groups, or generated file pages, write explicit `@ref`
  references or HTML links. Do not rely on m.css to autolink bare comma-
  separated names such as `MeshBase, Mesh.h`.
- Prefer plain code text over `@ref` for names that are intentionally not part
  of the generated public reference.
- Keep generated files and local build output out of `doc/`.
- Resource payloads should live under `resources/` and be tracked with Git LFS.
  Documentation images and guide assets outside `resources/` stay in regular
  Git unless they become large enough to affect ordinary checkout cost.

## Agent Knowledge Base

`doc/agents/` is the entry point for AI coding agents and maintainers who need
architectural context. Start with `doc/agents/README.md`, then read the
philosophy and conventions pages before touching code. Domain-specific pages
cover geometry, variational forms, solvers, PETSc, physics, and testing.

## Module and Template Navigation

The Modules index groups public Rodin components using `RODIN_MODULES` and
`RODIN_MODULE_CATEGORIES` in `conf.py.in`. Each category starts with a description.
A table of contents links to each category and module. Each module is framed
in a single-column list with a linked name, an italic brief in dimmed text without a block,
one compact representative formula or usage example, and one introductory
description paragraph in a flat block. Complete descriptions remain on the linked module pages.
These descriptions come from the namespace documentation in `Namespaces.dox`; they are not maintained
in a separate overview. Public module namespace pages use the title “module”;
nested implementation namespaces retain their C++ namespace labels. The first
displayed formula or code block supplies the representative item; the opening
detailed paragraph supplies the description. Keep both concise and self-contained.

`rodin_add_mcss` invokes `dev/build_documentation.py` with the repository-pinned
m.css renderer and `doc/templates`. The helper renders the complete API and search
index once, then renders the concise Modules page from those same parsed module
descriptions. It scopes copied formula and section identifiers to keep SVG glyphs
and local anchors distinct, and preserves links to anchors on API pages. The
upstream templates remain the fallback for other pages.

Rebuilds remove generated XML and HTML pages while retaining formula caches and
assets, so deleted source pages disappear from navigation and search as well.
Local stylesheets use content-based filenames, so browsers refresh the palette
when CSS changes without discarding unchanged cached styles.

Before rendering, `dev/documentation_hierarchy.py` completes navigation
relationships in generated XML. Directory descriptions keep documented header
files reachable through the Files tree. Missing ownership links for nested types
are resolved against known enclosing classes and namespaces, ignoring scope
separators inside template arguments. Existing ownership links are preserved. The class index retains documented types
through undocumented enclosing scopes, keeping the omitted scope in their
displayed names. Private types excluded by the renderer are also excluded from
the index, preventing links to pages that are not generated.
This step changes generated documentation only; it does not change declarations
or the m.css submodule.

Template-specialization groups retain their Doxygen identifiers and are linked
from the generic class's `@see` documentation. They are technical references,
not entries in the Modules index. When adding a specialization family, include
that backlink on the primary template. Preprocessor macros stay attached to
their declaring header; the Configuration page links to `Configure.h.in` and
explains the generated `Rodin/Configure.h` include.

After rendering the full m.css site, run:

```sh
python3 dev/check_documentation_navigation.py
```

This checks every catalogue entry, representative formula/example, configuration
header link and visible primary-template specialization backlink. The
Documentation workflow runs the helper regression tests and the same check before
deployment. Adding a public
module requires an entry in `RODIN_MODULES` and a documented namespace with a
brief, a detailed description, and a representative formula or usage example.

## In-page Sidebar

The styled site provides an “On this page” sidebar on wide screens. The two
placement icons beside its heading choose the left or right gutter without
moving the centered content. The choice is saved in local storage; left is
the default.
`rodin-navigation.js` builds links from the page's section headings and marks the
current section while scrolling. Nested lists and connector lines show heading
hierarchy under Overview. Existing contents menus are excluded. The
sidebar scrolls independently and supports normal keyboard navigation.
Smaller screens and printed pages retain the original layout and inline contents.
The sidebar uses no external JavaScript dependencies; without JavaScript the
original page remains usable. Its layout is defined in `rodin.css`, and
`conf.py.in` loads the script for every generated page, and `Doxygen.mcss.in`
registers it as a copied site asset.

## Main-page Presentation

`Images/Architecture.svg` presents the core stack in blue, model modules on the
left in green, and third-party interfaces on the right in yellow. All blocks
have transparent backgrounds. `Mainpage.dox` embeds this vector diagram.

`dev/documentation_highlighting.py` supplements Pygments with documented C++ type
names and enum declarations from Doxygen XML. The configured post-filter colors
types blue and enum values red, preserving variables, strings and comments.
Qualified enum values are matched against their declaring scopes to distinguish
them from identically named classes. The shared palette in `rodin.css` applies
to code throughout the site and uses normal font weight for every token.

## Color Themes

The site uses a dark base palette and follows the system's light preference
until a theme is selected explicitly. The icon at the right of the top bar
opens a System / Light / Dark menu below the icon, inside the page. Arrow keys,
Home/End, Enter/Space, and Escape provide keyboard navigation; clicking outside
closes the menu. Its choice is saved
in local storage; selecting System resumes automatic changes. The two existing
transparent logos switch with the active theme.

`rodin-theme.css` owns both palettes and the section-heading colors. Page titles
and ordinary subsection titles stay neutral. Main sections use blue; module
category and module headings use blue for Core, green for Models, gold for
integrations, and purple for Infrastructure. Status blocks retain semantic
colors. Light-mode accents are darker to preserve readability.

`rodin-theme.js` applies the saved choice before rendering and adds the selector.
`RODIN_MCSS_STYLESHEETS` in `conf.py.in` lists the pinned upstream CSS components.
The documentation builder bundles them with their variables intact and loads
Rodin's overrides afterwards. Upstream m.css is not modified. Update palettes
in Rodin's stylesheet rather than changing the vendor files or generated HTML.
