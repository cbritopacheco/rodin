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
- `../brand/` is the canonical source for the logo, wordmark, SVG favicon and bundled Manrope fonts.
- `Images/`, `rodin.css` and `rodin-theme.css` provide illustrations and site styling.
  Brand assets and font licenses are copied into the generated site; no external font service is required.
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

- Doxygen 1.15+ (CI pins 1.15.0)

Additional requirements for m.css:

- Python 3.9+
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

The documentation build also exports `rodin-module-catalogue.json` with the same
rendered module overviews, category descriptions, ownership paths and stylesheet
order. The benchmark navigation consumes this catalogue rather than duplicating
namespace prose. Namespace-backed entries derive their `Rodin/…` ownership path
from the C++ namespace. Documentation groups specify that path in
`RODIN_MODULE_PATHS` in `conf.py.in`; their Doxygen group identifiers are not
source directories or physical contexts. See [benchmark registration and display](../tests/benchmarks/README.md#registration-and-display).

## Benchmark documentation and publication

The m.css build also renders `benchmarks/index.html` and its module/class and
physical-context pages. It reads the same explicit CMake registrations as the
benchmark build, without compiling or executing numerical workloads. Authored
HTML fragments live beside benchmark drivers; canonical module descriptions
come from the rendered API module catalogue. The common theme and in-page
sidebar remain shared with the API documentation. Benchmark pages use the same
content grid and heading sizes; sidebar placement never changes the content
column. Top-bar branding aligns with the content headings.

On main-branch pushes the Benchmarks workflow builds documentation and a complete
benchmark campaign, then publishes both in one Pages commit under
`docs/refs/heads/<branch>/`. Component histories are stored separately in
`benchmark-data/refs/heads/<branch>/` and connected to the generated leaf pages.
A failed or incomplete campaign cannot replace the public documentation or
publish partial timing histories. Pull-request documentation builds upload an
artifact for validation without publishing it.

See [benchmark registration, shared builds and publication](../tests/benchmarks/README.md)
for the full graph, naming rules, timing contracts and extension instructions.

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

## Footer

The shared footer is configured through `FINE_PRINT` in `conf.py.in`. It keeps
the copyright notice separate from the Doxygen/m.css credits, including the
detected Doxygen version. `HTML_EXTRA_FILES` in `Doxygen.in` bundles the Doxygen
wordmark; a CSS mask gives it the surrounding text colour. The m.css credit is a
bold text link. Artwork provenance is in
[`Images/README.md`](Images/README.md). The footer stacks its two parts on narrow
screens, without changing upstream m.css templates.

## Class and Specialization References

`dev/documentation_reference.py` presents reference pages after the pinned m.css
renderer has produced the API and search index. Class titles use the local
class name; namespace breadcrumbs and a linked **Defined in** line establish
the scope and source header. The complete declaration appears separately in a
highlighted code block. Template parameters and existing member documentation
remain available below it.

Specializations list their argument bindings against the primary template's
parameter names. Parameters that remain generic are preserved in the declaration;
concrete instantiations are never inferred. Primary-template pages and class
lists group the available specializations into a two-column table. Every row
links to its original page and expands to show the complete declaration.
In class lists, this table is initially collapsed under “See specializations”,
with a plus/minus indicator. The primary declaration and brief remain visible;
keyboard activation expands the full list without JavaScript.
The Doxygen configuration enables `RODIN_USE_OPENMP` for extraction so the
parallel assembly primary declarations are documented even in serial builds.
If a primary declaration is unavailable, the renderer groups specializations by
fully qualified family name without inventing template parameters or a primary
page. Families in different namespaces remain separate.

Overview declarations longer than 120 characters (with whitespace normalized)
show only the linked symbol and brief, with the original linked declaration under
“See declaration”. This is a presentation threshold, independent of viewport
width; dedicated reference pages retain their complete declarations.
Source-header annotations consistently use “Defined in” with a linked header,
below the declaration. Example code retains its literal `#include` directives.
Deduction-guide details use the linked class name as their heading and keep the
full original declaration behind “See declaration”. Source-header links, the
description, and parameter documentation remain visible. Guides are identified
from their Doxygen member metadata, not from descriptive prose or name length.

Member signatures omit the current class or namespace qualification and avoid
repeating the enclosing class template declaration. Return types, member template
parameters, qualifiers, external namespace qualifications and original member
anchors are retained. Search continues to use the complete Doxygen names.
Included API headers without authored file prose receive a generated file brief
so **Defined in** links resolve to actual file pages. C++ sources and upstream
m.css remain unchanged.

### Search Ranking

Class-name searches promote matching general class and primary-template pages
before specializations and members, across every module. Among equally close
matches, shorter qualified names come first. The requested result limit is applied
after promotion, so a primary excluded by the initial broad lookup can still be
returned. Deprecated and deleted declarations retain their normal placement.
Queries containing template arguments, or ending in `::` or `/`, retain the
upstream behavior for specialized types and scope-member lists.

`dev/documentation_search.py` derives the general-class catalogue from Doxygen
XML. Family landing pages backed by otherwise undocumented forward declarations
receive a generated brief so that they are indexed as well. `rodin-search.js`
uses the existing search index to retrieve and rank results; complete names,
links, autocomplete, and keyboard navigation remain available. The builder
appends this extension to the unchanged upstream search engine and publishes a
content-versioned bundle. Vendor files are not modified. Search-ranking regression tests use Node.js in
addition to the Python documentation test suite.

### Inheritance and Type Aliases

Reference pages show **Inherits from** and **Derived classes** tables using the
actual Doxygen inheritance edges. Linked briefs supply the responsibility or
description; full type arguments, access and virtual-inheritance qualifiers remain
available. Template specialization alone does not create an inheritance edge.

A small **Direct inheritance** diagram shows the current class and at most three
distinct base and child class families on either side, grouping repeated
specializations. Whenever more than three underlying base or derived entries
exist, a **+N more** link leads to the corresponding complete table, including
when specialization grouping leaves fewer than three diagram nodes. Derived-class
tables initially show at most three entries, prioritizing distinct families.
A native **See more / See less** disclosure button reveals the remaining entries;
all exact types and links remain in the generated page. The diagram uses ordinary linked HTML and CSS, with no Graphviz or
JavaScript dependency. It wraps on small screens.

**Overrides** links the overriding members to the base declarations identified by
Doxygen. **Inherited interface** groups public and protected members by their
actual declaring base, inside expandable lists. Constructors, destructors and
private members and friends are excluded from those inherited lists. These sections preserve
the original member documentation and anchors rather than duplicate contracts.
Improve the shared contract in a base class's own Doxygen comment; explain
special behavior in its overrides.

Type-alias detail headings omit the current enclosing scope before the alias
name, while retaining the aliased type, external qualifiers, literals and links.
These transformations run inside `RodinDoxygen` and the documentation CI build.

Standalone **See also** class links already represented by inheritance or
specialization navigation are omitted. Other cross-references, explanatory prose,
technical specialization references and member-specific contracts are retained.

### Complete reference links

Doxygen 1.15.0 is pinned in documentation and warning CI. Version 1.14 can omit
explicit template specializations when they have aliases (upstream issue
[doxygen/doxygen#11761](https://github.com/doxygen/doxygen/issues/11761)).
Use the same parser for rendering and the zero-warning ratchet.

After rendering, references to undocumented forward declarations or private
enclosing scopes retain their text but have no hyperlink if no page exists.
Documented classes and stale hand-written links are never silently unlinked.
The navigation check scans every HTML page, including benchmark pages, and
rejects missing relative page, image, script and stylesheet destinations.
Prefer `@ref` to manually spelling generated specialization filenames.

Before m.css runs, `prepare_member_definitions()` resolves Doxygen's referenced
namespace/file members from their authoritative XML definitions. This preserves
their descriptions, source locations and anchor IDs instead of dropping a whole
section when it contains a reference. The upstream m.css checkout is unchanged.
Enum values without detailed sections receive their XML anchor IDs at their
existing declaration entries. The complete-site check also rejects missing
fragment destinations; the navigation control `#!` is intentionally exempt.

Doxygen `@warning` notes use the native m.css `m-note m-warning` component and
the theme's warning palette. Other danger notes retain `m-danger`. First note
headings have no inherited section-heading top margin. Unboxed declaration
qualifiers use their semantic text colours rather than filled-badge colours.
