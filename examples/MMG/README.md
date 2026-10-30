# MMG examples

Build Rodin with MMG support and run each executable in a separate scratch
directory. Remeshing changes topology and numbering; construct new finite
element spaces and problems before using the output mesh in another solve.

| Target | Operation | Checks |
| :--- | :--- | :--- |
| `MeshIO` | Write and reload a generated unit-square mesh | Mesh counts, corner set and required vertices |
| `MeshOptimization` | Optimize a perturbed triangle mesh | Positive triangle quality and unit area |
| `MeshAdaptation` | Remesh from the prescribed scalar field `0.05 + 0.15*x` | Positive triangle areas and unit area |

The last two programs export `Before.csv` and `After.csv` alongside MEDIT
meshes. Their figures can be regenerated with:

```sh
python examples/MMG/visualize_meshes.py OUTPUT_DIRECTORY FIGURE.png
```

The plots use the actual exported triangle coordinates, a common domain scale,
Manrope labels and transparent backgrounds. Shape quality is reported as
$4\sqrt{3}A/(\ell_1^2+\ell_2^2+\ell_3^2)$ in the optimization example.
