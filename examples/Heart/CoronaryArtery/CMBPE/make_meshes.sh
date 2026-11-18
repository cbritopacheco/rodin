#!/usr/bin/env bash
# Regenerates the coarse and fine coronary meshes from the medium one with
# mmg3d (third-party/mmg: cmake -S third-party/mmg -B mmg-build && cmake
# --build mmg-build --target mmg3d). Only Vertices, Triangles and Tetrahedra
# are kept, as in the medium mesh.
#   MMG3D=/path/to/mmg3d_O3 ./make_meshes.sh
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/meshes"
MMG3D="${MMG3D:-mmg3d_O3}"

clean() {   # in out
  python3 - "$1" "$2" <<'EOF'
import sys
keep = {'Vertices': 4, 'Triangles': 4, 'Tetrahedra': 5}
skip = {'Edges': 3, 'Normals': 3, 'Tangents': 3, 'NormalAtVertices': 2, 'TangentAtVertices': 2,
        'Corners': 1, 'RequiredVertices': 1, 'Ridges': 1, 'RequiredEdges': 1, 'RequiredTriangles': 1}
f = open(sys.argv[1]).read().split()
out = ["MeshVersionFormatted 2", "", "Dimension 3", ""]
i = 0
while i < len(f):
    w = f[i]
    if w in keep:
        n, per = int(f[i + 1]), keep[w]
        out += [w, str(n)] + [" ".join(f[i + 2 + k * per:i + 2 + (k + 1) * per]) for k in range(n)] + [""]
        i += 2 + n * per
    elif w in skip:
        i += 2 + skip[w] * int(f[i + 1])
    else:
        i += 1
out.append("End")
open(sys.argv[2], "w").write("\n".join(out) + "\n")
EOF
}

tmp="$(mktemp -d)"
"$MMG3D" -in "$DIR/coronary_medium.mesh" -out "$tmp/coarse.mesh" -hmax 1.0 -hausd 0.3 -v 0
"$MMG3D" -in "$DIR/coronary_medium.mesh" -out "$tmp/fine.mesh" -hmax 0.45 -hausd 0.1 -v 0
clean "$tmp/coarse.mesh" "$DIR/coronary_coarse.mesh"
clean "$tmp/fine.mesh" "$DIR/coronary_fine.mesh"
rm -rf "$tmp"
