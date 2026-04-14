#!/usr/bin/env bash
set -euo pipefail

root_dir="${1:-/home/wjj/Project/stage/levelset_stoke/build/out}"
surface_ref="${2:-13}"

shopt -s nullglob
all_meshes=("${root_dir}"/Omega.*.mesh)
old_outputs=("${root_dir}"/Omega.*.obstacle.mesh)
inputs=()

for mesh in "${all_meshes[@]}"; do
  case "${mesh}" in
    *.obstacle.mesh) ;;
    *) inputs+=("${mesh}") ;;
  esac
done

if [ ${#inputs[@]} -eq 0 ]; then
  echo "No Omega.*.mesh files found under ${root_dir}" >&2
  exit 1
fi

if [ ${#old_outputs[@]} -gt 0 ]; then
  echo "Removing existing obstacle surfaces under ${root_dir}"
  rm -f "${old_outputs[@]}"
fi

for mesh in "${inputs[@]}"; do
  out="${mesh%.mesh}.obstacle.mesh"
  echo "Extracting ref ${surface_ref}: ${mesh} -> ${out}"

  awk -v target_ref="${surface_ref}" '
    function trim(s) {
      sub(/^[[:space:]]+/, "", s)
      sub(/[[:space:]]+$/, "", s)
      return s
    }

    BEGIN {
      state = ""
      num_vertices = 0
      num_triangles = 0
      vertex_idx = 0
      triangle_idx = 0
      used_count = 0
      kept_triangles = 0
    }

    /^ *Vertices$/ {
      state = "vertices_count"
      next
    }

    /^ *Triangles$/ {
      state = "triangles_count"
      next
    }

    /^ *Edges$/ {
      state = "edges_count"
      next
    }

    /^ *Tetrahedra$/ {
      state = "stop"
      next
    }

    /^ *End$/ {
      state = "stop"
      next
    }

    state == "vertices_count" {
      num_vertices = int($1)
      state = "vertices"
      next
    }

    state == "vertices" {
      vertex_idx++
      vertices[vertex_idx] = $0
      if (vertex_idx >= num_vertices)
        state = ""
      next
    }

    state == "triangles_count" {
      num_triangles = int($1)
      state = "triangles"
      next
    }

    state == "edges_count" {
      num_edges = int($1)
      edges_read = 0
      state = "edges"
      next
    }

    state == "edges" {
      edges_read++
      if (edges_read >= num_edges)
        state = ""
      next
    }

    state == "triangles" {
      triangle_idx++
      if ($4 == target_ref) {
        tri1[triangle_idx] = $1
        tri2[triangle_idx] = $2
        tri3[triangle_idx] = $3
        tri_ref[triangle_idx] = $4
        if (!($1 in used)) {
          used[$1] = 1
          used_order[++used_count] = $1
        }
        if (!($2 in used)) {
          used[$2] = 1
          used_order[++used_count] = $2
        }
        if (!($3 in used)) {
          used[$3] = 1
          used_order[++used_count] = $3
        }
        keep[++kept_triangles] = triangle_idx
      }
      if (triangle_idx >= num_triangles)
        state = ""
      next
    }

    END {
      print "MeshVersionFormatted"
      print "2"
      print ""
      print "Dimension"
      print "3"
      print ""
      print "Vertices"
      print used_count

      for (i = 1; i <= used_count; ++i) {
        old_id = used_order[i]
        remap[old_id] = i
        print vertices[old_id]
      }

      print ""
      print "Triangles"
      print kept_triangles

      for (i = 1; i <= kept_triangles; ++i) {
        idx = keep[i]
        print remap[tri1[idx]], remap[tri2[idx]], remap[tri3[idx]], tri_ref[idx]
      }

      print ""
      print "End"
    }
  ' "${mesh}" > "${out}"
done

echo "Done."
