# Shared settings of the CMBPE campaign. Sourced by the numbered scripts;
# every variable can be overridden from the environment.

CMBPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$CMBPE_DIR/../../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT_ROOT="${OUT_ROOT:-$PWD/cmbpe_runs}"

# Core budget shared by the concurrent jobs.
CORES="${CORES:-60}"

RMU_EXE="${RMU_EXE:-$BUILD/examples/Heart/CoronaryArtery/RmuAxi/RmuAxi}"
COR_EXE="${COR_EXE:-$BUILD/examples/Heart/CoronaryArtery/CFD/CoronaryArtery}"

# OpenMPI: --bind-to none lets several mpirun share the node without pinning
# their ranks to the same cores. MPICH/Intel MPI: MPIRUN="mpiexec -bind-to none".
MPIRUN="${MPIRUN:-mpirun --bind-to none}"

# One thread per rank: several jobs share the node.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

MESH_DIR="${MESH_DIR:-$CMBPE_DIR/meshes}"

STATES="${STATES:-healthy diabetic hyperviscous}"

# Carreau-Yasuda fits of table 1 (SI): mu0 muinf lambda n a.
cy_values() {
  case "$1" in
    healthy)      echo "0.0712616 0.0041 9.2780 0.3572 1.8773" ;;
    diabetic)     echo "0.0973528 0.0049622 15.7825 0.2670 1.0829" ;;
    hyperviscous) echo "0.3467269 0.0053264 17.4172 0.2200 0.6957" ;;
    *) echo "unknown state $1" >&2; return 1 ;;
  esac
}

cy_rmu() {
  local v; read -r -a v <<< "$(cy_values "$1")"
  echo "-rmu_cy_mu0 ${v[0]} -rmu_cy_muinf ${v[1]} -rmu_cy_lambda ${v[2]} -rmu_cy_n ${v[3]} -rmu_cy_a ${v[4]}"
}

cy_cor() {
  local v; read -r -a v <<< "$(cy_values "$1")"
  echo "-coronary_cy_mu0 ${v[0]} -coronary_cy_muinf ${v[1]} -coronary_cy_lambda ${v[2]} -coronary_cy_n ${v[3]} -coronary_cy_a ${v[4]}"
}

# Appends one job to a queue file: name np workdir after command.
# "after" is a comma-separated list of job names, or "-".
add_job() {
  local file="$1" name="$2" np="$3" dir="$4" after="$5"; shift 5
  printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$np" "$dir" "$after" "$*" >> "$file"
}
