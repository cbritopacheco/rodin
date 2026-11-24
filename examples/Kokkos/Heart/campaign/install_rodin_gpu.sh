#!/usr/bin/env bash
# =============================================================================
#  install_rodin_gpu.sh
#
#  Rodin (rama model/AtriumKokkos) en una máquina Linux con GPU NVIDIA:
#  CUDA Toolkit, Kokkos y Kokkos-Kernels 5.2.2, cuDSS, PETSc 3.26.0 con CUDA y
#  Kokkos (MUMPS, ScaLAPACK, METIS, ParMETIS, PT-Scotch, hypre), Rodin, la
#  malla Coarse.mesh y una prueba corta en la GPU. La arquitectura se lee de la
#  tarjeta (RTX 5090: compute capability 12.0, Kokkos_ARCH_BLACKWELL120, sm_120).
#
#     bash install_rodin_gpu.sh deps    # con sudo: paquetes y driver (interactivo)
#     bash install_rodin_gpu.sh check   # solo diagnostica, no instala nada
#     bash install_rodin_gpu.sh build   # sin sudo: todo lo demás (con nohup)
#
#  Todo queda en $PREFIX (~/gpu). No toca ~/.bashrc ni otras instalaciones.
#  Se puede relanzar después de un fallo: salta lo que ya está hecho.
#
#  Variables opcionales: PREFIX, JOBS, BRANCH, TARGETS ("all" compila todo),
#  CXX_HOST (p. ej. g++-13), ARCH (p. ej. BLACKWELL120), CUDA_VERSION
#  (p. ej. 12.8.1), NO_HYPRE=1, SKIP_TEST=1.
# =============================================================================
set -Eeuo pipefail
MODE=${1:-build}
PREFIX=${PREFIX:-$HOME/gpu}
BRANCH=${BRANCH:-model/AtriumKokkos}
REPO=${REPO:-https://github.com/cbritopacheco/rodin}
TARGETS=${TARGETS:-"KokkosAtriumFlow KokkosAtriumFlowPETSc KokkosAtriumRheologyPETSc"}
KOKKOS_VERSION=5.2.2
PETSC_VERSION=3.26.0
MESH_URL=https://media.githubusercontent.com/media/cbritopacheco/rodin/model/Atrium/resources/examples/Heart/Atrium/Coarse.mesh
MESH_SHA=5a62edfd31040df475008a96722335e6a3e821096a6c84b56d9711c147ffb4d3
LOG=$PREFIX/logs
mkdir -p "$PREFIX" "$LOG" "$PREFIX/src" "$PREFIX/tmp"
cd "$PREFIX"

# Nada de conda ni de otros Kokkos/PETSc en el entorno: CMake los tomaría.
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -v -i -E 'conda|mamba|miniforge' | paste -sd: -)
unset Kokkos_ROOT Kokkos_DIR KokkosKernels_ROOT KokkosKernels_DIR CMAKE_PREFIX_PATH \
      PETSC_DIR PETSC_ARCH cudss_DIR CONDA_PREFIX 2>/dev/null || true
export GIT_LFS_SKIP_SMUDGE=1 TMPDIR=$PREFIX/tmp

STEP=inicio
say()  { STEP="$*"; printf '\n\033[1m[%s] ==> %s\033[0m\n' "$(date +%H:%M)" "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[31mERROR (%s):\033[0m %s\n' "$STEP" "$*" >&2; exit 1; }
trap 'printf "\n\033[31mFalló en: %s  (logs en %s)\033[0m\n" "$STEP" "$LOG" >&2' ERR
try()  { local log=$LOG/$1; shift; "$@" > "$log" 2>&1; }
run()  { local log=$LOG/$1; shift; "$@" > "$log" 2>&1 || { tail -30 "$log" >&2; die "ver $log"; }; }
ver_ge() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -1)" = "$2" ]; }
min() { if [ "$1" -lt "$2" ]; then echo "$1"; else echo "$2"; fi; }

# =============================================================================
#  deps: paquetes del sistema y driver (necesita sudo)
# =============================================================================
if [ "$MODE" = deps ]; then
  command -v apt-get > /dev/null || die "solo sé instalar con apt (Ubuntu/Debian); instala a mano: compiladores C/C++/Fortran, cmake, ninja, git, curl, Boost, Eigen, OpenMPI, HDF5, SuiteSparse y el driver NVIDIA."
  say "Paquetes del sistema (pide la clave de sudo una vez)"
  sudo -v
  sudo apt-get update
  PKGS=(build-essential gfortran cmake ninja-build git curl wget ca-certificates pkg-config
        python3 xz-utils bzip2 patch file flex bison zlib1g-dev libnuma-dev pciutils
        libboost-all-dev libeigen3-dev libopenmpi-dev openmpi-bin libhdf5-dev libsuitesparse-dev)
  GV=$( (gcc -dumpversion 2>/dev/null || echo 0) | cut -d. -f1)
  # Rodin necesita g++ >= 12; CUDA 12.x acepta hasta g++ 14.
  if [ "$GV" -lt 12 ]; then PKGS+=(gcc-12 g++-12 gfortran-12); fi
  if [ "$GV" -gt 14 ]; then PKGS+=(gcc-14 g++-14 gfortran-14); fi
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y "${PKGS[@]}"

  say "Driver NVIDIA"
  lspci | grep -i -E 'vga|3d' | grep -i nvidia | sed 's/^/    /' || die "no veo una tarjeta NVIDIA en lspci."
  DRV=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1 || true)
  # Las RTX 50 (Blackwell) necesitan el driver 570 o más nuevo, con los módulos "open".
  if [ -n "$DRV" ] && ver_ge "$DRV" 570; then
    info "driver $DRV: sirve"
  else
    if [ -n "$DRV" ]; then info "driver $DRV: muy antiguo para una RTX 50"; else info "no hay driver cargado"; fi
    PKG=$(apt-cache search --names-only '^nvidia-driver-[0-9]+-open$' | awk '{print $1}' | sort -V | tail -1)
    if [ -n "$PKG" ] && ver_ge "$(echo "$PKG" | tr -dc 0-9)" 570; then
      info "instalando $PKG"
      sudo DEBIAN_FRONTEND=noninteractive apt-get install -y "$PKG"
    elif command -v ubuntu-drivers > /dev/null; then
      sudo ubuntu-drivers install
    else
      die "no encontré un driver >= 570 en apt: instálalo desde nvidia.com y reinicia."
    fi
    say "Driver instalado: REINICIA el PC y vuelve a pegar el mismo comando"
    info "(si el PC tiene Secure Boot, al reiniciar te pedirá registrar la clave MOK que eligió apt)"
    exit 3
  fi
  say "Dependencias del sistema listas"
  exit 0
fi

# =============================================================================
#  1. GPU, CUDA y compilador
# =============================================================================
say "GPU"
nvidia-smi > /dev/null 2>&1 || die "nvidia-smi no responde: corre primero 'bash $0 deps' y reinicia."
GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d ' ')
DRV=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1)
DRV_CUDA=$(nvidia-smi | sed -n 's/.*CUDA Version: *\([0-9.]*\).*/\1/p' | head -1)
MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -1)
if [ -z "${ARCH:-}" ]; then
  case "$CC" in
    7.0) ARCH=VOLTA70 ;;   7.5) ARCH=TURING75 ;;   8.0) ARCH=AMPERE80 ;;
    8.6) ARCH=AMPERE86 ;;  8.7) ARCH=AMPERE87 ;;   8.9) ARCH=ADA89 ;;
    9.0) ARCH=HOPPER90 ;;  10.0) ARCH=BLACKWELL100 ;; 12.0) ARCH=BLACKWELL120 ;;
    *) die "compute capability $CC sin mapeo en Kokkos $KOKKOS_VERSION: fija ARCH=..." ;;
  esac
fi
SM=$(echo "$CC" | tr -d .)
MIN_CUDA=12.2; if ver_ge "$CC" 10.0; then MIN_CUDA=12.8; fi      # Blackwell: CUDA >= 12.8
info "$GPU ($MEM), compute capability $CC -> Kokkos_ARCH_$ARCH, sm_$SM"
info "driver $DRV (CUDA hasta $DRV_CUDA)"
ver_ge "$DRV_CUDA" "$MIN_CUDA" || die "el driver solo llega a CUDA $DRV_CUDA y esta GPU necesita $MIN_CUDA: corre 'bash $0 deps'."

say "CUDA Toolkit"
cuda_ok() {   # <nvcc>: versión >= MIN_CUDA y <= la del driver
  local v; v=$("$1" --version 2> /dev/null | sed -n 's/.*release \([0-9]*\.[0-9]*\).*/\1/p')
  # Solo CUDA 12: es con la que se prueban Kokkos 5.2.2 y PETSc 3.26.
  [ -n "$v" ] && [ "${v%%.*}" = 12 ] && ver_ge "$v" "$MIN_CUDA" && ver_ge "$DRV_CUDA" "$v"
}
NVCC=""
if [ -z "${CUDA_VERSION:-}" ]; then
  for c in "$(command -v nvcc || true)" "$PREFIX"/cuda-*/bin/nvcc /usr/local/cuda/bin/nvcc /usr/local/cuda-*/bin/nvcc; do
    if [ -x "$c" ] && cuda_ok "$c"; then NVCC=$(readlink -f "$c"); break; fi
  done
fi
if [ -z "$NVCC" ]; then
  # Runfiles: versión y driver con que vienen. Se instala solo el toolkit, sin root.
  declare -A RUNFILE=([12.9.1]=575.57.08 [12.8.1]=570.124.06 [12.6.3]=560.35.05 [12.4.1]=550.54.15)
  if [ -z "${CUDA_VERSION:-}" ]; then
    for v in 12.9.1 12.8.1 12.6.3 12.4.1; do
      if ver_ge "$DRV_CUDA" "${v%.*}" && ver_ge "${v%.*}" "$MIN_CUDA"; then CUDA_VERSION=$v; break; fi
    done
  fi
  [ -n "${CUDA_VERSION:-}" ] && [ -n "${RUNFILE[$CUDA_VERSION]:-}" ] || die "no hay un CUDA que calce con el driver (CUDA $DRV_CUDA)."
  CUDA_ROOT=$PREFIX/cuda-$CUDA_VERSION
  if [ "$MODE" = check ] && [ ! -x "$CUDA_ROOT/bin/nvcc" ]; then
    info "no hay nvcc que sirva: se instalaría CUDA $CUDA_VERSION en $CUDA_ROOT"
    say "Diagnóstico terminado (el resto se revisa después de instalar CUDA)"; exit 0
  fi
  if [ ! -x "$CUDA_ROOT/bin/nvcc" ]; then
    F=cuda_${CUDA_VERSION}_${RUNFILE[$CUDA_VERSION]}_linux.run
    info "descargando CUDA $CUDA_VERSION (~5 GB)"
    [ -s "$PREFIX/tmp/$F" ] || run cuda-download.log curl -fL -o "$PREFIX/tmp/$F" "https://developer.download.nvidia.com/compute/cuda/$CUDA_VERSION/local_installers/$F"
    info "instalando en $CUDA_ROOT"
    run cuda-install.log env -u DISPLAY -u WAYLAND_DISPLAY sh "$PREFIX/tmp/$F" --silent --toolkit --installpath="$CUDA_ROOT" \
      --no-opengl-libs --no-drm --no-man-page --tmpdir="$PREFIX/tmp"
    rm -f "$PREFIX/tmp/$F"
  fi
  NVCC=$CUDA_ROOT/bin/nvcc
fi
CUDA_ROOT=$(dirname "$(dirname "$NVCC")")
CUDA_VER=$("$NVCC" --version | sed -n 's/.*release \([0-9]*\.[0-9]*\).*/\1/p')
CUDA_MAJOR=${CUDA_VER%%.*}
export PATH=$CUDA_ROOT/bin:$PATH CUDA_HOME=$CUDA_ROOT
export LD_LIBRARY_PATH=$CUDA_ROOT/lib64:${LD_LIBRARY_PATH:-}
info "CUDA $CUDA_VER en $CUDA_ROOT"

say "Compilador y dependencias"
missing=()
for t in git curl cmake make mpicxx mpicc python3; do command -v "$t" > /dev/null || missing+=("$t"); done
if command -v cmake > /dev/null; then ver_ge "$(cmake --version | sed -n 's/cmake version //p')" 3.22 || missing+=("cmake>=3.22"); fi
[ -f /usr/include/boost/version.hpp ] || missing+=("Boost")
[ -d /usr/include/eigen3 ] || missing+=("Eigen")
compgen -G "/usr/include/hdf5/*/H5public.h" > /dev/null || [ -f /usr/include/H5public.h ] || missing+=("HDF5")
[ ${#missing[@]} -eq 0 ] || die "falta: ${missing[*]}. Corre 'bash $0 deps'."
if [ -z "${CXX_HOST:-}" ]; then
  printf '#include <concepts>\nint main(){return 0;}\n' > "$PREFIX/tmp/probe.cu"
  for c in g++ g++-13 g++-12 g++-14; do
    command -v "$c" > /dev/null || continue
    [ "$("$c" -dumpversion | cut -d. -f1)" -ge 12 ] || continue
    if nvcc -std=c++20 -ccbin "$(command -v "$c")" "$PREFIX/tmp/probe.cu" -o "$PREFIX/tmp/probe" > "$LOG/probe.log" 2>&1; then
      CXX_HOST=$c; break
    fi
  done
fi
[ -n "${CXX_HOST:-}" ] || die "nvcc $CUDA_VER no acepta ningún g++ >= 12 en C++20 (ver $LOG/probe.log). Corre 'bash $0 deps'."
CXX_HOST=$(command -v "$CXX_HOST")
CC_HOST=$(command -v "$(basename "$CXX_HOST" | sed 's/g++/gcc/')")
FC_HOST=$(command -v "$(basename "$CXX_HOST" | sed 's/g++/gfortran/')" 2> /dev/null || command -v gfortran || true)
# Los wrappers de OpenMPI y nvcc_wrapper usan el mismo compilador.
export OMPI_CC=$CC_HOST OMPI_CXX=$CXX_HOST NVCC_WRAPPER_DEFAULT_COMPILER=$CXX_HOST
if [ -n "$FC_HOST" ]; then export OMPI_FC=$FC_HOST; fi
export NVCC_APPEND_FLAGS="${NVCC_APPEND_FLAGS:-} -Wno-deprecated-gpu-targets"
CPUS=$(nproc); RAM=$(awk '/MemTotal/ {print int($2/1048576)}' /proc/meminfo)
JOBS=${JOBS:-$CPUS}
HEAVY=$(min "$JOBS" $(( RAM / 5 > 0 ? RAM / 5 : 1 )))    # nvcc con Kokkos: ~5 GB por proceso
GEN=(-G "Unix Makefiles"); if command -v ninja > /dev/null; then GEN=(-G Ninja); fi
info "g++ $($CXX_HOST -dumpfullversion) ($CXX_HOST), Fortran: ${FC_HOST:-no}"
info "$CPUS núcleos, $RAM GB de RAM: -j $JOBS, -j $HEAVY en lo que compila nvcc"
info "MPI: $(mpicxx --showme:command 2> /dev/null || mpicxx -show | cut -d' ' -f1)"
FREE=$(df -BG --output=avail "$PREFIX" | tail -1 | tr -dc 0-9)
[ "$FREE" -ge 30 ] || die "hay $FREE GB libres en $PREFIX; hacen falta unos 30."
if [ "$MODE" = check ]; then say "Diagnóstico terminado: se puede instalar"; exit 0; fi

# =============================================================================
#  2. Kokkos y Kokkos-Kernels
# =============================================================================
K=$PREFIX/kokkos-install
say "Kokkos $KOKKOS_VERSION (CUDA, $ARCH)"
if grep -qs "Kokkos_ARCH .*$ARCH" "$K"/lib*/cmake/Kokkos/KokkosConfigCommon.cmake; then
  info "ya instalado"
else
  [ -d "src/kokkos-$KOKKOS_VERSION" ] || run kokkos-download.log bash -c \
    "curl -fsSL https://github.com/kokkos/kokkos/releases/download/$KOKKOS_VERSION/kokkos-$KOKKOS_VERSION.tar.gz | tar xz -C src"
  rm -rf build-kokkos
  run kokkos-cmake.log cmake -S "src/kokkos-$KOKKOS_VERSION" -B build-kokkos "${GEN[@]}" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$K" \
    -DCMAKE_CXX_COMPILER="$PREFIX/src/kokkos-$KOKKOS_VERSION/bin/nvcc_wrapper" \
    -DCMAKE_CXX_STANDARD=20 -DBUILD_SHARED_LIBS=ON -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_COMPLEX_ALIGN=ON "-DKokkos_ARCH_$ARCH=ON" \
    -DCUDAToolkit_ROOT="$CUDA_ROOT"
  info "compilando (~5 min)"
  run kokkos-build.log cmake --build build-kokkos -j "$HEAVY"
  run kokkos-install.log cmake --install build-kokkos
fi
say "Kokkos-Kernels $KOKKOS_VERSION"
if compgen -G "$K/lib*/cmake/KokkosKernels/KokkosKernelsConfig.cmake" > /dev/null; then
  info "ya instalado"
else
  [ -d "src/kokkos-kernels-$KOKKOS_VERSION" ] || run kk-download.log bash -c \
    "curl -fsSL https://github.com/kokkos/kokkos-kernels/archive/refs/tags/$KOKKOS_VERSION.tar.gz | tar xz -C src"
  rm -rf build-kk
  run kk-cmake.log cmake -S "src/kokkos-kernels-$KOKKOS_VERSION" -B build-kk "${GEN[@]}" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$K" -DCMAKE_CXX_COMPILER="$CXX_HOST" \
    -DBUILD_SHARED_LIBS=ON -DKokkos_ROOT="$K" -DKokkosKernels_ENABLE_TESTS=OFF \
    -DKokkosKernels_ENABLE_EXAMPLES=OFF -DKokkosKernels_ENABLE_PERFTESTS=OFF
  info "compilando (20-60 min, lo más largo)"
  run kk-build.log cmake --build build-kk -j "$HEAVY"
  run kk-install.log cmake --install build-kk
fi
export LD_LIBRARY_PATH=$K/lib:$K/lib64:$LD_LIBRARY_PATH

# =============================================================================
#  3. cuDSS
# =============================================================================
say "cuDSS (CUDA $CUDA_MAJOR)"
CUDSS=$PREFIX/cudss
if compgen -G "$CUDSS/lib*/cmake/cudss/cudss-config.cmake" > /dev/null; then
  info "ya instalado"
else
  case "$(uname -m)" in x86_64) PLAT=linux-x86_64 ;; aarch64) PLAT=linux-sbsa ;; *) die "CPU $(uname -m) sin cuDSS" ;; esac
  IDX=https://developer.download.nvidia.com/compute/cudss/redist/libcudss/$PLAT/
  V=$(curl -fsSL "$IDX" | grep -o "libcudss-$PLAT-[0-9.]*_cuda$CUDA_MAJOR-archive.tar.xz" \
      | sed "s/libcudss-$PLAT-\(.*\)_cuda.*/\1/" | sort -uV | tail -1 || true)
  V=${V:-1.0.0.35}
  T=libcudss-$PLAT-${V}_cuda$CUDA_MAJOR-archive
  info "versión $V"
  run cudss-download.log bash -c "curl -fsSL '$IDX$T.tar.xz' | tar xJ -C '$PREFIX/src'"
  rm -rf "$CUDSS" && mv "$PREFIX/src/$T" "$CUDSS"
fi
export LD_LIBRARY_PATH=$CUDSS/lib:$LD_LIBRARY_PATH

# =============================================================================
#  4. PETSc con CUDA y Kokkos
# =============================================================================
P=$PREFIX/petsc
say "PETSc $PETSC_VERSION (CUDA sm_$SM, Kokkos)"
if [ -f "$P/lib/petsc/conf/petscvariables" ]; then
  info "ya instalado"
else
  [ -d "src/petsc-$PETSC_VERSION" ] || run petsc-download.log bash -c \
    "curl -fsSL https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSC_VERSION.tar.gz | tar xz -C src"
  cd "src/petsc-$PETSC_VERSION"
  OPTS=(PETSC_ARCH=arch-cuda --prefix="$P" --with-cc=mpicc --with-cxx=mpicxx
        --with-cxx-dialect=20 --with-debugging=0 --with-x=0
        COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 CUDAOPTFLAGS=-O3
        --with-cuda=1 --with-cuda-dir="$CUDA_ROOT" --with-cuda-arch="$SM"
        --with-kokkos-dir="$K" --with-kokkos-kernels-dir="$K"
        --download-metis --download-parmetis --with-make-np="$JOBS")
  if [ -n "$FC_HOST" ] && command -v mpif90 > /dev/null; then
    OPTS+=(--with-fc=mpif90 --download-fblaslapack --download-scalapack --download-mumps --download-ptscotch)
  else
    info "sin Fortran: PETSc va sin MUMPS"
    OPTS+=(--with-fc=0 --download-f2cblaslapack --download-ptscotch)
  fi
  info "configurando (10-20 min; baja MUMPS, hypre, etc.)"
  if [ "${NO_HYPRE:-0}" = 1 ] || ! try petsc-configure.log ./configure "${OPTS[@]}" --download-hypre; then
    if [ "${NO_HYPRE:-0}" != 1 ]; then info "falló con hypre (ver $LOG/petsc-configure.log); de nuevo sin hypre"; fi
    run petsc-configure-nohypre.log ./configure "${OPTS[@]}"
  fi
  info "compilando (20-40 min)"
  run petsc-make.log make PETSC_DIR="$PWD" PETSC_ARCH=arch-cuda all
  run petsc-install.log make PETSC_DIR="$PWD" PETSC_ARCH=arch-cuda install
  cd "$PREFIX"
fi
grep -q -- "--with-kokkos-dir=$K" "$P"/lib/petsc/conf/reconfigure-*.py \
  || die "el PETSc de $P no usa el Kokkos de $K."
export LD_LIBRARY_PATH=$P/lib:$LD_LIBRARY_PATH

# =============================================================================
#  5. Rodin
# =============================================================================
say "Rodin, rama $BRANCH"
if [ -d rodin/.git ]; then
  run rodin-fetch.log git -C rodin fetch origin "$BRANCH"
  run rodin-checkout.log git -C rodin checkout --detach FETCH_HEAD
else
  run rodin-clone.log git clone --single-branch --branch "$BRANCH" "$REPO" rodin
fi
run rodin-submodules.log git -C rodin submodule update --init --recursive
info "commit $(git -C rodin log -1 --format='%h %s' | cut -c1-90)"
info "configurando"
run rodin-cmake.log env PETSC_DIR="$P" cmake -S rodin -B rodin/build-cuda "${GEN[@]}" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$CC_HOST" -DCMAKE_CXX_COMPILER="$CXX_HOST" \
  -DCMAKE_PREFIX_PATH="$K;$CUDSS;$P" -DPETSc_DIR="$P" -DCUDAToolkit_ROOT="$CUDA_ROOT" \
  -DRODIN_USE_KOKKOS=ON -DRODIN_USE_CUDSS=ON -DRODIN_USE_PETSC=ON -DRODIN_USE_MPI=ON \
  -DRODIN_BUILD_EXAMPLES=ON -DRODIN_BUILD_UNIT_TESTS=OFF -DRODIN_BUILD_MANUFACTURED_TESTS=OFF \
  -DRODIN_BUILD_BENCHMARKS=OFF -DRODIN_INSTALL_RESOURCES=OFF -DRODIN_MULTITHREADED=OFF \
  -DRODIN_CODE_COVERAGE=OFF
BT=()
if [ "$TARGETS" != all ]; then read -r -a T <<< "$TARGETS"; BT=(--target "${T[@]}"); fi
info "compilando ${TARGETS} (20-40 min)"
run rodin-build.log cmake --build rodin/build-cuda -j "$HEAVY" "${BT[@]}"
BIN=$PREFIX/rodin/build-cuda/examples/Kokkos/Heart
if ldd "$BIN/KokkosAtriumRheologyPETSc" 2>&1 | grep 'not found'; then die "faltan bibliotecas compartidas."; fi

# =============================================================================
#  6. Entorno, malla y prueba
# =============================================================================
cat > "$PREFIX/env.sh" << ENVEOF
# source $PREFIX/env.sh antes de correr
export PATH=$CUDA_ROOT/bin:\$PATH
export LD_LIBRARY_PATH=$K/lib:$K/lib64:$CUDSS/lib:$P/lib:$CUDA_ROOT/lib64:\${LD_LIBRARY_PATH:-}
export PETSC_OPTIONS="-use_gpu_aware_mpi 0 -ksp_type bcgsl -ksp_bcgsl_ell 4 -pc_type jacobi"
export RODIN_BIN=$BIN
export ATRIUM_MESH=$PREFIX/meshes/Coarse.mesh
ENVEOF
# shellcheck disable=SC1091
source "$PREFIX/env.sh"

say "Malla Coarse.mesh"
mkdir -p meshes test
if echo "$MESH_SHA  meshes/Coarse.mesh" | sha256sum -c --status 2> /dev/null; then
  info "ya está"
elif try mesh-download.log curl -fL -o meshes/Coarse.mesh "$MESH_URL" \
     && echo "$MESH_SHA  meshes/Coarse.mesh" | sha256sum -c --status; then
  info "descargada ($(du -h meshes/Coarse.mesh | cut -f1))"
else
  rm -f meshes/Coarse.mesh
  info "no pude descargarla: cópiala a $PREFIX/meshes/Coarse.mesh (no se hace la prueba)"
  SKIP_TEST=1
fi

if [ "${SKIP_TEST:-0}" != 1 ]; then
  for f in N TVE; do
    say "Prueba: $f, 20 pasos en Coarse.mesh"
    nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 500 > "$LOG/test_${f}_mem.txt" &
    MON=$!
    if ! (cd test && timeout 1800 "$BIN/KokkosAtriumRheologyPETSc" "$ATRIUM_MESH" "fluid=$f" t_end=0.02 \
          output_every=10 "output=test_$f" > "$LOG/test_$f.log" 2>&1); then
      kill $MON 2> /dev/null || true
      tail -20 "$LOG/test_$f.log"; die "la prueba $f falló (ver $LOG/test_$f.log)"
    fi
    kill $MON 2> /dev/null || true
    grep -E '^Step 20|^Rheology' "$LOG/test_$f.log" | cut -c1-230 | sed 's/^/    /'
    info "pico de memoria GPU: $(sort -n "$LOG/test_${f}_mem.txt" | tail -1) MiB"
  done
fi

say "Listo en $((SECONDS / 60)) min"
cat << EOF
    GPU $GPU, Kokkos_ARCH_$ARCH, CUDA $CUDA_VER, g++ $($CXX_HOST -dumpfullversion)
    Para correr:
      source $PREFIX/env.sh
      mkdir -p ~/runs && cd ~/runs
      nohup \$RODIN_BIN/KokkosAtriumRheologyPETSc \$ATRIUM_MESH fluid=TVE cycles=3 output=TVE > TVE.log 2>&1 &
EOF
