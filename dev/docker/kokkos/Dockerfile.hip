# syntax=docker/dockerfile:1
#
# Dependencies of the Kokkos CI lane on AMD GPUs: ROCm 10.1, UCX and Open MPI
# with ROCm, Boost with MPI, Kokkos 5 with HIP, PETSc with Kokkos,
# Kokkos-Kernels, hypre, MUMPS, SuperLU_DIST, STRUMPACK and SLATE on the
# device, and UMFPACK for the first factorization of RocSolverRF. Rodin
# itself is built by the CI job. See README.md.
#
# ROCm 10.1 packages its device libraries per GPU target, and Kokkos compiles
# for one architecture: GFX selects both. The default is the Radeon AI PRO
# R9700 of zephyr (RDNA 4, gfx1201).

ARG UBUNTU_IMAGE=ubuntu:24.04
FROM ${UBUNTU_IMAGE}

ARG JOBS=36
ARG GFX=gfx1201
ARG KOKKOS_ARCH=AMD_GFX1201
# RDNA runs waves of 32 lanes, CDNA of 64. hypre detects it with rocm-smi,
# which ROCm 10.1 lacks, and falls back to 64: BoomerAMG then faults on RDNA.
ARG GPU_WARP_SIZE=32
ARG ROCM_VERSION=10.1
ARG UCX_VERSION=1.22.0
ARG OMPI_VERSION=5.0.11
ARG BOOST_VERSION=1.90.0
ARG KOKKOS_VERSION=5.2.2
ARG PETSC_VERSION=3.26.0

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      autoconf automake bison bzip2 ca-certificates cmake curl file flex g++ gcc \
      gfortran git gpg libhdf5-dev libnuma-dev libsuitesparse-dev libeigen3-dev \
      libtool make ninja-build patch pkg-config python3 xz-utils zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- ROCm ------------------------------------------------------------------
RUN install -d /etc/apt/keyrings \
    && curl -fsSL https://stable.repo.amd.com/rocm/gpg/packages.gpg \
      | gpg --dearmor -o /etc/apt/keyrings/amdrocm.gpg \
    && printf '%s\n' 'Types: deb' \
         'URIs: https://stable.repo.amd.com/rocm/core/packages/ubuntu2404/' \
         'Suites: stable' 'Components: main' 'Architectures: amd64' \
         'Signed-By: /etc/apt/keyrings/amdrocm.gpg' \
         > /etc/apt/sources.list.d/amdrocm.sources \
    && apt-get update \
    && apt-get install -y --no-install-recommends amdrocm-core-dev${ROCM_VERSION}-${GFX} \
    && rm -rf /var/lib/apt/lists/*

ENV ROCM_PATH=/opt/rocm
ENV PATH=/opt/rocm/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/rocm/lib
# Kokkos-Kernels finds rocPRIM and rocThrust through it.
ENV CMAKE_PREFIX_PATH=/opt/rocm

# ---- UCX and Open MPI, ROCm aware ----------------------------------------
RUN curl -fsSL https://github.com/openucx/ucx/releases/download/v${UCX_VERSION}/ucx-${UCX_VERSION}.tar.gz \
      | tar xz -C /tmp \
    && cd /tmp/ucx-${UCX_VERSION} \
    && ./contrib/configure-release --prefix=/opt/ompi --with-rocm=${ROCM_PATH} \
         --without-cuda --without-verbs --without-rdmacm --without-knem --without-xpmem \
    && make -j${JOBS} && make install \
    && rm -rf /tmp/ucx-${UCX_VERSION}

RUN curl -fsSL https://download.open-mpi.org/release/open-mpi/v${OMPI_VERSION%.*}/openmpi-${OMPI_VERSION}.tar.bz2 \
      | tar xj -C /tmp \
    && cd /tmp/openmpi-${OMPI_VERSION} \
    && ./configure --prefix=/opt/ompi --with-ucx=/opt/ompi --with-rocm=${ROCM_PATH} \
         --without-verbs --without-ofi --with-hwloc=internal --with-libevent=internal \
         --with-pmix=internal --with-prrte=internal --enable-mpi-fortran=usempif08 \
    && make -j${JOBS} && make install \
    && rm -rf /tmp/openmpi-${OMPI_VERSION}

ENV PATH=/opt/ompi/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/ompi/lib:${LD_LIBRARY_PATH}
# GitHub runs container jobs as root.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
# The CI processes share one GPU.
ENV OMPI_MCA_rmaps_base_oversubscribe=1
# The PCI component of hwloc loads the system zlib after the one of ROCm and
# randomly aborts MPI_Init.
ENV HWLOC_COMPONENTS=-pci

# ---- Boost with MPI --------------------------------------------------------
RUN curl -fsSL https://archives.boost.io/release/${BOOST_VERSION}/source/boost_$(echo ${BOOST_VERSION} | tr . _).tar.bz2 \
      | tar xj -C /tmp \
    && cd /tmp/boost_$(echo ${BOOST_VERSION} | tr . _) \
    && ./bootstrap.sh --prefix=/opt/boost \
         --with-libraries=atomic,container,filesystem,graph,mpi,serialization,system \
    && echo "using mpi : /opt/ompi/bin/mpicxx ;" >> project-config.jam \
    && ./b2 -j${JOBS} cxxstd=20 variant=release link=shared threading=multi install \
    && rm -rf /tmp/boost_*

ENV LD_LIBRARY_PATH=/opt/boost/lib:${LD_LIBRARY_PATH}

# ---- Kokkos 5 --------------------------------------------------------------
RUN curl -fsSL https://github.com/kokkos/kokkos/releases/download/${KOKKOS_VERSION}/kokkos-${KOKKOS_VERSION}.tar.gz \
      | tar xz -C /tmp \
    && cmake -S /tmp/kokkos-${KOKKOS_VERSION} -B /tmp/kokkos-build -G Ninja \
         -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/kokkos \
         -DCMAKE_CXX_COMPILER=/opt/rocm/bin/hipcc \
         -DCMAKE_CXX_STANDARD=20 -DBUILD_SHARED_LIBS=ON \
         -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_ROCTHRUST=ON \
         -DKokkos_ENABLE_COMPLEX_ALIGN=ON -DKokkos_ARCH_${KOKKOS_ARCH}=ON \
    && cmake --build /tmp/kokkos-build -j${JOBS} && cmake --install /tmp/kokkos-build \
    && rm -rf /tmp/kokkos-build /tmp/kokkos-${KOKKOS_VERSION}

ENV LD_LIBRARY_PATH=/opt/kokkos/lib:${LD_LIBRARY_PATH}

# ---- PETSc -----------------------------------------------------------------
# Kokkos-Kernels 5.2.2 names rocsparse_indextype_u16, which ROCm 10 removed.
COPY kokkos-kernels-rocsparse-u16.patch /tmp/
RUN curl -fsSL https://github.com/kokkos/kokkos-kernels/archive/refs/tags/${KOKKOS_VERSION}.tar.gz \
      | tar xz -C /tmp \
    && patch -d /tmp/kokkos-kernels-${KOKKOS_VERSION} -p1 < /tmp/kokkos-kernels-rocsparse-u16.patch \
    && tar czf /tmp/kokkos-kernels.tar.gz -C /tmp kokkos-kernels-${KOKKOS_VERSION} \
    && rm -rf /tmp/kokkos-kernels-${KOKKOS_VERSION}

# STRUMPACK links roc::rocthrust, whose hip::device adds -x hip to every C++
# source: it needs the clang++ of ROCm. SLATE must use the same compiler, as
# its headers embed omp_nest_lock_t, whose size differs between libgomp and
# libomp.
RUN curl -fsSL https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz \
      | tar xz -C /tmp \
    && cd /tmp/petsc-${PETSC_VERSION} \
    && ./configure PETSC_ARCH=arch-hip --prefix=/opt/petsc \
         --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
         --with-hip=1 --with-hip-dir=${ROCM_PATH} --with-hip-arch=${GFX} \
         --with-cxx-dialect=20 --with-debugging=0 --with-openmp=1 \
         --with-kokkos-dir=/opt/kokkos --download-kokkos-kernels=/tmp/kokkos-kernels.tar.gz \
         --download-fblaslapack --download-hypre --download-metis --download-mumps \
         --download-parmetis --download-scalapack --download-superlu_dist \
         --download-strumpack --download-slate \
         --download-hypre-configure-arguments=--with-gpu-warp-size=${GPU_WARP_SIZE} \
         --download-strumpack-cmake-arguments=-DCMAKE_CXX_COMPILER=/opt/rocm/lib/llvm/bin/clang++ \
         --download-slate-cmake-arguments=-DCMAKE_CXX_COMPILER=/opt/rocm/lib/llvm/bin/clang++ \
         --with-make-np=${JOBS} \
    && make PETSC_DIR=$PWD PETSC_ARCH=arch-hip all \
    && make PETSC_DIR=$PWD PETSC_ARCH=arch-hip install \
    && rm -rf /tmp/petsc-${PETSC_VERSION} /tmp/kokkos-kernels.tar.gz

ENV PETSC_DIR=/opt/petsc
ENV LD_LIBRARY_PATH=/opt/petsc/lib:${LD_LIBRARY_PATH}
ENV CMAKE_PREFIX_PATH=/opt/petsc:/opt/kokkos:/opt/boost:/opt/rocm
ENV RODIN_KOKKOS_BACKEND=hip
