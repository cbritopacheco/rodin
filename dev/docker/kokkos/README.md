# Kokkos CI images

Two images hold the dependencies of the Kokkos CI lane, one per GPU vendor.
Kokkos compiles for one GPU architecture and PETSc for one of CUDA and HIP,
so one image cannot cover both cards of a host.

| Image | GPU | Toolchain |
|---|---|---|
| `rodin-kokkos-ci:cuda-sm52` | NVIDIA GeForce GTX 980 Ti (Maxwell, sm_52) | CUDA 12.9, GCC 13 |
| `rodin-kokkos-ci:hip-gfx1201` | AMD Radeon AI PRO R9700 (RDNA 4) | ROCm 10.1, ROCm clang |

Both contain UCX 1.22 and Open MPI 5.0 aware of the GPU, Boost 1.90 with
MPI, Kokkos 5.2.2, and PETSc 3.26 with Kokkos-Kernels, hypre, MUMPS,
SuperLU_DIST, ParMETIS and ScaLAPACK. The CUDA image adds cuDSS; the HIP
image adds rocSOLVER, UMFPACK, STRUMPACK and SLATE. Rodin itself is built by
the CI job with `ci-test.sh`.

The default targets are the two cards of zephyr. PETSc and hypre embed the
machine code of one architecture and no PTX, so the image runs on that
architecture only: another GPU needs its own image, with the build
arguments `KOKKOS_ARCH` (for instance `TURING75`) and `CUDA_ARCH` (`75`); below sm_60, SuperLU_DIST runs on the host. The HIP image takes `GFX` (`gfx90a`), `KOKKOS_ARCH`
(`AMD_GFX90A`) and `GPU_WARP_SIZE` (64 on CDNA, 32 on RDNA).

## Build

The CUDA image needs no GPU to build. A PETSc build takes most of the time;
`JOBS` sets the parallelism (36 by default).

```sh
cd dev/docker/kokkos
docker build -f Dockerfile.cuda -t rodin-kokkos-ci:cuda-sm52 .
docker build -f Dockerfile.hip -t rodin-kokkos-ci:hip-gfx1201 .
```

To move an image to the CI host without a registry:

```sh
docker save rodin-kokkos-ci:cuda-sm52 | zstd | ssh zephyr 'zstd -d | docker load'
```

## Run

NVIDIA needs the NVIDIA Container Toolkit on the host; AMD needs the kernel
devices and membership of the `render` and `video` groups.

```sh
docker run --rm --gpus all -v "$PWD:/rodin" -w /rodin \
  rodin-kokkos-ci:cuda-sm52 dev/docker/kokkos/ci-test.sh
docker run --rm --device /dev/kfd --device /dev/dri --group-add render \
  --group-add video -v "$PWD:/rodin" -w /rodin \
  rodin-kokkos-ci:hip-gfx1201 dev/docker/kokkos/ci-test.sh
```

On a host with several GPUs, `--gpus device=0`, or `HIP_VISIBLE_DEVICES`,
selects one. `ci-test.sh [source-dir] [build-dir]` configures Rodin for the
backend of the image, builds the Kokkos unit tests with the cuDSS or
RocSolverRF test and the Kokkos examples, and runs the tests; the MPI tests run on 1, 2 and 4 processes,
which share the GPU. Logs land in `<build-dir>/test-logs`.

`ci-test.sh --build-only` stops after the build, on a host without GPU,
and compiles at -O0 on the host and the device.
The workflow `KokkosBuild.yml` runs it on GitHub hosted runners, with the
images pulled from `ghcr.io/<owner>/rodin-kokkos-ci`, or from the registry
named by the repository variable `RODIN_KOKKOS_CI_IMAGES`. To publish them:

```sh
docker tag rodin-kokkos-ci:cuda-sm52 ghcr.io/<owner>/rodin-kokkos-ci:cuda-sm52
docker push ghcr.io/<owner>/rodin-kokkos-ci:cuda-sm52
```

The package must be public, or readable by the repository.

The submodules must be checked out (`git submodule update --init
--recursive`), GoogleTest among them.
