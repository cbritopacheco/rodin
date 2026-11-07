/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Context module.
 */
#ifndef RODIN_CONTEXT_FORWARDDECLS_H
#define RODIN_CONTEXT_FORWARDDECLS_H

#include <type_traits>

namespace Rodin::Context
{
  class Base;

  class Local;

  class MPI;

  class CUDA;

  class HIP;

  class SYCL;

  template <class Distribution, class Backend>
  class Kokkos;

  /// @brief Whether a context executes on a Kokkos device.
  template <class C>
  struct IsKokkos : std::false_type
  {};

  /// @brief Kokkos contexts execute on a Kokkos device.
  template <class Distribution, class Backend>
  struct IsKokkos<Kokkos<Distribution, Backend>> : std::true_type
  {};

  /// @brief Whether a context runs its mesh in a single process.
  template <class C>
  struct IsLocal : std::is_same<C, Local>
  {};

  /// @brief Kokkos contexts of a local mesh run in a single process.
  template <class Backend>
  struct IsLocal<Kokkos<Local, Backend>> : std::true_type
  {};

  /// @brief Whether a context distributes its mesh over MPI processes.
  template <class C>
  struct IsMPI : std::is_same<C, MPI>
  {};

  /// @brief Kokkos contexts of a distributed mesh run over MPI processes.
  template <class Backend>
  struct IsMPI<Kokkos<MPI, Backend>> : std::true_type
  {};
}

#endif
