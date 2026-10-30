/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_FORWARDDECLS_H
#define RODIN_KOKKOS_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations for Kokkos execution representations.
 */

#include "Rodin/Context/ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @brief Resident mesh traversal for Kokkos assembly.
   * @tparam Mesh Mesh type supplying resident entities.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref KokkosIteration "KokkosIteration<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>" | Indexed traversal of a local resident mesh region. |
   */
  template <class Mesh>
  class KokkosIteration;

  /**
   * @brief Device assembly on resident Kokkos meshes.
   * @tparam LinearAlgebraType Target linear algebra type.
   * @tparam Object Variational form type to assemble.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Kokkos "Kokkos<Math::SparseMatrix<Real>, BilinearForm<...>>" | Device assembly of one sparse bilinear form. |
   */
  template <class LinearAlgebraType, class Object>
  class Kokkos;
}

namespace Rodin::Kokkos::Geometry
{
  class Vector;

  class Matrix;

  struct IndexView;

  struct IncidenceDescriptor;

  struct IncidenceView;

  struct ConnectivityView;

  struct MeshView;

  struct PolytopeView;

  class PointView;
}

#endif
