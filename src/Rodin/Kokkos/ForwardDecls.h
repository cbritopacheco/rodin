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
   * @defgroup KokkosAssemblySpecializations Kokkos Assembly Template Specializations
   * @brief Template specializations of the Kokkos assembly class.
   * @see <a href="class_rodin_1_1_assembly_1_1_kokkos.html">Kokkos</a>
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Kokkos "Kokkos<::Rodin::Math::SparseMatrix<Real>, BilinearForm<...>>" | Sparse bilinear form into a host matrix. |
   * | @ref Kokkos "Kokkos<::Rodin::Math::Vector<Real>, LinearForm<...>>" | Linear form into a host vector. |
   * | @ref Kokkos "Kokkos<::Rodin::Math::LinearSystem<Matrix, ::Rodin::Math::Vector<Real>>, Problem<...>>" | Single field or mixed problem into a host or device CSR system. |
   * | @ref Kokkos "Kokkos<IndexMap<...>, DirichletBC<...>>" | Value and identification Dirichlet conditions. |
   * | @ref Kokkos "Kokkos<::Mat, BilinearForm<...>>" | Bilinear form into a PETSc Kokkos matrix. |
   * | @ref Kokkos "Kokkos<::Vec, LinearForm<...>>" | Linear form into a PETSc Kokkos vector. |
   * | @ref Kokkos "Kokkos<PETSc::Math::LinearSystem, Problem<...>>" | Single field or multi-field problem into a PETSc Kokkos system. |
   */

  /**
   * @brief Device assembly on resident Kokkos meshes.
   * @tparam LinearAlgebraType Target linear algebra type.
   * @tparam Object Variational form type to assemble.
   *
   * @see @ref KokkosAssemblySpecializations "Kokkos Assembly Template Specializations"
   */
  template <class LinearAlgebraType, class Object>
  class Kokkos;
}

namespace Rodin::Kokkos::Math
{
  /// @brief Owning backend coefficient vector.
  template <class Scalar, class Backend>
  class Vector;

  /// @brief Owning backend compressed row matrix.
  template <class Scalar, class Backend>
  class SparseMatrix;
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
