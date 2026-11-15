/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VIEW_H
#define RODIN_KOKKOS_VIEW_H

/**
 * @file
 * @brief Device-callable representations of Rodin objects.
 */

namespace Rodin::Kokkos
{
  /**
   * @brief Execution-facing representation of a Rodin object.
   * @tparam T Represented Rodin type.
   *
   * Specializations live in the corresponding Kokkos module. Geometry views
   * borrow resident arrays and expire on owner replacement; constitutive views
   * copy material parameters. Construction does not synchronize or transfer an
   * owning object. Unsupported types have no default implementation.
   * This object interface is distinct from ::Kokkos::View array storage.
   *
   * | Represented type | Representation |
   * |------------------|----------------|
   * | Geometry::Mesh<Context::Kokkos<Distribution, Backend>> | Typed borrowed mesh. |
   * | Geometry::Connectivity<Context::Kokkos<Context::Local, Backend>> | Typed borrowed connectivity. |
   * | Geometry::MeshBase, Geometry::ConnectivityBase | Backend-independent borrowed kernel interfaces. |
   * | Geometry::Polytope, Geometry::Point | Borrowed entity and geometric evaluation. |
   * | IndexVector, Geometry::Incidence | Read-only index range and CSR incidence. |
   * | Solid::SaintVenantKirchhoff, Solid::Yeoh, Solid::MooneyRivlin | Copied material parameters and constitutive evaluation. |
   */
  template <class T>
  struct View;
}

#endif
