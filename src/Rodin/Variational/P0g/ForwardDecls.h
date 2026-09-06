/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the global-constant (P0g) finite element
 * space.
 */
#ifndef RODIN_VARIATIONAL_P0G_FORWARDDECLS_H
#define RODIN_VARIATIONAL_P0G_FORWARDDECLS_H

namespace Rodin::Variational
{
  /**
   * @brief Forward declaration of global P0 finite elements.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P0gElement "P0gElement<Scalar>" | Scalar global constant element with one basis function @f$1@f$. |
   * | @ref P0gElement "P0gElement<Math::SpatialVector<Scalar>>" | Vector-valued global constant element with one constant basis vector per component. |
   *
   * @tparam Range Value range of the element.
   */
  template <class Range>
  class P0gElement;

  /**
   * @brief Forward declaration of global constant finite element spaces.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P0g "P0g<Real, Mesh<Context::Local>>" | Scalar global constant space on a local mesh with one degree of freedom. |
   * | @ref P0g "P0g<Math::SpatialVector<Real>, Mesh<Context::Local>>" | Vector-valued global constant space on a local mesh with one degree of freedom per component. |
   *
   * @tparam Range Value range of the finite element space.
   * @tparam Mesh Mesh type on which the space is defined.
   */
  template <class Range, class Mesh>
  class P0g;
}

#endif
