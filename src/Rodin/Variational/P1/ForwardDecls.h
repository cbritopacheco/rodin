/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_FORWARDDECLS_H
#define RODIN_VARIATIONAL_P1_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Degree 1 Lagrange element
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * P1 class, please see
   * <a href="_variational_2_p1_8h.html">P1Specializations</a>.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P1Element "P1Element<Scalar>" | Scalar-valued continuous piecewise linear element. |
   * | @ref P1Element "P1Element<SpatialVector<Scalar>>" | Vector-valued continuous piecewise linear element. |
   *
   * @see <a href="_p1_element_8h.html">P1ElementSpecializations</a>
   */
  template <class Range>
  class P1Element;

  /**
   * @brief Degree 1 Lagrange finite element space
   * @tparam Range Range value type
   * @tparam Context Context type
   * @tparam Args Additional arguments
   *
   * Represents the finite element space composed of continuous, piecewise
   * linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h)^d = \{ v \in C^0(\mathcal{T}_h)^d \mid v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ ,
   * @f]
   * for a given vector dimension @f$ d \in \mathbb{N} @f$.
   *
   * @note For an overview of all the possible specializations of the
   * P1 class, please see
   * <a href="_variational_2_p1_8h.html">P1Specializations</a>.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P1 "P1<Scalar, Mesh<Context::Local>>" | Scalar-valued local-mesh continuous piecewise linear space. |
   * | @ref P1 "P1<SpatialVector<Scalar>, Mesh<Context::Local>>" | Vector-valued local-mesh continuous piecewise linear space. |
   * | @ref P1 "P1<Range, Mesh<Context::MPI>>" | Scalar or vector-valued distributed-mesh continuous piecewise linear space. |
   *
   * @see <a href="_variational_2_p1_8h.html">P1Specializations</a>
   */
  template <class Range, class Mesh>
  class P1;

  /**
   * @brief Alias for P1Element<Real>
   */
  using RealP1Element = P1Element<Real>;

  /**
   * @brief Alias for P1Element<Complex>
   */
  using ComplexP1Element = P1Element<Complex>;

  /**
   * @brief Alias for P1Element<Math::SpatialVector<Scalar>>
   * @tparam Scalar Scalar type for vector components
   */
  template <class Scalar>
  using VectorP1Element = P1Element<Math::SpatialVector<Scalar>>;

  /**
   * @brief Alias for real vector P1Element
   */
  using RealVectorP1Element = P1Element<Math::SpatialVector<Real>>;
}

#endif
