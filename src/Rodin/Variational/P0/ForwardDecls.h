/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_FORWARDDECLS_H
#define RODIN_VARIATIONAL_P0_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Variational
{

  /**
   * @brief Degree 0 Lagrange element
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * P0 class, please see
   * <a href="_variational_2_p0_8h.html">P0Specializations</a>.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P0Element "P0Element<Scalar>" | Scalar-valued discontinuous piecewise constant element. |
   * | @ref P0Element "P0Element<SpatialVector<Scalar>>" | Vector-valued discontinuous piecewise constant element. |
   *
   * @see <a href="_p0_element_8h.html">P0ElementSpecializations</a>
   */
  template <class Range>
  class P0Element;

  /**
   * @brief Degree 0 Lagrange finite element space
   * @tparam Range Range value type
   * @tparam Context Context type
   * @tparam Args Additional arguments
   *
   * Represents the finite element space composed of discontinuous, piecewise
   * constant functions:
   * @f[
   *  \mathbb{P}_0 (\mathcal{T}_h)^d = \{ v : \mathcal{T}_h \rightarrow \mathbb{R}^d : v|_{\tau} \in \mathbb{P}_0(\tau), \ \tau \in \mathcal{T}_h \} \ ,
   * @f]
   * for a given vector dimension @f$ d \in \mathbb{N} @f$.
   *
   * @note For an overview of all the possible specializations of the
   * P0 class, please see
   * <a href="_variational_2_p0_8h.html">P0Specializations</a>.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P0 "P0<Real, Mesh<Context::Local>>" | Scalar-valued local-mesh discontinuous piecewise constant space. |
   * | @ref P0 "P0<Range, Mesh<Context::MPI>>" | Scalar or vector-valued distributed-mesh discontinuous piecewise constant space. |
   *
   * @see <a href="_variational_2_p0_8h.html">P0Specializations</a>
   */
  template <class Range, class Mesh>
  class P0;

  /**
   * @brief Alias for P0Element<Real>
   */
  using RealP0Element = P0Element<Real>;

  /**
   * @brief Alias for P0Element<Complex>
   */
  using ComplexP0Element = P0Element<Complex>;

  /**
   * @brief Alias for P0Element<Math::SpatialVector<Scalar>>
   */
  template <class Scalar>
  using VectorP0Element = P0Element<Math::SpatialVector<Scalar>>;

  /**
   * @brief Alias for real vector P0Element
   */
  using RealVectorP0Element = VectorP0Element<Real>;

  /**
   * @brief Alias for complex vector P0Element
   */
  using ComplexVectorP0Element = VectorP0Element<Complex>;
}

#endif
