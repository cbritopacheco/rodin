/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include <unordered_map>

#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Simplex.h"
#include "Rodin/Geometry/SimplexTransformation.h"

#include "TensorBasis.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FiniteElements Supported finite elements
   * @brief List of finite elements already implemented.
   */

  enum FiniteElementMapping
  {
    Identity,
    CovariantPiola,
    ContravariantPiola,
    DoubleCovariantPiola,
    DoubleContravariantPiola
  };

  /**
   * @brief Base class for finite elements.
   */
  template <class Derived>
  class FiniteElementBase
  {
    public:
      constexpr
      FiniteElementBase(Geometry::Polytope::Geometry g)
        : m_g(g)
      {}

      virtual ~FiniteElementBase() = default;

      inline
      constexpr
      Geometry::Polytope::Geometry getGeometry() const
      {
        return m_g;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       */
      inline
      constexpr
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      /**
       * @brief Gets the i-th degree of freedom on the finite element.
       */
      inline
      constexpr
      auto getNode(size_t i) const
      {
        return getNodes().col(i);
      }

      inline
      constexpr
      const Math::Matrix& getNodes() const
      {
        return static_cast<const Derived&>(*this).getNodes();
      }

      inline
      constexpr
      const auto& getBasis(size_t i) const
      {
        return static_cast<const Derived&>(*this).getBasis(i);
      }

      inline
      constexpr
      auto getBasis(const Math::Vector& r) const
      {
        return static_cast<const Derived&>(*this).getBasis(r);
      }

      inline
      constexpr
      const auto& getLinearForm(size_t i) const
      {
        return static_cast<const Derived&>(*this).getLinearForm(i);
      }

      inline
      constexpr
      FiniteElementMapping getMapping() const
      {
        return static_cast<const Derived&>(*this).getMapping();
      }

    private:
      const Geometry::Polytope::Geometry m_g;
  };
}

#endif

