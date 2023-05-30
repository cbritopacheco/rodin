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
#include "MFEM.h"

namespace Rodin::Variational
{
  class FiniteElementOrder
  {
    public:
      explicit
      constexpr
      FiniteElementOrder(size_t v)
        : m_v(v)
      {}

      constexpr
      FiniteElementOrder(const FiniteElementOrder&) = default;

      constexpr
      FiniteElementOrder(FiniteElementOrder&&) = default;

      inline
      constexpr
      operator size_t() const
      {
        return m_v;
      }

    private:
      const size_t m_v;
  };

  template <class Derived>
  class FiniteElementBase
  {
    public:
      constexpr
      FiniteElementBase(Geometry::Polytope::Geometry g)
        : m_g(g)
      {}

      inline
      auto getBasis(size_t local) const
      {
        const size_t count = getCount();
        return Math::Matrix::Identity(count, count).row(local);
      }

      inline
      constexpr
      Geometry::Polytope::Geometry getGeometry() const
      {
        return m_g;
      }

      inline
      auto getDOF(size_t i) const
      {
        return getDOFs().col(i);
      }

      virtual ~FiniteElementBase() = default;

      inline
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      inline
      auto getBasis(const Math::Vector& rc) const
      {
        return static_cast<const Derived&>(*this).getBasis(rc);
      }

      inline
      const Math::Matrix& getDOFs() const
      {
        return static_cast<const Derived&>(*this).getDOFs();
      }

    private:
      const Geometry::Polytope::Geometry m_g;
  };

  class FiniteElementProduct
  {
    public:
    private:
  };
}

#endif

