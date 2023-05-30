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
      constexpr
      Geometry::Polytope::Geometry getGeometry() const
      {
        return m_g;
      }

      inline
      constexpr
      auto getDOF(size_t i) const
      {
        return getDOFs().col(i);
      }

      virtual ~FiniteElementBase() = default;

      inline
      constexpr
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      inline
      constexpr
      auto getBasis(size_t local) const
      {
        return static_cast<const Derived&>(*this).getBasis(local);
      }

      inline
      constexpr
      auto getBasis(const Math::Vector& rc) const
      {
        return static_cast<const Derived&>(*this).getBasis(rc);
      }

      inline
      constexpr
      const Math::Matrix& getDOFs() const
      {
        return static_cast<const Derived&>(*this).getDOFs();
      }

    private:
      const Geometry::Polytope::Geometry m_g;
  };

  template <class FE>
  class FiniteElementProduct : public FiniteElementBase<FiniteElementProduct<FE>>
  {
    public:
      using Parent = FiniteElementBase<FiniteElementProduct<FE>>;

      constexpr
      FiniteElementProduct(Geometry::Polytope::Geometry g, size_t vdim)
        : Parent(g), m_vdim(vdim), m_fe(g)
      {}

      constexpr
      FiniteElementProduct(const FiniteElementProduct& other)
        : Parent(other), m_vdim(other.m_vdim), m_fe(other.m_fe)
      {}

      constexpr
      FiniteElementProduct(FiniteElementProduct&& other)
        : Parent(std::move(other)), m_vdim(other.m_vdim), m_fe(other.m_fe)
      {}

      inline
      constexpr
      size_t getCount() const
      {
        return m_fe.getCount() * m_vdim;
      }

      inline
      auto getBasis(size_t local) const
      {
        return 0;
      }

      inline
      auto getBasis(const Math::Vector& rc) const
      {
        const size_t rdim = Geometry::Polytope::getGeometryDimension(this->getGeometry());
        Math::Matrix res(rdim, getCount());
        for (size_t i = 0; i < getCount(); i++)
          res.col(i) = m_fe.getBasis(rc);
        return res;
      }

      inline
      const Math::Matrix& getDOFs() const
      {
      }

    private:
      const size_t m_vdim;
      const FE m_fe;
  };
}

#endif

