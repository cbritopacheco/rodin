/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

#include <Eigen/QR>

#include "Rodin/Geometry/Simplex.h"

#include "SimplexTransformation.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Polytope isoparametric transformation.
   */
  template <class FE>
  class IsoparametricTransformation final : public PolytopeTransformation
  {
    static_assert(std::is_same_v<typename FE::RangeType, Scalar>,
        "Type of finite element must be scalar valued.");

    public:
      using Parent = PolytopeTransformation;
      using Parent::transform;
      using Parent::jacobian;
      using Parent::inverse;

      IsoparametricTransformation(Math::SpatialMatrix&& pm, FE&& fe)
        : Parent(Polytope::getGeometryDimension(fe.getGeometry()), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * pm : sdim x dof
       */
      IsoparametricTransformation(const Math::SpatialMatrix& pm, const FE& fe)
        : Parent(Polytope::getGeometryDimension(fe.getGeometry()), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(Math::SpatialMatrix&& pm, const FE& fe)
        : Parent(Polytope::getGeometryDimension(fe.getGeometry()), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const Math::SpatialMatrix& pm, FE&& fe)
        : Parent(Polytope::getGeometryDimension(fe.getGeometry()), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const IsoparametricTransformation& other)
        : Parent(other),
          m_fe(other.m_fe),
          m_pm(other.m_pm)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(IsoparametricTransformation&& other)
        : Parent(std::move(other)),
          m_fe(std::move(other.m_fe)),
          m_pm(std::move(other.m_pm))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      inline
      void transform(const Math::SpatialVector& rc, Math::SpatialVector& pc) const override
      {
        const size_t pdim = getPhysicalDimension();
        assert(rc.size() >= 0);
        assert(static_cast<size_t>(rc.size()) == getReferenceDimension());
        pc.resize(pdim);
        pc.setZero();
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          assert(pc.size() == m_pm.col(local).size());
          pc.noalias() += m_pm.col(local) * m_fe.getBasis(local)(rc);
        }
      }

      inline
      void jacobian(const Math::SpatialVector& rc, Math::SpatialMatrix& res) const override
      {
        const size_t rdim = getReferenceDimension();
        assert(rc.size() >= 0);
        assert(static_cast<size_t>(rc.size()) == rdim);
        const size_t pdim = getPhysicalDimension();
        res.resize(pdim, rdim);
        res.setZero();
        Math::SpatialVector gradient;
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          m_fe.getGradient(local)(gradient, rc);
          for (size_t i = 0; i < rdim; i++)
          {
            assert(res.col(i).size() == m_pm.col(local).size());
            res.col(i).noalias() += m_pm.col(local) * gradient.coeff(i);
          }
        }
      }


      inline
      const Math::SpatialMatrix& getPointMatrix() const
      {
        return m_pm;
      }

    private:
      Math::SpatialMatrix m_pm;
      FE m_fe;
  };
}

#endif
