/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "PolytopeTransformation.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Polytope isoparametric transformation.
   */
  template <class FE>
  class IsoparametricTransformation final : public PolytopeTransformation
  {
    static_assert(std::is_same_v<typename FE::RangeType, Real>,
        "Type of finite element must be scalar valued.");

    friend class boost::serialization::access;

    public:
      using Parent = PolytopeTransformation;
      using Parent::transform;
      using Parent::jacobian;
      using Parent::inverse;

      IsoparametricTransformation(Math::PointMatrix&& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * pm : sdim x dof
       */
      IsoparametricTransformation(const Math::PointMatrix& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(Math::PointMatrix&& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const Math::PointMatrix& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const IsoparametricTransformation& other)
        : Parent(other),
          m_pm(other.m_pm),
          m_fe(other.m_fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(IsoparametricTransformation&& other)
        : Parent(std::move(other)),
          m_pm(std::move(other.m_pm)),
          m_fe(std::move(other.m_fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      size_t getOrder() const override
      {
        return m_fe.getOrder();
      }

      size_t getJacobianOrder() const override
      {
        return m_fe.getOrder();
      }

      void transform(const Math::SpatialVector<Real>& rc, Math::SpatialVector<Real>& pc) const override
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

      void jacobian(const Math::SpatialVector<Real>& rc, Math::SpatialMatrix<Real>& res) const override
      {
        const size_t rdim = getReferenceDimension();
        assert(rc.size() >= 0);
        assert(static_cast<size_t>(rc.size()) == rdim);
        const size_t pdim = getPhysicalDimension();
        res.resize(pdim, rdim);
        res.setZero();
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          const auto basis = m_fe.getBasis(local);
          for (size_t i = 0; i < rdim; i++)
          {
            const auto derivative = basis.template getDerivative<1>(i);
            assert(res.col(i).size() == m_pm.col(local).size());
            res.col(i).noalias() += m_pm.col(local) * derivative(rc);
          }
        }
      }

      const Math::PointMatrix& getPointMatrix() const
      {
        return m_pm;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<PolytopeTransformation>(*this);
        ar & m_pm;
        ar & m_fe;
      }

      IsoparametricTransformation* copy() const noexcept override
      {
        return new IsoparametricTransformation(*this);
      }

    private:
      Math::PointMatrix m_pm;
      FE m_fe;
  };
}

#endif
