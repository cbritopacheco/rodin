/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_PARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_PARAMETRICTRANSFORMATION_H

/**
 * @file
 * @brief Parametric transformation for finite-element mesh geometry.
 */

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PointCloud.h"

#include "PolytopeTransformation.h"

#include "ForwardDecls.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Geometry
{
  /**
   * @brief Parametric finite-element transformation for polytopes.
   *
   * Given geometry-node coordinates @f$P_i@f$ and scalar finite-element basis
   * functions @f$\phi_i@f$, the physical coordinate map is
   * @f[
   *   x(r) = \sum_i P_i \phi_i(r).
   * @f]
   *
   * The same class covers affine P1 geometry and curved Pk geometry.
   */
  template <class FE>
  class ParametricTransformation final : public PolytopeTransformation
  {
      static_assert(std::is_same_v<typename FE::RangeType, Real>,
        "Type of finite element must be scalar valued.");

      friend class boost::serialization::access;

    public:
      /// @brief Parent class type.
      using Parent = PolytopeTransformation;
      using Parent::transform;
      using Parent::jacobian;
      using Parent::inverse;

      /// @brief Constructs the transformation from a control-point cloud and
      /// a finite element (the control-point count must match the element's
      /// node count).
      ParametricTransformation(Geometry::PointCloud&& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /// @brief Constructs the transformation from a control-point cloud and
      /// a finite element.
      ParametricTransformation(const Geometry::PointCloud& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /// @brief Constructs the transformation from a control-point cloud and
      /// a finite element.
      ParametricTransformation(Geometry::PointCloud&& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /// @brief Constructs the transformation from a control-point cloud and
      /// a finite element.
      ParametricTransformation(const PointCloud& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /// @brief Copy constructor.
      ParametricTransformation(const ParametricTransformation& other)
        : Parent(other),
          m_pm(other.m_pm),
          m_fe(other.m_fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /// @brief Move constructor.
      ParametricTransformation(ParametricTransformation&& other)
        : Parent(std::move(other)),
          m_pm(std::move(other.m_pm)),
          m_fe(std::move(other.m_fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      size_t getOrder() const override
      {
        return m_fe.getOrder();
      }

      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const override
      {
        const size_t pdim = getPhysicalDimension();
        assert(static_cast<size_t>(rc.size()) == getReferenceDimension());
        pc.resize(pdim);
        pc.setZero();
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          assert(pc.size() == m_pm[local].size());
          pc += m_pm[local] * m_fe.getBasis(local)(rc);
        }
      }

      void jacobian(
        Math::SpatialMatrix<Real>& pc, const Math::SpatialPoint& rc) const override
      {
        const size_t rdim = getReferenceDimension();
        assert(static_cast<size_t>(rc.size()) == rdim);
        const size_t pdim = getPhysicalDimension();
        pc.resize(pdim, rdim);
        pc.setZero();
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          const auto& basis = m_fe.getBasis(local);
          for (size_t i = 0; i < rdim; i++)
          {
            const auto derivative = basis.template getDerivative<1>(i);
            for (size_t j = 0; j < pdim; j++)
              pc(j, i) += m_pm(j, local) * derivative(rc);
          }
        }
      }

      /// @brief Gets the control-point coordinate matrix.
      const PointCloud& getPointMatrix() const
      {
        return m_pm;
      }

      /// @brief Serializes the transformation (for boost::serialization).
      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar& boost::serialization::base_object<PolytopeTransformation>(*this);
        ar & m_pm;
        ar & m_fe;
      }

      ParametricTransformation* copy() const noexcept override
      {
        return new ParametricTransformation(*this);
      }

    private:
      PointCloud m_pm;
      FE m_fe;
  };
}

#endif
