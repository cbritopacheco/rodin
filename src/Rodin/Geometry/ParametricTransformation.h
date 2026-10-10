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
   *
   * ## Quadrature evaluation
   *
   * A bulk Jacobian evaluation borrows the element's existing reference
   * tabulation for the duration of the call, then contracts it with this
   * transformation's control points. It retains no table reference or new
   * cache. Elements without tabulation use the pointwise interface. The
   * scalar accumulation order is unchanged in either case.
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

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element (the control-point count must match the element's
       * node count).
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(Geometry::PointCloud&& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(const Geometry::PointCloud& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(Geometry::PointCloud&& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs the transformation from a control-point cloud and
       * a finite element.
       * @param pm Point at which the operation is evaluated.
       * @param fe Finite element used by the operation.
       */
      ParametricTransformation(const PointCloud& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      ParametricTransformation(const ParametricTransformation& other)
        : Parent(other),
          m_pm(other.m_pm),
          m_fe(other.m_fe)
      {
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
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

      /**
       * @brief Returns the element's factor degree when available, otherwise
       * the conservative total-degree bound of the transformation interface.
       * @returns The element's factor degree when available, otherwise the conservative total-degree bound of the transformation interface.
       */
      size_t getFactorOrder() const override
      {
        if constexpr (requires { m_fe.getFactorOrder(); })
          return m_fe.getFactorOrder();
        else
          return Parent::getFactorOrder();
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
            // The reference derivative is shared by all physical components.
            const Real value = derivative(rc);
            // Expose the spatial storage bound to keep this a small loop.
            for (size_t j = 0; j < pdim && j < Math::SpatialMatrix<Real>::MaxSize; j++)
              pc(j, i) += m_pm(j, local) * value;
          }
        }
      }

      /**
       * @brief Computes quadrature Jacobians using the element's reference table.
       *
       * For control coordinates @f$P_{ja}@f$, each matrix entry is
       * @f$J_{ji}(\hat x_q)=\sum_a P_{ja}\partial_i\phi_a(\hat x_q)@f$.
       * Only reference derivatives are shared across cells; mapped matrices
       * remain specific to the current transformation. Generic scalar elements
       * without a table retain the default pointwise evaluation.
       *
       * @param[out] jacobians Jacobian matrices in quadrature-point order.
       * @param[in] qf Reference formula identifying the element tabulation.
       */
      void jacobian(std::vector<Math::SpatialMatrix<Real>>& jacobians,
        const QF::QuadratureFormulaBase& qf) const override
      {
        if constexpr (requires { m_fe.getTabulation(qf); })
        {
          const auto& table = m_fe.getTabulation(qf);
          const size_t rdim = getReferenceDimension(), pdim = getPhysicalDimension();
          assert(table.dim == rdim);
          assert(table.ndof == m_fe.getCount());
          jacobians.resize(table.nqp);
          for (size_t qp = 0; qp < table.nqp; ++qp)
          {
            auto& matrix = jacobians[qp];
            matrix.resize(pdim, rdim);
            matrix.setZero();
            for (size_t local = 0; local < m_fe.getCount(); ++local)
            {
              const auto gradient = table.getGradient(qp, local);
              for (size_t i = 0; i < rdim; ++i)
              {
                const Real value = gradient[i];
                for (size_t j = 0; j < pdim && j < Math::SpatialMatrix<Real>::MaxSize;
                     ++j)
                  matrix(j, i) += m_pm(j, local) * value;
              }
            }
          }
        }
        else
          Parent::jacobian(jacobians, qf);
      }

      /**
       * @brief Gets the control-point coordinate matrix.
       * @returns The control-point coordinate matrix.
       */
      const PointCloud& getPointMatrix() const
      {
        return m_pm;
      }

      /**
       * @brief Serializes the transformation (for boost::serialization).
       * @param ar Serialization archive.
       * @param version Boost.Serialization class version; unused by this implementation.
       */
      template <class Archive>
      void serialize(Archive& ar, [[maybe_unused]] const unsigned int version)
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
      // Boost constructs this empty state before loading all serialized fields.
      ParametricTransformation()
        : Parent(0, 0)
      {}

      PointCloud m_pm;
      FE m_fe;
  };
}

#endif
