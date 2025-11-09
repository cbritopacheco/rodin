/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

/**
 * @file
 * @brief Isoparametric transformation for finite elements.
 */

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "PolytopeTransformation.h"

#include "ForwardDecls.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Geometry
{
  /**
   * @brief Isoparametric transformation for polytopes.
   *
   * An isoparametric transformation uses finite element basis functions to
   * map reference coordinates to physical coordinates. Given a point matrix
   * @f$ P @f$ containing the physical coordinates of degrees of freedom and
   * a finite element @f$ \text{FE} @f$ with basis functions @f$ \phi_i @f$,
   * the transformation is:
   * @f[
   *    x(r) = \sum_{i=1}^{n} P_i \phi_i(r)
   * @f]
   * where @f$ P_i @f$ is the @f$ i @f$-th column of the point matrix and
   * @f$ n @f$ is the number of degrees of freedom.
   *
   * The Jacobian matrix is computed as:
   * @f[
   *    \mathbf{J}_x(r) = \sum_{i=1}^{n} P_i \otimes \nabla \phi_i(r)
   * @f]
   *
   * This is the standard transformation used for curved finite elements and
   * enables higher-order geometric approximations.
   *
   * @tparam FE Finite element type (must have scalar-valued basis functions)
   *
   * @note The finite element type must provide scalar-valued basis functions
   * (FE::RangeType must be Real).
   *
   * @see PolytopeTransformation, IdentityTransformation
   */
  template <class FE>
  class IsoparametricTransformation final : public PolytopeTransformation
  {
    static_assert(std::is_same_v<typename FE::RangeType, Real>,
        "Type of finite element must be scalar valued.");

    friend class boost::serialization::access;

    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeTransformation;
      using Parent::transform;
      using Parent::jacobian;
      using Parent::inverse;

      /**
       * @brief Constructs an isoparametric transformation (move semantics).
       * @param[in] pm Point matrix of size @f$ s \times n @f$ where @f$ s @f$
       *            is the spatial dimension and @f$ n @f$ is the number of DOFs
       * @param[in] fe Finite element providing basis functions
       */
      IsoparametricTransformation(Math::PointMatrix&& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs an isoparametric transformation (copy semantics).
       * @param[in] pm Point matrix of size @f$ s \times n @f$
       * @param[in] fe Finite element providing basis functions
       */
      IsoparametricTransformation(const Math::PointMatrix& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs an isoparametric transformation (mixed semantics).
       * @param[in] pm Point matrix (move)
       * @param[in] fe Finite element (copy)
       */
      IsoparametricTransformation(Math::PointMatrix&& pm, const FE& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(std::move(pm)),
          m_fe(fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Constructs an isoparametric transformation (mixed semantics).
       * @param[in] pm Point matrix (copy)
       * @param[in] fe Finite element (move)
       */
      IsoparametricTransformation(const Math::PointMatrix& pm, FE&& fe)
        : Parent(Polytope::Traits(fe.getGeometry()).getDimension(), pm.rows()),
          m_pm(pm),
          m_fe(std::move(fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Copy constructor.
       */
      IsoparametricTransformation(const IsoparametricTransformation& other)
        : Parent(other),
          m_pm(other.m_pm),
          m_fe(other.m_fe)
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Move constructor.
       */
      IsoparametricTransformation(IsoparametricTransformation&& other)
        : Parent(std::move(other)),
          m_pm(std::move(other.m_pm)),
          m_fe(std::move(other.m_fe))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * @brief Gets the polynomial order of the transformation.
       * @returns Order of the finite element basis functions
       */
      size_t getOrder() const override
      {
        return m_fe.getOrder();
      }

      /**
       * @brief Gets the polynomial order of the Jacobian.
       * @returns Order of the Jacobian (same as basis function order)
       */
      size_t getJacobianOrder() const override
      {
        return m_fe.getOrder();
      }

      /**
       * @brief Applies the isoparametric transformation.
       * @param[out] pc Physical coordinates
       * @param[in] rc Reference coordinates
       *
       * Computes @f$ pc = x(rc) = \sum_{i=1}^{n} P_i \phi_i(rc) @f$.
       */
      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const override
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

      /**
       * @brief Computes the Jacobian matrix.
       * @param[out] pc Jacobian matrix of size @f$ s \times d @f$
       * @param[in] rc Reference coordinates
       *
       * Computes @f$ \mathbf{J}_x(rc) = \sum_{i=1}^{n} P_i \otimes \nabla \phi_i(rc) @f$.
       */
      void jacobian(Math::SpatialMatrix<Real>& pc, const Math::SpatialPoint& rc) const override
      {
        const size_t rdim = getReferenceDimension();
        assert(rc.size() >= 0);
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
            assert(pc.col(i).size() == m_pm.col(local).size());
            pc.col(i).noalias() += m_pm.col(local) * derivative(rc);
          }
        }
      }

      /**
       * @brief Gets the point matrix.
       * @returns Reference to the point matrix containing DOF coordinates
       *
       * The returned matrix has size @f$ s \times n @f$ where @f$ s @f$ is
       * the spatial dimension and @f$ n @f$ is the number of degrees of freedom.
       */
      const Math::PointMatrix& getPointMatrix() const
      {
        return m_pm;
      }

      /**
       * @brief Serialization method for Boost.Serialization.
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version (unused)
       */
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<PolytopeTransformation>(*this);
        ar & m_pm;
        ar & m_fe;
      }

      /**
       * @brief Creates a copy of this transformation.
       * @returns Pointer to a new IsoparametricTransformation object
       */
      IsoparametricTransformation* copy() const noexcept override
      {
        return new IsoparametricTransformation(*this);
      }

    private:
      Math::PointMatrix m_pm; ///< Point matrix (spatial_dim x num_dofs)
      FE m_fe;                ///< Finite element providing basis functions
  };
}

#endif
