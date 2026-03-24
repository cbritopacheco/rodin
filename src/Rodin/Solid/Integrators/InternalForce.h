/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalForce.h
 * @brief Internal force vector integrator for hyperelastic formulations.
 *
 * Assembles the nonlinear residual vector (internal force) contribution:
 * @f[
 *   R_{\text{int}}(\mathbf{v}) = \int_\Omega \mathbf{P}(\mathbf{u}) : \nabla \mathbf{v} \, dX
 * @f]
 * where @f$ \mathbf{P} @f$ is the first Piola-Kirchhoff stress and
 * @f$ \mathbf{v} @f$ is the test function.
 *
 * The integrator is generic: it supports arbitrary finite element spaces
 * and quadrature rules (not limited to P1/centroid).  The constitutive law
 * receives a ConstitutivePoint at each quadrature point, which bundles
 * kinematics, coordinates, and region id.
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H
#define RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Inputs/ConstitutivePoint.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Linear form integrator for the internal force vector in
   * hyperelastic problems.
   *
   * Computes the element-level contribution to the nonlinear residual:
   * @f[
   *   R^e_i = \int_{K} \mathbf{P} : \nabla \phi_i \, dX
   * @f]
   *
   * Uses generic quadrature (order 2 by default) and builds a
   * ConstitutivePoint at each quadrature point for constitutive evaluation.
   *
   * @tparam LawDerived The hyperelastic constitutive law type
   * @tparam FES The finite element space type
   */
  template <class LawDerived, class FES>
  class InternalForce final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using FESType = FES;

      /**
       * @brief Constructs the internal force integrator.
       * @param law The constitutive law (stored by value)
       * @param v The test function
       */
      template <class TestFES>
      InternalForce(const LawDerived& law, const Variational::TestFunction<TestFES>& v)
        : Parent(v),
          m_law(law),
          m_fes(v.getFiniteElementSpace()),
          m_linData(nullptr),
          m_quadOrder(2)
      {
        static_assert(std::is_same_v<TestFES, FES>);
      }

      InternalForce(const InternalForce& other)
        : Parent(other),
          m_law(other.m_law),
          m_fes(other.m_fes),
          m_linData(other.m_linData),
          m_quadOrder(other.m_quadOrder)
      {}

      /**
       * @brief Sets the linearization point (current displacement DOF vector).
       * @param data Reference to the coefficient vector of the displacement GridFunction
       * @returns Reference to this object for chaining
       */
      InternalForce& setLinearizationPoint(const Math::Vector<ScalarType>& data)
      {
        m_linData = &data;
        return *this;
      }

      /**
       * @brief Sets the quadrature order.
       * @param order Polynomial order for exact integration
       * @returns Reference to this object for chaining
       */
      InternalForce& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      InternalForce& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linData);
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fes = m_fes.get();
        const size_t vdim = fes.getVectorDimension();

        const auto& qf = QF::GenericPolytopeQuadrature::get(m_quadOrder, polytope.getGeometry());
        const size_t nqp = qf.getSize();

        // Get scalar FE for basis gradients (P1 for now, extensible)
        const Variational::P1Element<ScalarType> fe_scalar(polytope.getGeometry());
        const size_t nv = fe_scalar.getCount();
        const size_t ndof = nv * vdim;

        // Zero element vector
        m_elemVec.resize(ndof);
        m_elemVec.setZero();

        // Loop over quadrature points
        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& rc = qf.getPoint(q);
          const ScalarType wq = qf.getWeight(q);

          Geometry::Point pt(polytope, rc);
          const ScalarType distortion = pt.getDistortion();
          const auto& JacInv = pt.getJacobianInverse();

          // Physical gradients at this quadrature point
          std::vector<Math::SpatialVector<ScalarType>> physGrads(nv);
          for (size_t v = 0; v < nv; ++v)
          {
            Math::SpatialVector<ScalarType> ghat(static_cast<std::uint8_t>(d));
            for (size_t k = 0; k < d; ++k)
              ghat(k) = fe_scalar.getBasis(v).template getDerivative<1>(k)(rc);
            physGrads[v] = JacInv.transpose() * ghat;
          }

          // Evaluate displacement gradient H from DOF values
          Math::SpatialMatrix<ScalarType> H;
          H.resize(static_cast<std::uint8_t>(vdim), static_cast<std::uint8_t>(d));
          H.setZero();
          for (size_t v = 0; v < nv; ++v)
            for (size_t c = 0; c < vdim; ++c)
            {
              const size_t local = v * vdim + c;
              const ScalarType uc = (*m_linData)(fes.getGlobalIndex({d, idx}, local));
              for (size_t col = 0; col < d; ++col)
                H(c, col) += uc * physGrads[v](col);
            }

          // Build ConstitutivePoint
          KinematicState state(d);
          state.setDisplacementGradient(H);

          ConstitutivePoint cp(state);
          Math::SpatialVector<ScalarType> xiVec(static_cast<std::uint8_t>(d));
          Math::SpatialVector<ScalarType> xVec(static_cast<std::uint8_t>(d));
          for (size_t k = 0; k < d; ++k)
          {
            xiVec(static_cast<std::uint8_t>(k)) = rc(static_cast<std::uint8_t>(k));
            xVec(static_cast<std::uint8_t>(k)) = pt.getPhysicalCoordinates()(static_cast<std::uint8_t>(k));
          }
          cp.setReferenceCoordinates(xiVec);
          cp.setPhysicalCoordinates(xVec);
          if (polytope.getAttribute())
            cp.setRegionId(*polytope.getAttribute());

          typename LawType::Cache cache;
          m_law.setCache(cache, cp);

          Math::SpatialMatrix<ScalarType> P;
          m_law.getFirstPiolaKirchhoffStress(P, cache, cp);

          // Accumulate into element vector
          for (size_t te = 0; te < ndof; ++te)
          {
            const size_t node_te = te / vdim;
            const size_t comp_te = te % vdim;
            ScalarType val = 0;
            for (size_t k = 0; k < d; ++k)
              val += P(comp_te, k) * physGrads[node_te](k);
            m_elemVec(te) += wq * distortion * val;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t te) final override
      {
        return m_elemVec(te);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      InternalForce* copy() const noexcept final override
      {
        return new InternalForce(*this);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      LawType m_law;
      std::reference_wrapper<const FESType> m_fes;
      const Math::Vector<ScalarType>* m_linData;
      size_t m_quadOrder;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /// CTAD deduction guide for InternalForce
  template <class LawDerived, class TestFES>
  InternalForce(const LawDerived&, const Variational::TestFunction<TestFES>&)
    -> InternalForce<LawDerived, TestFES>;
}

#endif
