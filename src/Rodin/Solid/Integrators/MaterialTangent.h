/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MaterialTangent.h
 * @brief Material tangent stiffness integrator for Newton-Raphson linearization.
 *
 * Evaluates the bilinear form arising from the linearization of the
 * internal virtual work:
 * @f[
 *   a(\delta\mathbf{u}, \mathbf{v})
 *     = \int_\Omega D\mathbf{P}[\nabla \delta\mathbf{u}]
 *       : \nabla \mathbf{v} \, dX
 * @f]
 * where @f$ D\mathbf{P}[\cdot] @f$ denotes the directional derivative
 * of the first Piola-Kirchhoff stress.
 *
 * The integrator is generic: it supports arbitrary finite element spaces
 * and quadrature rules (not limited to P1/centroid).  The constitutive law
 * receives a ConstitutivePoint at each quadrature point.
 */
#ifndef RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H
#define RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/TrialFunction.h"
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
   * @brief Local bilinear form integrator for the material tangent stiffness
   * in hyperelastic problems.
   *
   * Computes the element-level tangent stiffness matrix:
   * @f[
   *   K^e_{ij} = \int_{K} D\mathbf{P}[\nabla \phi_j] : \nabla \phi_i \, dX
   * @f]
   *
   * Uses generic quadrature (order 2 by default) and builds a
   * ConstitutivePoint at each quadrature point for constitutive evaluation.
   *
   * @tparam LawDerived The hyperelastic constitutive law type
   * @tparam Solution The solution type of the trial function
   * @tparam FES The finite element space type
   */
  template <class LawDerived, class Solution, class FES>
  class MaterialTangent final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using FESType = FES;

      /**
       * @brief Constructs the material tangent integrator.
       * @param law The constitutive law (stored by value)
       * @param u The trial function
       * @param v The test function
       */
      template <class S, class TrialFES, class TestFES>
      MaterialTangent(
          const LawDerived& law,
          const Variational::TrialFunction<S, TrialFES>& u,
          const Variational::TestFunction<TestFES>& v)
        : Parent(u, v),
          m_law(law),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace()),
          m_linData(nullptr),
          m_quadOrder(2)
      {
        static_assert(std::is_same_v<TrialFES, FES>);
        static_assert(std::is_same_v<TestFES, FES>);
      }

      MaterialTangent(const MaterialTangent& other)
        : Parent(other),
          m_law(other.m_law),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_linData(other.m_linData),
          m_quadOrder(other.m_quadOrder)
      {}

      /**
       * @brief Sets the linearization point (current displacement DOF vector).
       * @param data Reference to the coefficient vector of the displacement GridFunction
       * @returns Reference to this object for chaining
       */
      MaterialTangent& setLinearizationPoint(const Math::Vector<ScalarType>& data)
      {
        m_linData = &data;
        return *this;
      }

      /**
       * @brief Sets the quadrature order.
       * @param order Polynomial order for exact integration
       * @returns Reference to this object for chaining
       */
      MaterialTangent& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      MaterialTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linData);
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const auto d_u8 = static_cast<std::uint8_t>(d);
        const Index idx = polytope.getIndex();
        const auto& fes = m_trialfes.get();
        const size_t vdim = fes.getVectorDimension();

        const auto& qf = QF::GenericPolytopeQuadrature::get(m_quadOrder, polytope.getGeometry());
        const size_t nqp = qf.getSize();

        // Get scalar FE for basis gradients (P1 for now, extensible)
        const Variational::P1Element<ScalarType> fe_scalar(polytope.getGeometry());
        const size_t nv = fe_scalar.getCount();
        const size_t ndof = nv * vdim;

        // Zero element stiffness matrix
        m_matrix.resize(ndof, ndof);
        m_matrix.setZero();

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
            Math::SpatialVector<ScalarType> ghat(d_u8);
            for (size_t k = 0; k < d; ++k)
              ghat(k) = fe_scalar.getBasis(v).template getDerivative<1>(k)(rc);
            physGrads[v] = JacInv.transpose() * ghat;
          }

          // Evaluate displacement gradient H from DOF values
          Math::SpatialMatrix<ScalarType> H;
          H.resize(static_cast<std::uint8_t>(vdim), d_u8);
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
          Math::SpatialVector<ScalarType> xiVec(d_u8);
          Math::SpatialVector<ScalarType> xVec(d_u8);
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

          // Build element stiffness at this quadrature point
          for (size_t tr = 0; tr < ndof; ++tr)
          {
            const size_t node_tr = tr / vdim;
            const size_t comp_tr = tr % vdim;

            // Construct dF = e_{comp_tr} ⊗ phys_grad[node_tr]
            Math::SpatialMatrix<ScalarType> dF;
            dF.resize(static_cast<std::uint8_t>(vdim), d_u8);
            dF.setZero();
            for (size_t k = 0; k < d; ++k)
              dF(comp_tr, k) = physGrads[node_tr](k);

            // Compute material tangent action dP = DP[dF]
            Math::SpatialMatrix<ScalarType> dP;
            m_law.getMaterialTangent(dP, cache, cp, dF);

            for (size_t te = 0; te < ndof; ++te)
            {
              const size_t node_te = te / vdim;
              const size_t comp_te = te % vdim;
              ScalarType val = 0;
              for (size_t k = 0; k < d; ++k)
                val += dP(comp_te, k) * physGrads[node_te](k);
              m_matrix(te, tr) += wq * distortion * val;
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      MaterialTangent* copy() const noexcept final override
      {
        return new MaterialTangent(*this);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      LawType m_law;
      std::reference_wrapper<const FESType> m_trialfes;
      std::reference_wrapper<const FESType> m_testfes;
      const Math::Vector<ScalarType>* m_linData;
      size_t m_quadOrder;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  /// CTAD deduction guide for MaterialTangent
  template <class LawDerived, class S, class TrialFES, class TestFES>
  MaterialTangent(const LawDerived&,
                  const Variational::TrialFunction<S, TrialFES>&,
                  const Variational::TestFunction<TestFES>&)
    -> MaterialTangent<LawDerived, S, TrialFES>;
}

#endif
