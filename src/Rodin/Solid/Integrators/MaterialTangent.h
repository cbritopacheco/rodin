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
 * The integrator is generic: it obtains the finite element basis from the
 * FE space (not hardcoded to P1), supports arbitrary quadrature rules, and
 * builds a ConstitutivePoint (composed over Geometry::Point) at each
 * quadrature point for constitutive evaluation.
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
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"
#include "Rodin/Solid/Local/Input.h"
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
   * Obtains the finite element basis from the FE space via
   * @c getFiniteElement(), supports configurable quadrature order, and
   * builds a ConstitutivePoint (composed over Geometry::Point) at each
   * quadrature point for constitutive evaluation.  An optional
   * Input can inject auxiliary data (fiber directions, activation,
   * etc.) into the ConstitutivePoint at each quadrature point.
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
      using GridFunctionType = Variational::GridFunction<FESType, Math::Vector<ScalarType>>;

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
          m_linGf(nullptr),
          m_quadOrder(0)
      {
        static_assert(std::is_same_v<TrialFES, FES>);
        static_assert(std::is_same_v<TestFES, FES>);
      }

      MaterialTangent(const MaterialTangent& other)
        : Parent(other),
          m_law(other.m_law),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_linGf(other.m_linGf),
          m_quadOrder(other.m_quadOrder),
          m_input(other.m_input)
      {}

      /**
       * @brief Sets the linearization point (current displacement GridFunction).
       * @param gf Reference to the displacement GridFunction
       * @returns Reference to this object for chaining
       */
      MaterialTangent& setLinearizationPoint(const GridFunctionType& gf)
      {
        m_linGf = &gf;
        return *this;
      }

      /**
       * @brief Sets the quadrature order.
       *
       * When set to 0 (the default), the quadrature order is determined
       * automatically from the finite element approximation order as
       * @c 2 * fe.getOrder().
       *
       * @param order Polynomial order for exact integration (0 = auto)
       * @returns Reference to this object for chaining
       */
      MaterialTangent& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /**
       * @brief Sets an input for auxiliary constitutive data.
       *
       * The input is called at each quadrature point after the
       * ConstitutivePoint has been constructed with geometric context and
       * kinematics, allowing injection of fiber directions, activation
       * parameters, region-wise material properties, etc.
       *
       * @param input A callable with signature void(ConstitutivePoint&)
       * @returns Reference to this object for chaining
       */
      MaterialTangent& setInput(InputFunction input)
      {
        m_input = std::move(input);
        return *this;
      }

      MaterialTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linGf);
        m_polytope = polytope;

        const auto& linData = m_linGf->getData();

        const size_t d = polytope.getDimension();
        const auto d_u8 = static_cast<std::uint8_t>(d);
        const Index idx = polytope.getIndex();
        const auto& fes = m_trialfes.get();
        const size_t vdim = fes.getVectorDimension();

        // Get element from the FE space (not hardcoded to any element type)
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t ndof = fe.getCount();

        // Determine effective quadrature order
        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * fe.getOrder();
        const auto& qf = QF::GenericPolytopeQuadrature::get(effectiveOrder, polytope.getGeometry());
        const size_t nqp = qf.getSize();

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

          // Precompute physical Jacobians for each DOF via the FE element.
          // physJacs[dof] = refJac(dof) * JacInv, where refJac is the
          // vdim × d reference Jacobian of the basis function.
          std::vector<Math::SpatialMatrix<ScalarType>> physJacs(ndof);
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            Math::SpatialMatrix<ScalarType> refJac = fe.getBasis(dof).getJacobian()(rc);
            physJacs[dof] = refJac * JacInv;
          }

          // Evaluate displacement gradient H from DOF values
          Math::SpatialMatrix<ScalarType> H;
          H.resize(static_cast<std::uint8_t>(vdim), d_u8);
          H.setZero();
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            const ScalarType u_dof = linData(fes.getGlobalIndex({d, idx}, dof));
            for (size_t c = 0; c < vdim; ++c)
              for (size_t k = 0; k < d; ++k)
                H(c, k) += u_dof * physJacs[dof](c, k);
          }

          // Build ConstitutivePoint composed over Geometry::Point
          KinematicState state(d);
          state.setDisplacementGradient(H);

          ConstitutivePoint cp(pt, state);

          // Invoke input for auxiliary data injection
          if (m_input)
            m_input(cp);

          typename LawType::Cache cache;
          m_law.setCache(cache, cp);

          // Build element stiffness at this quadrature point.
          // For trial DOF tr, the deformation gradient perturbation dF equals
          // physJacs[tr] (the physical Jacobian of the basis function).
          for (size_t tr = 0; tr < ndof; ++tr)
          {
            const auto& dF = physJacs[tr];

            // Compute material tangent action dP = DP[dF]
            Math::SpatialMatrix<ScalarType> dP;
            m_law.getMaterialTangent(dP, cache, cp, dF);

            // K_{te,tr} += wq * distortion * (dP : physJac_te)
            for (size_t te = 0; te < ndof; ++te)
            {
              ScalarType val = 0;
              for (size_t c = 0; c < vdim; ++c)
                for (size_t k = 0; k < d; ++k)
                  val += dP(c, k) * physJacs[te](c, k);
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
      const GridFunctionType* m_linGf;
      size_t m_quadOrder;
      InputFunction m_input;

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
