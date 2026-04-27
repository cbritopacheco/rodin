/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ViscousTangent.h
 * @brief Viscous tangent stiffness integrator for Newton linearization of
 * generalized Newtonian fluid problems.
 *
 * Evaluates the bilinear form arising from the linearization of the viscous
 * virtual work:
 * @f[
 *   a(\delta\mathbf{u}, \mathbf{v})
 *     = \int_\Omega \mathrm{D}\boldsymbol{\tau}[\nabla\delta\mathbf{u}]
 *       : \mathbf{D}(\mathbf{v}) \, dx
 * @f]
 * where @f$ \mathrm{D}\boldsymbol{\tau}[\cdot] @f$ denotes the directional
 * derivative of the deviatoric Cauchy stress.  Since
 * @f$\mathrm{D}\boldsymbol{\tau}@f$ is symmetric, the contraction against
 * @f$\mathbf{D}(\mathbf{v})@f$ equals the contraction against the full
 * Jacobian @f$\nabla\mathbf{v}@f$, so the implementation contracts the
 * tangent action against the physical Jacobian of the test basis.
 */
#ifndef RODIN_FLUID_INTEGRATORS_VISCOUSTANGENT_H
#define RODIN_FLUID_INTEGRATORS_VISCOUSTANGENT_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Fluid/Local/FlowPoint.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Local bilinear form integrator for the viscous tangent stiffness
   * in generalized Newtonian fluid problems.
   *
   * Computes the element-level tangent stiffness matrix:
   * @f[
   *   K^e_{ij}
   *     = \int_{K} \mathrm{D}\boldsymbol{\tau}[\nabla\phi_j]
   *       : \nabla\phi_i \, dx
   * @f]
   *
   * Obtains the finite element basis from the FE space via
   * @c getFiniteElement(), supports configurable quadrature order, and
   * builds a FlowPoint (composed over @c Geometry::Point) at each quadrature
   * point for constitutive evaluation.
   *
   * @tparam LawDerived The rheology law type (must satisfy the RheologyLaw CRTP contract)
   * @tparam Solution The solution type of the trial function
   * @tparam FES The finite element space type
   */
  template <class LawDerived, class Solution, class FES>
  class ViscousTangent final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using FESType = FES;
      using GridFunctionType = Variational::GridFunction<FESType, Math::Vector<ScalarType>>;

      /**
       * @brief Constructs the viscous tangent integrator.
       * @param law  The rheology law (stored by value)
       * @param u    The trial function (velocity increment @f$\delta\mathbf{u}@f$)
       * @param v    The test function
       */
      template <class S, class TrialFES, class TestFES>
      ViscousTangent(
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

      ViscousTangent(const ViscousTangent& other)
        : Parent(other),
          m_law(other.m_law),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_linGf(other.m_linGf),
          m_quadOrder(other.m_quadOrder)
      {}

      /**
       * @brief Sets the linearization point (current velocity @c GridFunction).
       *
       * The rheology cache is precomputed at this velocity field at each
       * quadrature point, and the tangent action is evaluated relative to it.
       * For linear laws (e.g., Newtonian) this point does not affect the
       * resulting stiffness matrix.
       *
       * @param gf Reference to the velocity GridFunction
       * @returns Reference to this object for chaining
       */
      ViscousTangent& setVelocity(const GridFunctionType& gf)
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
      ViscousTangent& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      ViscousTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const auto d_u8 = static_cast<std::uint8_t>(d);
        const Index idx = polytope.getIndex();
        const auto& fes = m_trialfes.get();
        const size_t vdim = fes.getVectorDimension();

        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t ndof = fe.getCount();

        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * fe.getOrder();
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        m_matrix.resize(ndof, ndof);
        m_matrix.setZero();

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const auto& rc = qf.getPoint(q);
          const ScalarType wq = qf.getWeight(q);
          const ScalarType distortion = pt.getDistortion();
          const auto& JacInv = pt.getJacobianInverse();

          // Precompute physical Jacobians of all DOF basis functions.
          std::vector<Math::SpatialMatrix<ScalarType>> physJacs(ndof);
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            Math::SpatialMatrix<ScalarType> refJac = fe.getBasis(dof).getJacobian()(rc);
            physJacs[dof] = refJac * JacInv;
          }

          // Interpolate velocity and velocity gradient at the linearization
          // point (from m_linGf if set, zero otherwise).
          Math::SpatialVector<ScalarType> vel(vdim);
          vel.setZero();
          Math::SpatialMatrix<ScalarType> gradU(vdim, d_u8);
          gradU.setZero();

          if (m_linGf)
          {
            const auto& linData = m_linGf->getData();
            for (size_t dof = 0; dof < ndof; ++dof)
            {
              const ScalarType u_dof = linData(fes.getGlobalIndex({d, idx}, dof));
              for (size_t c = 0; c < vdim; ++c)
                for (size_t k = 0; k < d; ++k)
                  gradU(c, k) += u_dof * physJacs[dof](c, k);
            }
          }

          // Build FlowPoint at the linearization point.
          FlowPoint fp(pt, vel, gradU);

          // Precompute constitutive cache at this linearization point.
          typename LawType::Cache cache;
          m_law.setCache(cache, fp);

          // Build element stiffness at this quadrature point.
          // For trial DOF tr, the velocity gradient perturbation is physJacs[tr].
          for (size_t tr = 0; tr < ndof; ++tr)
          {
            const auto& dGradU = physJacs[tr];

            // Compute tangent action: dtau = D_tau[dGradU]
            Math::SpatialMatrix<ScalarType> dtau(vdim, d_u8);
            m_law.getTangent(dtau, cache, fp, dGradU);

            // K_{te,tr} += wq * distortion * (dtau : physJacs[te])
            for (size_t te = 0; te < ndof; ++te)
            {
              ScalarType val = 0;
              for (size_t c = 0; c < vdim; ++c)
                for (size_t k = 0; k < d; ++k)
                  val += dtau(c, k) * physJacs[te](c, k);
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

      ViscousTangent* copy() const noexcept final override
      {
        return new ViscousTangent(*this);
      }

      const LawType& getLaw() const { return m_law; }

    private:
      LawType m_law;
      std::reference_wrapper<const FESType> m_trialfes;
      std::reference_wrapper<const FESType> m_testfes;
      const GridFunctionType* m_linGf;
      size_t m_quadOrder;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  /// CTAD deduction guide for ViscousTangent
  template <class LawDerived, class S, class TrialFES, class TestFES>
  ViscousTangent(
      const LawDerived&,
      const Variational::TrialFunction<S, TrialFES>&,
      const Variational::TestFunction<TestFES>&)
    -> ViscousTangent<LawDerived, S, TrialFES>;
}

#endif
