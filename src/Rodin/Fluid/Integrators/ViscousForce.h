/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ViscousForce.h
 * @brief Internal viscous force vector integrator for fluid formulations.
 *
 * Assembles the nonlinear residual vector (viscous force) contribution:
 * @f[
 *   R_{\text{visc}}(\mathbf{v})
 *     = \int_\Omega \boldsymbol{\tau}(\mathbf{u}) : \mathbf{D}(\mathbf{v}) \, dx
 * @f]
 * where @f$\boldsymbol{\tau}@f$ is the deviatoric Cauchy stress computed by
 * the rheology law and
 * @f$\mathbf{D}(\mathbf{v}) = \tfrac12(\nabla\mathbf{v}+\nabla\mathbf{v}^T)@f$
 * is the symmetric part of the velocity gradient.
 *
 * Since @f$\boldsymbol{\tau}@f$ is symmetric,
 * @f$\boldsymbol{\tau}:\mathbf{D}(\mathbf{v}) = \boldsymbol{\tau}:\nabla\mathbf{v}@f$,
 * so the implementation contracts @f$\boldsymbol\tau@f$ against the full
 * physical Jacobian of the test basis.
 */
#ifndef RODIN_FLUID_INTEGRATORS_VISCOUSFORCE_H
#define RODIN_FLUID_INTEGRATORS_VISCOUSFORCE_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Fluid/Local/FlowPoint.h"
#include "Rodin/Fluid/Constitutive/RheologyLaw.h"

namespace Rodin::Fluid
{
  /**
   * @brief Linear form integrator for the viscous force in generalized
   * Newtonian fluid problems.
   *
   * Computes the element-level contribution to the viscous residual:
   * @f[
   *   R^e_i = \int_{K} \boldsymbol{\tau} : \nabla \phi_i \, dx
   * @f]
   *
   * Obtains the finite element basis from the FE space via
   * @c getFiniteElement(), supports configurable quadrature order, and
   * builds a FlowPoint (composed over @c Geometry::Point) at each
   * quadrature point for constitutive evaluation.
   *
   * @tparam LawDerived The rheology law type (must satisfy the RheologyLaw CRTP contract)
   * @tparam FES The finite element space type
   */
  template <class LawDerived, class FES>
  class ViscousForce final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using FESType = FES;
      using GridFunctionType = Variational::GridFunction<FESType, Math::Vector<ScalarType>>;

      /**
       * @brief Constructs the viscous force integrator.
       * @param law  The rheology law (stored by value)
       * @param v    The test function
       */
      template <class TestFES>
      ViscousForce(const LawDerived& law, const Variational::TestFunction<TestFES>& v)
        : Parent(v),
          m_law(law),
          m_fes(v.getFiniteElementSpace()),
          m_linGf(nullptr),
          m_quadOrder(0)
      {
        static_assert(std::is_same_v<TestFES, FES>);
      }

      ViscousForce(const ViscousForce& other)
        : Parent(other),
          m_law(other.m_law),
          m_fes(other.m_fes),
          m_linGf(other.m_linGf),
          m_quadOrder(other.m_quadOrder)
      {}

      /**
       * @brief Sets the linearization point (current velocity @c GridFunction).
       * @param gf Reference to the velocity GridFunction
       * @returns Reference to this object for chaining
       */
      ViscousForce& setVelocity(const GridFunctionType& gf)
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
      ViscousForce& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      ViscousForce& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linGf);
        m_polytope = polytope;

        const auto& linData = m_linGf->getData();

        const size_t d = polytope.getDimension();
        const auto d_u8 = static_cast<std::uint8_t>(d);
        const Index idx = polytope.getIndex();
        const auto& fes = m_fes.get();
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

        m_elemVec.resize(ndof);
        m_elemVec.setZero();

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const auto& rc = qf.getPoint(q);
          const ScalarType wq = qf.getWeight(q);
          const ScalarType distortion = pt.getDistortion();
          const auto& JacInv = pt.getJacobianInverse();

          // Precompute physical Jacobians of all DOF basis functions.
          // For vector P1: physJacs[dof] is a vdim×d matrix with only row
          // (dof % vdim) nonzero, equal to the physical gradient of the scalar
          // nodal basis function indexed by (dof / vdim).
          std::vector<Math::SpatialMatrix<ScalarType>> physJacs(ndof);
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            Math::SpatialMatrix<ScalarType> refJac = fe.getBasis(dof).getJacobian()(rc);
            physJacs[dof] = refJac * JacInv;
          }

          // Interpolate velocity gradient at this quadrature point.
          // Each DOF carries a scalar coefficient for a single vector component;
          // physJacs[dof] is zero in all rows except (dof % vdim), so the loop
          // correctly accumulates only the contribution of the relevant component.
          // The velocity value is not required by current rheology laws
          // (all implemented laws depend only on the velocity gradient), so a
          // zero placeholder is passed to FlowPoint.
          Math::SpatialVector<ScalarType> velZero(vdim);
          velZero.setZero();

          Math::SpatialMatrix<ScalarType> gradU(vdim, d_u8);
          gradU.setZero();
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            const ScalarType u_dof = linData(fes.getGlobalIndex({d, idx}, dof));
            for (size_t c = 0; c < vdim; ++c)
              for (size_t k = 0; k < d; ++k)
                gradU(c, k) += u_dof * physJacs[dof](c, k);
          }

          // Build FlowPoint at this quadrature point.
          FlowPoint fp(pt, velZero, gradU);

          // Evaluate the deviatoric stress tau from the law.
          typename LawType::Cache cache;
          m_law.setCache(cache, fp);

          Math::SpatialMatrix<ScalarType> tau(vdim, d_u8);
          m_law.getDeviatoricStress(tau, cache, fp);

          // Accumulate: R_te += wq * distortion * (tau : physJacs[te])
          for (size_t te = 0; te < ndof; ++te)
          {
            ScalarType val = 0;
            for (size_t c = 0; c < vdim; ++c)
              for (size_t k = 0; k < d; ++k)
                val += tau(c, k) * physJacs[te](c, k);
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

      ViscousForce* copy() const noexcept final override
      {
        return new ViscousForce(*this);
      }

      const LawType& getLaw() const { return m_law; }

    private:
      LawType m_law;
      std::reference_wrapper<const FESType> m_fes;
      const GridFunctionType* m_linGf;
      size_t m_quadOrder;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /// CTAD deduction guide for ViscousForce
  template <class LawDerived, class TestFES>
  ViscousForce(const LawDerived&, const Variational::TestFunction<TestFES>&)
    -> ViscousForce<LawDerived, TestFES>;
}

#endif
