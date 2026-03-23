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
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H
#define RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
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
          m_linData(nullptr)
      {
        static_assert(std::is_same_v<TestFES, FES>);
      }

      InternalForce(const InternalForce& other)
        : Parent(other),
          m_law(other.m_law),
          m_fes(other.m_fes),
          m_linData(other.m_linData)
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

      InternalForce& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linData);
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fes = m_fes.get();
        const size_t vdim = fes.getVectorDimension();

        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, m_qf->getPoint(0));
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();

        const auto& rc = m_qf->getPoint(0);
        const auto& JacInv = m_p->getJacobianInverse();

        // Compute physical gradients of scalar P1 basis functions
        const Variational::P1Element<ScalarType> fe_scalar(polytope.getGeometry());
        const size_t nv = fe_scalar.getCount();

        m_physGrads.resize(nv);
        for (size_t v = 0; v < nv; ++v)
        {
          Math::SpatialVector<ScalarType> ghat(d);
          for (size_t k = 0; k < d; ++k)
            ghat(k) = fe_scalar.getBasis(v).template getDerivative<1>(k)(rc);
          m_physGrads[v] = JacInv.transpose() * ghat;
        }

        // Evaluate displacement gradient H from DOF values
        Math::Matrix<ScalarType> H(vdim, d);
        H.setZero();
        for (size_t v = 0; v < nv; ++v)
          for (size_t c = 0; c < vdim; ++c)
          {
            const size_t local = v * vdim + c;
            const ScalarType uc = (*m_linData)(fes.getGlobalIndex({d, idx}, local));
            for (size_t col = 0; col < d; ++col)
              H(c, col) += uc * m_physGrads[v](col);
          }

        // Compute kinematic state, cache, and stress
        KinematicState state(d);
        state.setDisplacementGradient(H).update();

        typename LawType::Cache cache;
        m_law.setCache(cache, state);

        Math::Matrix<ScalarType> P;
        m_law.getFirstPiolaKirchhoffStress(P, cache, state);

        // Precompute element vector
        const size_t ndof = nv * vdim;
        m_elemVec.resize(ndof);
        for (size_t te = 0; te < ndof; ++te)
        {
          const size_t node_te = te / vdim;
          const size_t comp_te = te % vdim;
          ScalarType val = 0;
          for (size_t k = 0; k < d; ++k)
            val += P(comp_te, k) * m_physGrads[node_te](k);
          m_elemVec(te) = val;
        }

        return *this;
      }

      ScalarType integrate(size_t te) final override
      {
        return m_weight * m_distortion * m_elemVec(te);
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

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      ScalarType m_weight;
      ScalarType m_distortion;
      std::vector<Math::SpatialVector<ScalarType>> m_physGrads;
      Math::Vector<ScalarType> m_elemVec;
  };
}

#endif
