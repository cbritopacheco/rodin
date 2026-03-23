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
 */
#ifndef RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H
#define RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H

#include <cassert>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
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
          m_linData(nullptr)
      {
        static_assert(std::is_same_v<TrialFES, FES>);
        static_assert(std::is_same_v<TestFES, FES>);
      }

      MaterialTangent(const MaterialTangent& other)
        : Parent(other),
          m_law(other.m_law),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_linData(other.m_linData)
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

      MaterialTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        assert(m_linData);
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fes = m_trialfes.get();
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

        // Compute kinematic state and cache
        KinematicState state(d);
        state.setDisplacementGradient(H).update();

        typename LawType::Cache cache;
        m_law.setCache(cache, state);

        // Precompute element stiffness matrix
        const size_t ndof = nv * vdim;
        m_matrix.resize(ndof, ndof);

        for (size_t tr = 0; tr < ndof; ++tr)
        {
          const size_t node_tr = tr / vdim;
          const size_t comp_tr = tr % vdim;

          // Construct dF = e_{comp_tr} ⊗ phys_grad[node_tr]
          Math::Matrix<ScalarType> dF(vdim, d);
          dF.setZero();
          for (size_t k = 0; k < d; ++k)
            dF(comp_tr, k) = m_physGrads[node_tr](k);

          // Compute material tangent action dP = DP[dF]
          Math::Matrix<ScalarType> dP;
          m_law.getMaterialTangent(dP, cache, state, dF);

          for (size_t te = 0; te < ndof; ++te)
          {
            const size_t node_te = te / vdim;
            const size_t comp_te = te % vdim;
            ScalarType val = 0;
            for (size_t k = 0; k < d; ++k)
              val += dP(comp_te, k) * m_physGrads[node_te](k);
            m_matrix(te, tr) = val;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
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

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Optional<QF::Centroid> m_qf;
      Optional<Geometry::Point> m_p;

      ScalarType m_weight;
      ScalarType m_distortion;
      std::vector<Math::SpatialVector<ScalarType>> m_physGrads;
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
