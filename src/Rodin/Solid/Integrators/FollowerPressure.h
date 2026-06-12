/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FollowerPressure.h
 * @brief Follower (deformation-dependent) pressure load with consistent
 *        tangent, for total-Lagrangian hyperelastic formulations.
 *
 * A pressure @f$ p @f$ acting on the CURRENT configuration of a boundary
 * surface contributes the external virtual work
 * @f[
 *   \delta W_{\text{ext}}
 *     = -p \int_{\Gamma(t)} \mathbf{n} \cdot \mathbf{w} \, da
 *     = -p \int_{\hat{K}}
 *         (\mathbf{x}_{,\xi} \times \mathbf{x}_{,\eta}) \cdot \mathbf{w}
 *         \, d\xi \, d\eta ,
 * @f]
 * where @f$ \mathbf{x}(\xi,\eta) = \mathbf{X}(\xi,\eta) + \mathbf{d} @f$ is
 * the deformed surface and the cross product carries BOTH the outward normal
 * direction and the current area element (Nanson).  In residual convention
 * (R = internal - external), the load term and its EXACT linearization are
 * @f[
 *   R_p(\mathbf{w}) = +p \int_{\hat{K}}
 *       (\mathbf{x}_{,\xi} \times \mathbf{x}_{,\eta}) \cdot \mathbf{w} ,
 *   \qquad
 *   K_p(\delta\mathbf{d}, \mathbf{w}) = +p \int_{\hat{K}}
 *       (\delta\mathbf{d}_{,\xi} \times \mathbf{x}_{,\eta}
 *        + \mathbf{x}_{,\xi} \times \delta\mathbf{d}_{,\eta})
 *       \cdot \mathbf{w} ,
 * @f]
 * (cf. Bonet & Wood, Nonlinear Continuum Mechanics for Finite Element
 * Analysis, ch. 8).  Only SURFACE-parametric derivatives appear: the
 * operators are assembled face-locally, with no volume trace required.  The
 * load-stiffness @f$ K_p @f$ is generally NONSYMMETRIC; use a direct solver
 * or a nonsymmetric Krylov method.
 *
 * At @f$ \mathbf{d} = 0 @f$ the residual reduces EXACTLY to the dead load
 * @f$ +p \int_{\Gamma_0} \mathbf{N} \cdot \mathbf{w} \, dA @f$ (the classic
 * @c BoundaryIntegral(p * Dot(w, N)) term), which both fixes the sign
 * convention (positive p pushes AGAINST the outward normal of the solid,
 * as a pressure must; for a lumen pressure on an inner wall, whose outward
 * normal points INTO the lumen, this pushes the wall outward) and gives
 * a built-in verification limit.
 *
 * Both integrators support the standard @c .over(attribute) restriction.
 * Supported spaces: vector P1 (the local basis is evaluated analytically
 * with the framework's @f$ \text{local} = \text{vertex} \cdot d_v +
 * \text{component} @f$ layout). Supported dimensions: 3D (triangular faces)
 * and 2D (segment faces, where the "cross product" is the +90-degree
 * rotation of the tangent).
 *
 * Usage (static prestress with FULL Newton convergence):
 * @code
 * Real p = 0.0; // ramped externally
 * Solid::FollowerPressureForce   load(p, w, dState);    // residual, "+"
 * Solid::FollowerPressureTangent loadK(p, d, w, dState); // tangent,  "+"
 * problem = tangent + internal
 *         + loadK + load // note: both PLUS, matching internal/tangent
 *         + DirichletBC(...);
 * @endcode
 * The pressure is captured BY REFERENCE (the integrators are copied into the
 * Problem, and all copies read the same variable), so ramping the original
 * variable between solves updates the load with no reassembly bookkeeping.
 */
#ifndef RODIN_SOLID_INTEGRATORS_FOLLOWERPRESSURE_H
#define RODIN_SOLID_INTEGRATORS_FOLLOWERPRESSURE_H

#include <cassert>
#include <functional>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/TrialFunction.h"
#include "Rodin/Variational/TestFunction.h"

namespace Rodin::Solid
{
  namespace Internal
  {
    /// Scalar P1 basis on the reference segment/triangle.
    inline Real p1FaceBasis(size_t a, const Math::SpatialVector<Real>& rc,
                            size_t faceDim)
    {
      if (faceDim == 1)
        return (a == 0) ? (1.0 - rc(0)) : rc(0);
      // triangle
      switch (a)
      {
        case 0:  return 1.0 - rc(0) - rc(1);
        case 1:  return rc(0);
        default: return rc(1);
      }
    }

    /// Constant reference gradient component d(phi_a)/d(xi_j).
    inline Real p1FaceBasisGrad(size_t a, size_t j, size_t faceDim)
    {
      if (faceDim == 1)
        return (a == 0) ? -1.0 : 1.0;
      // triangle
      if (a == 0) return -1.0;
      return ((a == 1 && j == 0) || (a == 2 && j == 1)) ? 1.0 : 0.0;
    }

    /// "Cross product" carrying normal direction and area element:
    /// 3D: t1 x t2; 2D: +90-degree rotation of t1 (t2 ignored).
    inline void surfaceCross(Math::SpatialVector<Real>& out,
                             const Math::SpatialVector<Real>& t1,
                             const Math::SpatialVector<Real>& t2,
                             size_t sdim)
    {
      out.resize(sdim);
      if (sdim == 3)
      {
        out(0) = t1(1) * t2(2) - t1(2) * t2(1);
        out(1) = t1(2) * t2(0) - t1(0) * t2(2);
        out(2) = t1(0) * t2(1) - t1(1) * t2(0);
      }
      else
      {
        out(0) =  t1(1);
        out(1) = -t1(0);
      }
    }

    /// Deformed surface tangents x_{,j} = X_{,j} + sum_a d_a dphi_a/dxi_j and
    /// the outward-orientation sign of surfaceCross(x_{,0}, x_{,1}) computed
    /// from the REFERENCE geometry against the incident cell centroid.
    template <class StateFES, class State>
    struct FaceKinematics
    {
      size_t sdim;
      size_t faceDim;
      size_t nv;
      std::vector<Index> vertices;
      Math::SpatialMatrix<Real> Xjac; ///< Reference geometric Jacobian X_{,xi}
      Real orientation;               ///< +-1: makes the cross point OUTWARD

      void setFace(const Geometry::Polytope& face)
      {
        const auto& mesh = face.getMesh();
        sdim = mesh.getSpaceDimension();
        faceDim = face.getDimension();
        const auto& vs = face.getVertices();
        vertices.assign(vs.begin(), vs.end());
        nv = vertices.size();

        // Reference tangents from the (affine) vertex coordinates:
        // X_{,xi_j} = X_{j+1} - X_0 for the unit simplex parametrization,
        // consistent with p1FaceBasis/p1FaceBasisGrad above.
        Xjac.resize(sdim, faceDim);
        const auto X0 = mesh.getVertexCoordinates(vertices[0]);
        for (size_t j = 0; j < faceDim; ++j)
        {
          const auto Xj = mesh.getVertexCoordinates(vertices[j + 1]);
          for (size_t c = 0; c < sdim; ++c)
            Xjac(c, j) = Xj(c) - X0(c);
        }

        // Orientation: sign such that surfaceCross points away from the
        // incident cell centroid (outward of the solid), evaluated on the
        // REFERENCE configuration (orientation is topological and constant).
        Math::SpatialVector<Real> t1(sdim), t2(sdim), cr;
        for (size_t c = 0; c < sdim; ++c)
        {
          t1(c) = Xjac(c, 0);
          t2(c) = (faceDim == 2) ? Xjac(c, 1) : 0.0;
        }
        surfaceCross(cr, t1, t2, sdim);

        const auto& conn = mesh.getConnectivity();
        const auto& inc = conn.getIncidence({ faceDim, faceDim + 1 },
                                            face.getIndex());
        assert(inc.size() == 1);
        const auto cellIt = mesh.getPolytope(faceDim + 1, *inc.begin());
        const auto rcCell =
          Geometry::Polytope::Traits(cellIt->getGeometry()).getCentroid();
        const Geometry::Point pc(*cellIt, rcCell);
        const auto& xc = pc.getPhysicalCoordinates();

        Real dot = 0.0;
        for (size_t c = 0; c < sdim; ++c)
          dot += cr(c) * (xc(c) - X0(c));
        orientation = (dot > 0.0) ? -1.0 : 1.0;
      }

      /// Nodal displacement of face vertex a, component c.
      Real nodal(const StateFES& fes, const State& d, size_t a, size_t c) const
      {
        return d[fes.getGlobalIndex({ 0, vertices[a] }, c)];
      }

      /// Deformed tangents at this face (P1: constant over the face).
      void deformedTangents(Math::SpatialMatrix<Real>& xjac,
                            const StateFES& fes, const State& d) const
      {
        xjac = Xjac;
        for (size_t j = 0; j < faceDim; ++j)
          for (size_t a = 0; a < nv; ++a)
          {
            const Real g = p1FaceBasisGrad(a, j, faceDim);
            if (g == 0.0)
              continue;
            for (size_t c = 0; c < sdim; ++c)
              xjac(c, j) += g * nodal(fes, d, a, c);
          }
      }
    };
  }

  /**
   * @brief Linear-form integrator: follower pressure residual
   * @f$ +p \int (\mathbf{x}_{,\xi} \times \mathbf{x}_{,\eta})
   *      \cdot \mathbf{w} @f$. Add it with "+" to the Problem (it already
   * carries the residual-convention sign of a pressure pushing along the
   * outward normal).
   */
  template <class TestFunctionType, class DisplacementType>
  class FollowerPressureForce final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;
      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::FollowerPressureForce expects a Rodin test function.");

      /// Pressure captured BY REFERENCE: all copies follow the variable.
      FollowerPressureForce(const Real& pressure, const TestType& v,
                            const StateType& displacement)
        : Parent(v),
          m_pressure(pressure),
          m_test(v),
          m_displacement(displacement),
          m_testfes(v.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace())
      {}

      FollowerPressureForce(const FollowerPressureForce&) = default;

      FollowerPressureForce& setPolytope(const Geometry::Polytope& face)
        final override
      {
        m_polytope = face;
        m_kin.setFace(face);
        const size_t sdim = m_kin.sdim;
        const size_t vdim = m_testfes.get().getVectorDimension();
        assert(vdim == sdim);

        const size_t testDofs = m_kin.nv * vdim;
        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        // Deformed tangents and oriented current normal-area vector
        // (constant over a P1 face).
        Math::SpatialMatrix<Real> xjac;
        m_kin.deformedTangents(xjac, m_statefes.get(),
                               m_displacement.get());
        Math::SpatialVector<Real> t1(sdim), t2(sdim), cr;
        for (size_t c = 0; c < sdim; ++c)
        {
          t1(c) = xjac(c, 0);
          t2(c) = (m_kin.faceDim == 2) ? xjac(c, 1) : 0.0;
        }
        Internal::surfaceCross(cr, t1, t2, sdim);

        const Real p = m_pressure.get();
        const auto& qf = QF::PolytopeQuadratureFormula::get(
            2, face.getGeometry());
        const size_t nqp = qf.getSize();

        for (size_t q = 0; q < nqp; ++q)
        {
          const Real wq = qf.getWeight(q);
          const auto& rc = qf.getPoint(q);
          for (size_t a = 0; a < m_kin.nv; ++a)
          {
            const Real phi = Internal::p1FaceBasis(a, rc, m_kin.faceDim);
            for (size_t c = 0; c < vdim; ++c)
              m_elemVec(a * vdim + c) +=
                wq * m_kin.orientation * p * cr(c) * phi;
          }
        }
        return *this;
      }

      ScalarType integrate(size_t te) final override { return m_elemVec(te); }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Boundary;
      }

      FollowerPressureForce* copy() const noexcept final override
      {
        return new FollowerPressureForce(*this);
      }

    private:
      std::reference_wrapper<const Real> m_pressure;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;

      Internal::FaceKinematics<StateFESType, StateType> m_kin;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  template <class TestFunctionType, class DisplacementType>
  FollowerPressureForce(const Real&, const TestFunctionType&,
                        const DisplacementType&)
    -> FollowerPressureForce<std::decay_t<TestFunctionType>,
                             std::decay_t<DisplacementType>>;

  /**
   * @brief Bilinear-form integrator: the EXACT load stiffness of the
   * follower pressure,
   * @f$ +p \int (\delta\mathbf{d}_{,\xi} \times \mathbf{x}_{,\eta}
   *    + \mathbf{x}_{,\xi} \times \delta\mathbf{d}_{,\eta})
   *    \cdot \mathbf{w} @f$.  Add it with "+" alongside the material
   * tangent.  NONSYMMETRIC.
   */
  template <class TrialFunctionType, class TestFunctionType,
            class DisplacementType>
  class FollowerPressureTangent final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
      using TrialType = TrialFunctionType;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;
      using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTrialFunction<TrialType>::Value,
        "Solid::FollowerPressureTangent expects a Rodin trial function.");
      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::FollowerPressureTangent expects a Rodin test function.");

      FollowerPressureTangent(const Real& pressure, const TrialType& u,
                              const TestType& v,
                              const StateType& displacement)
        : Parent(u, v),
          m_pressure(pressure),
          m_trial(u),
          m_test(v),
          m_displacement(displacement),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace())
      {}

      FollowerPressureTangent(const FollowerPressureTangent&) = default;

      FollowerPressureTangent& setPolytope(const Geometry::Polytope& face)
        final override
      {
        m_polytope = face;
        m_kin.setFace(face);
        const size_t sdim = m_kin.sdim;
        const size_t vdim = m_testfes.get().getVectorDimension();
        assert(vdim == sdim);
        const size_t nv = m_kin.nv;
        const size_t nDofs = nv * vdim;

        m_matrix.resize(nDofs, nDofs);
        m_matrix.setZero();

        Math::SpatialMatrix<Real> xjac;
        m_kin.deformedTangents(xjac, m_statefes.get(),
                               m_displacement.get());
        Math::SpatialVector<Real> x1(sdim), x2(sdim);
        for (size_t c = 0; c < sdim; ++c)
        {
          x1(c) = xjac(c, 0);
          x2(c) = (m_kin.faceDim == 2) ? xjac(c, 1) : 0.0;
        }

        const Real p = m_pressure.get();
        const auto& qf = QF::PolytopeQuadratureFormula::get(
            2, face.getGeometry());
        const size_t nqp = qf.getSize();

        Math::SpatialVector<Real> e(sdim), d1(sdim), d2(sdim), c1, c2, dcr(sdim);
        for (size_t q = 0; q < nqp; ++q)
        {
          const Real wq = qf.getWeight(q);
          const auto& rc = qf.getPoint(q);
          for (size_t b = 0; b < nv; ++b) // trial vertex
          {
            const Real g1 = Internal::p1FaceBasisGrad(b, 0, m_kin.faceDim);
            const Real g2 = (m_kin.faceDim == 2)
              ? Internal::p1FaceBasisGrad(b, 1, m_kin.faceDim) : 0.0;
            for (size_t cb = 0; cb < vdim; ++cb) // trial component
            {
              // delta x_{,1} = g1 e_cb ; delta x_{,2} = g2 e_cb
              e.setZero();
              e(cb) = 1.0;
              for (size_t c = 0; c < sdim; ++c)
              {
                d1(c) = g1 * e(c);
                d2(c) = g2 * e(c);
              }
              Internal::surfaceCross(c1, d1, x2, sdim);
              Internal::surfaceCross(c2, x1, d2, sdim);
              for (size_t c = 0; c < sdim; ++c)
                dcr(c) = c1(c) + ((m_kin.faceDim == 2) ? c2(c) : 0.0);

              const size_t tr = b * vdim + cb;
              for (size_t a = 0; a < nv; ++a) // test vertex
              {
                const Real phi =
                  Internal::p1FaceBasis(a, rc, m_kin.faceDim);
                for (size_t ca = 0; ca < vdim; ++ca)
                  m_matrix(a * vdim + ca, tr) +=
                    wq * m_kin.orientation * p * dcr(ca) * phi;
              }
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
        return Geometry::Region::Boundary;
      }

      FollowerPressureTangent* copy() const noexcept final override
      {
        return new FollowerPressureTangent(*this);
      }

    private:
      std::reference_wrapper<const Real> m_pressure;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const TrialFESType> m_trialfes;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;

      Internal::FaceKinematics<StateFESType, StateType> m_kin;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class TrialFunctionType, class TestFunctionType,
            class DisplacementType>
  FollowerPressureTangent(const Real&, const TrialFunctionType&,
                          const TestFunctionType&, const DisplacementType&)
    -> FollowerPressureTangent<std::decay_t<TrialFunctionType>,
                               std::decay_t<TestFunctionType>,
                               std::decay_t<DisplacementType>>;
}

#endif
