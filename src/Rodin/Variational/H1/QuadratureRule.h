#ifndef RODIN_VARIATIONAL_H1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_H1_QUADRATURERULE_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "H1.h"
#include "H1Element.h"

namespace Rodin::Variational
{
  template <
    size_t K,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
              H1<K, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
              H1<K, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, FESType, TrialSpace>>, FESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, FESType, TestSpace>>,  FESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point),
          m_dim(0),
          m_trCount(0),
          m_teCount(0)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_dim(std::exchange(other.m_dim, 0)),
          m_trCount(std::exchange(other.m_trCount, 0)),
          m_teCount(std::exchange(other.m_teCount, 0)),
          m_trRefGrad(std::move(other.m_trRefGrad)),
          m_teRefGrad(std::move(other.m_teRefGrad)),
          m_mat(std::move(other.m_mat))
      {}

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set = true;
          m_order = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        // Symmetry applicability for this specialization:
        const bool symmetric =
          (&trialfes == &testfes) && (ntr == nte);
        // If you can assert stronger (same FE object / same ordering), do it here.

        // Row-major local matrix for fast row writes.
        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data();

        const auto& trTab = trialfe.getTabulation(*m_qf);
        const auto& teTab = testfe .getTabulation(*m_qf);

        // Scratch for physical gradients at this qp.
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gtr;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> Gte;

        if (Gtr.size() < ntr) Gtr.resize(ntr);
        if (Gte.size() < nte) Gte.resize(nte);
        for (size_t a = 0; a < ntr; ++a) Gtr[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < nte; ++b) Gte[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Build physical gradients for trial and test
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gtr[a][0] = a00*gx + a10*gy + a20*gz;
              Gtr[a][1] = a01*gx + a11*gy + a21*gz;
              Gtr[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              Gte[b][0] = a00*gx + a10*gy + a20*gz;
              Gte[b][1] = a01*gx + a11*gy + a21*gz;
              Gte[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              Gtr[a][0] = a00*gx + a10*gy;
              Gtr[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              Gte[b][0] = a00*gx + a10*gy;
              Gte[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < ntr; ++a)
            {
              const auto g = trTab.getGradient(qp, a);
              Gtr[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < nte; ++b)
            {
              const auto g = teTab.getGradient(qp, b);
              Gte[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          if (symmetric)
          {
            // Accumulate only lower triangle (i,j with j <= i)
            // Use trial index as column, test index as row: A[row*ntr + col]
            for (size_t i = 0; i < nte; ++i)
            {
              ScalarType* row = A + i * ntr;
              const auto& gi = Gte[i];

              // diagonal
              row[i] += wdet * Math::dot(gi, Gtr[i]);

              for (size_t j = 0; j < i; ++j)
                row[j] += wdet * Math::dot(gi, Gtr[j]);
            }
          }
          else
          {
            // General (non-symmetric) accumulation
            for (size_t b = 0; b < nte; ++b)
            {
              ScalarType* row = A + b * ntr;
              const auto& gb = Gte[b];
              for (size_t a = 0; a < ntr; ++a)
                row[a] += wdet * Math::dot(gb, Gtr[a]);
            }
          }
        }

        if (symmetric)
        {
          // Mirror lower -> upper
          // m_mat is row-major, but Eigen triangular views still work.
          // We filled (i,i) and (i,j) for j<i. Mirror to upper.
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      // sizes
      size_t m_dim;
      size_t m_qps;
      size_t m_trCount;
      size_t m_teCount;

      // reference gradients (qp-major)
      std::vector<Math::SpatialVector<ScalarType>> m_trRefGrad;
      std::vector<Math::SpatialVector<ScalarType>> m_teRefGrad;

      // local matrix (rows=test, cols=trial)
      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
    };

  // CTAD helper
  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>&)
  -> QuadratureRule<
       Dot<
         ShapeFunctionBase<
           Grad<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
           H1<K, Scalar, Mesh>, TrialSpace>,
         ShapeFunctionBase<
           Grad<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
           H1<K, Scalar, Mesh>, TestSpace>>>;

  template <
    size_t K,
    class LHSDerived, class RHSDerived,
    class Scalar, class Mesh>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>>
    : public LocalBilinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Dot<
            ShapeFunctionBase<
              Jacobian<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
              H1<K, Scalar, Mesh>, TrialSpace>,
            ShapeFunctionBase<
              Jacobian<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
              H1<K, Scalar, Mesh>, TestSpace>>>::ScalarType>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;

      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, FESType, TrialSpace>>,
          FESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, FESType, TestSpace>>,
          FESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr),
          m_polytope(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point)),
          m_mat(std::move(other.m_mat))
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();

        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes .getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe .getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute_qf = (!m_set || m_order != order || m_geometry != geometry);

        if (recompute_qf)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);

        // --- infer vdim from your H1Element<K, Math::Vector<Scalar>> convention ---
        // scalarCount is purely geometry+K (no vdim)
        const H1Element<K, ScalarType> trialScalarElement(polytope.getGeometry());
        const H1Element<K, ScalarType> testScalarElement(polytope.getGeometry());
        const size_t scalarCount = trialScalarElement.getCount();

        assert(scalarCount > 0);
        assert(ntr % scalarCount == 0);
        assert(nte % scalarCount == 0);

        const size_t vdim_tr = ntr / scalarCount;
        // const size_t vdim_te = nte / scalarCount;
        assert(vdim_tr == nte / scalarCount);
        const size_t vdim = vdim_tr;

        const bool symmetric =
          (&trialfes == &testfes) && (ntr == nte);

        m_mat.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_mat.setZero();
        ScalarType* A = m_mat.data(); // row-major (rows=test, cols=trial)

        // Use scalar tabulations (fast, cached in your H1Element<K, Scalar>::getTabulation)
        const auto& trTabS = trialScalarElement.getTabulation(*m_qf);
        const auto& teTabS = testScalarElement.getTabulation(*m_qf);

        static thread_local std::vector<Math::SpatialVector<ScalarType>> GtrS;
        static thread_local std::vector<Math::SpatialVector<ScalarType>> GteS;

        if (GtrS.size() < scalarCount) GtrS.resize(scalarCount);
        if (GteS.size() < scalarCount) GteS.resize(scalarCount);
        for (size_t a = 0; a < scalarCount; ++a) GtrS[a].resize(static_cast<std::uint8_t>(d));
        for (size_t b = 0; b < scalarCount; ++b) GteS[b].resize(static_cast<std::uint8_t>(d));

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          // Map scalar reference gradients to physical gradients (once per scalar DOF)
          if (d == 3)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0), a20 = Jinv(2,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1), a21 = Jinv(2,1);
            const ScalarType a02 = Jinv(0,2), a12 = Jinv(1,2), a22 = Jinv(2,2);

            for (size_t a = 0; a < scalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GtrS[a][0] = a00*gx + a10*gy + a20*gz;
              GtrS[a][1] = a01*gx + a11*gy + a21*gz;
              GtrS[a][2] = a02*gx + a12*gy + a22*gz;
            }
            for (size_t b = 0; b < scalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1], gz = g[2];
              GteS[b][0] = a00*gx + a10*gy + a20*gz;
              GteS[b][1] = a01*gx + a11*gy + a21*gz;
              GteS[b][2] = a02*gx + a12*gy + a22*gz;
            }
          }
          else if (d == 2)
          {
            const ScalarType a00 = Jinv(0,0), a10 = Jinv(1,0);
            const ScalarType a01 = Jinv(0,1), a11 = Jinv(1,1);

            for (size_t a = 0; a < scalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              const ScalarType gx = g[0], gy = g[1];
              GtrS[a][0] = a00*gx + a10*gy;
              GtrS[a][1] = a01*gx + a11*gy;
            }
            for (size_t b = 0; b < scalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              const ScalarType gx = g[0], gy = g[1];
              GteS[b][0] = a00*gx + a10*gy;
              GteS[b][1] = a01*gx + a11*gy;
            }
          }
          else if (d == 1)
          {
            const ScalarType a00 = Jinv(0,0);
            for (size_t a = 0; a < scalarCount; ++a)
            {
              const auto g = trTabS.getGradient(qp, a);
              GtrS[a][0] = a00 * g[0];
            }
            for (size_t b = 0; b < scalarCount; ++b)
            {
              const auto g = teTabS.getGradient(qp, b);
              GteS[b][0] = a00 * g[0];
            }
          }
          else
          {
            assert(false);
          }

          // Assemble: for each component c, add the same scalar matrix into block (c,c)
          // Local DOF mapping is exactly your convention: dof = scalarIndex*vdim + c.
          if (symmetric)
          {
            // lower triangle in vector-dof space corresponds to scalar lower + component blocks
            for (size_t ib = 0; ib < scalarCount; ++ib)
            {
              const auto& gb = GteS[ib];

              for (size_t ia = 0; ia <= ib; ++ia)
              {
                const ScalarType kij = wdet * Math::dot(gb, GtrS[ia]);
                if (kij == ScalarType(0))
                  continue;

                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t row = ib * vdim + c;
                  const size_t col = ia * vdim + c;
                  A[row * ntr + col] += kij;
                }
              }
            }
          }
          else
          {
            for (size_t ib = 0; ib < scalarCount; ++ib)
            {
              const auto& gb = GteS[ib];

              for (size_t ia = 0; ia < scalarCount; ++ia)
              {
                const ScalarType kij = wdet * Math::dot(gb, GtrS[ia]);
                if (kij == ScalarType(0))
                  continue;

                for (size_t c = 0; c < vdim; ++c)
                {
                  const size_t row = ib * vdim + c;
                  const size_t col = ia * vdim + c;
                  A[row * ntr + col] += kij;
                }
              }
            }
          }
        }

        if (symmetric)
        {
          m_mat.template triangularView<Eigen::Upper>() =
            m_mat.transpose().template triangularView<Eigen::Upper>();
        }

        return *this;
      }

      inline ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_mat(te, tr);
      }

      virtual Geometry::Region getRegion() const override = 0;
      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;

      Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_mat;
  };

  // CTAD helper
  template <size_t K, class LHSDerived, class RHSDerived, class Scalar, class Mesh>
  QuadratureRule(
    const Dot<
      ShapeFunctionBase<
        Jacobian<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
        H1<K, Scalar, Mesh>, TrialSpace>,
      ShapeFunctionBase<
        Jacobian<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
        H1<K, Scalar, Mesh>, TestSpace>>&)
  -> QuadratureRule<
       Dot<
         ShapeFunctionBase<
           Jacobian<ShapeFunction<LHSDerived, H1<K, Scalar, Mesh>, TrialSpace>>,
           H1<K, Scalar, Mesh>, TrialSpace>,
         ShapeFunctionBase<
           Jacobian<ShapeFunction<RHSDerived, H1<K, Scalar, Mesh>, TestSpace>>,
           H1<K, Scalar, Mesh>, TestSpace>>>;
}

#endif
