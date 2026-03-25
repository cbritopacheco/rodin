#!/usr/bin/env python3
"""Replace 11 specializations in P1/QuadratureRule.h with GenericPolytopeQuadrature pattern."""

import sys

FILE = "/home/runner/work/rodin/rodin/src/Rodin/Variational/P1/QuadratureRule.h"

with open(FILE, 'r') as f:
    lines = f.readlines()

# Replacements: (start_line, end_line, new_content)  -- 1-indexed, inclusive
# We apply bottom-to-top so line numbers stay valid.

replacements = []

# ============================================================
# Spec 1 (lines 79-193): ∫ v dx
# ============================================================
spec1 = """\
    public:
      using FESType = P1<Range, Mesh>;

      using IntegrandType =
        ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>>;

      using IntegrandRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_vec(std::move(other.m_vec))
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
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order =
          integrand.getOrder(polytope).value_or(fe.getOrder());

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
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

        const size_t nte = integrand.getDOFs(polytope);

        const P1Element<ScalarType> scalarFE(geometry);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          if constexpr (std::is_same_v<IntegrandRangeType, ScalarType>)
          {
            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * scalarFE.getBasis(local)(m_qf->getPoint(qp));
          }
          else if constexpr (std::is_same_v<IntegrandRangeType, Math::Vector<ScalarType>>)
          {
            const size_t vdim = fes.getVectorDimension();
            assert(nte == scalarFE.getCount() * vdim);

            for (size_t local = 0; local < nte; ++local)
            {
              const size_t scalarLocal = local / vdim;
              m_vec(local) += wdet * scalarFE.getBasis(scalarLocal)(m_qf->getPoint(qp));
            }
          }
          else
          {
            static_assert(
              std::is_same_v<IntegrandRangeType, ScalarType>
              || std::is_same_v<IntegrandRangeType, Math::Vector<ScalarType>>,
              "Unsupported P1 Integral(v) range type.");
          }
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
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

      Math::Vector<ScalarType> m_vec;
"""
replacements.append((79, 193, spec1))

# ============================================================
# Spec 2 (lines 241-374): ∫ f·v dx
# ============================================================
spec2 = """\
    public:
      using FESType =
        P1<Range, Mesh>;

      using LHSType =
        FunctionBase<LHSDerived>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;

      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      using IntegrandType =
        ShapeFunctionBase<Dot<LHSType, RHSType>>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_vec(std::move(other.m_vec))
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
        static thread_local LHSRangeType s_v;
        m_polytope = &polytope;

        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        auto& integrand = *m_integrand;
        const auto& f   = integrand.getDerived().getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t order =
          integrand.getOrder(polytope).value_or(fe.getOrder());

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
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

        const size_t nte = integrand.getDOFs(polytope);

        const P1Element<ScalarType> scalarFE(geometry);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        if constexpr (std::is_same_v<RHSRangeType, ScalarType>)
        {
          for (size_t qp = 0; qp < m_ps.size(); ++qp)
          {
            const auto& p = m_ps[qp];
            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const ScalarType fval = f(p);

            for (size_t local = 0; local < nte; ++local)
              m_vec(local) += wdet * fval * scalarFE.getBasis(local)(m_qf->getPoint(qp));
          }
        }
        else if constexpr (std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>)
        {
          const size_t vdim = fes.getVectorDimension();
          const size_t scalarCount = scalarFE.getCount();
          assert(nte == scalarCount * vdim);

          for (size_t qp = 0; qp < m_ps.size(); ++qp)
          {
            const auto& p = m_ps[qp];
            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            s_v = f(p);

            for (size_t local = 0; local < nte; ++local)
            {
              const size_t scalarLocal = local / vdim;
              const size_t comp = local % vdim;
              m_vec(local) +=
                wdet * s_v.coeff(comp) * scalarFE.getBasis(scalarLocal)(m_qf->getPoint(qp));
            }
          }
        }
        else
        {
          static_assert(
            std::is_same_v<RHSRangeType, ScalarType>
            || std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>,
            "Unsupported P1 Integral(f.v) range type.");
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_vec(local);
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

      Math::Vector<ScalarType> m_vec;
"""
replacements.append((241, 374, spec2))

# ============================================================
# Spec 4 (lines 688-929): ∫ (Au)·v dx - Coefficient mass bilinear
# ============================================================
spec4 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType =
        typename FormLanguage::Traits<MultiplicandType>::RangeType;

      using LHSRangeType =
        typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType =
        typename FormLanguage::Traits<RHSType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t order = trialfe.getOrder() + testfe.getOrder();

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
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

        const P1Element<ScalarType> trialScalarFE(geometry);
        const P1Element<ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& rc = m_qf->getPoint(qp);

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);

            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = testScalarFE.getBasis(ib)(rc);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType kij =
                  wdet * csv * phi_te * trialScalarFE.getBasis(ia)(rc);
                for (size_t c = 0; c < vdim; ++c)
                  m_matrix(ib * vdim + c, ia * vdim + c) += kij;
              }
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            static_assert(std::is_same_v<MultiplicandRangeType, Math::Vector<ScalarType>>);
            static_assert(std::is_same_v<RHSRangeType, Math::Vector<ScalarType>>);

            static thread_local Math::Matrix<ScalarType> s_cmv;
            coeff.getValue(s_cmv, p);

            for (size_t ib = 0; ib < scalarCountTe; ++ib)
            {
              const ScalarType phi_te = testScalarFE.getBasis(ib)(rc);
              for (size_t ia = 0; ia < scalarCountTr; ++ia)
              {
                const ScalarType phi_tr = trialScalarFE.getBasis(ia)(rc);
                const ScalarType w = wdet * phi_te * phi_tr;
                for (size_t ci = 0; ci < vdim; ++ci)
                  for (size_t cj = 0; cj < vdim; ++cj)
                    m_matrix(ib * vdim + ci, ia * vdim + cj) += w * s_cmv(ci, cj);
              }
            }
          }
          else
          {
            assert(false);
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((688, 929, spec4))

# ============================================================
# Spec 5 (lines 994-1152): ∫ ∇u·∇v dx - Isotropic stiffness
# ============================================================
spec5 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_refGrad(std::move(other.m_refGrad)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const P1Element<ScalarType> fe(geometry);
          const size_t n = fe.getCount();

          m_refGrad.resize(n);

          const auto& rc = m_qf->getPoint(0);
          for (size_t local = 0; local < n; ++local)
          {
            auto& g = m_refGrad[local];
            g.resize(d);
            const auto& basis = fe.getBasis(local);
            for (size_t j = 0; j < d; ++j)
              g(j) = basis.template getDerivative<1>(j)(rc);
          }

          m_matrix.resize(n, n);
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t n = m_refGrad.size();

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto G = Jinv * Jinv.transpose();

          for (size_t i = 0; i < n; ++i)
          {
            const auto Ggi = G * m_refGrad[i];
            m_matrix(i, i) += wdet * Math::dot(m_refGrad[i], Ggi);
            for (size_t j = 0; j < i; ++j)
              m_matrix(i, j) += wdet * Math::dot(m_refGrad[j], Ggi);
          }
        }

        m_matrix.template triangularView<Eigen::Upper>() =
          m_matrix.template triangularView<Eigen::Lower>().transpose();

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<Math::SpatialVector<ScalarType>> m_refGrad;
      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((994, 1152, spec5))

# ============================================================
# Spec 6 (lines 1220-1481): ∫ (A∇u)·∇v dx - Anisotropic stiffness
# ============================================================
spec6 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using MultiplicandType =
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using MultiplicandOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandRangeType =
        typename FormLanguage::Traits<MultiplicandType>::RangeType;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using MultiplicandOperandRangeType =
        typename FormLanguage::Traits<MultiplicandOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<MultiplicandOperandRangeType, ScalarType>);
      static_assert(std::is_same_v<RHSOperandRangeType, ScalarType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_trialRefGrad(std::move(other.m_trialRefGrad)),
          m_testRefGrad(std::move(other.m_testRefGrad)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const P1Element<ScalarType> trialScalarFE(geometry);
          const P1Element<ScalarType> testScalarFE(geometry);

          const auto& rc = m_qf->getPoint(0);

          m_trialRefGrad.resize(trialScalarFE.getCount());
          for (size_t local = 0; local < trialScalarFE.getCount(); ++local)
          {
            auto& g = m_trialRefGrad[local];
            g.resize(d);
            const auto& basis = trialScalarFE.getBasis(local);
            for (size_t j = 0; j < d; ++j)
              g(j) = basis.template getDerivative<1>(j)(rc);
          }

          if (trialfes == testfes)
          {
            m_testRefGrad = m_trialRefGrad;
          }
          else
          {
            m_testRefGrad.resize(testScalarFE.getCount());
            for (size_t local = 0; local < testScalarFE.getCount(); ++local)
            {
              auto& g = m_testRefGrad[local];
              g.resize(d);
              const auto& basis = testScalarFE.getBasis(local);
              for (size_t j = 0; j < d; ++j)
                g(j) = basis.template getDerivative<1>(j)(rc);
            }
          }

          m_matrix.resize(m_testRefGrad.size(), m_trialRefGrad.size());
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto G = Jinv * Jinv.transpose();

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);

            if (trialfes == testfes)
            {
              const size_t n = m_trialRefGrad.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ggi = G * m_trialRefGrad[i];
                m_matrix(i, i) += wdet * csv * Math::dot(m_trialRefGrad[i], Ggi);
                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * csv * Math::dot(m_trialRefGrad[j], Ggi);
              }
            }
            else
            {
              const size_t ntr = m_trialRefGrad.size();
              const size_t nte = m_testRefGrad.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Ggte = G * m_testRefGrad[te];
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * csv * Math::dot(m_trialRefGrad[tr], Ggte);
              }
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            static thread_local Math::Matrix<ScalarType> s_cmv;
            coeff.getValue(s_cmv, p);

            if (trialfes == testfes)
            {
              const size_t n = m_trialRefGrad.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto AGgi = s_cmv * (G * m_trialRefGrad[i]);
                m_matrix(i, i) += wdet * Math::dot(AGgi, m_trialRefGrad[i]);
                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * Math::dot(AGgi, m_trialRefGrad[j]);
              }

              // Fill the upper triangle for the non-symmetric coefficient case
              for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                  m_matrix(i, j) += wdet * Math::dot(
                    s_cmv * (G * m_trialRefGrad[j]), m_trialRefGrad[i]);
            }
            else
            {
              const size_t ntr = m_trialRefGrad.size();
              const size_t nte = m_testRefGrad.size();

              for (size_t te = 0; te < nte; ++te)
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * Math::dot(
                    s_cmv * (G * m_trialRefGrad[tr]), m_testRefGrad[te]);
            }
          }
          else
          {
            assert(false);
          }
        }

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          if (trialfes == testfes)
          {
            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<Math::SpatialVector<ScalarType>> m_trialRefGrad;
      std::vector<Math::SpatialVector<ScalarType>> m_testRefGrad;

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((1220, 1481, spec6))

# ============================================================
# Spec 7 (lines 1543-1733): ∫ c(u·v) dx - Coefficient × mass
# ============================================================
spec7 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;
      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>;

      using RHSType =
        ShapeFunctionBase<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>;

      using InnerIntegrandType = Dot<LHSType, RHSType>;
      using IntegrandType = Mult<CoefficientType, InnerIntegrandType>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
      using ScalarType   = typename FormLanguage::Traits<InnerIntegrandType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getRHS().getLHS().getLeaf(), integrand.getRHS().getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& coeff = integrand.getLHS();
        const auto& inner = integrand.getRHS();
        const auto& lhs = inner.getLHS();
        const auto& rhs = inner.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t order = trialfe.getOrder() + testfe.getOrder();

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
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

        const P1Element<ScalarType> trialScalarFE(geometry);
        const P1Element<ScalarType> testScalarFE(geometry);
        const size_t scalarCountTr = trialScalarFE.getCount();
        const size_t scalarCountTe = testScalarFE.getCount();

        assert(scalarCountTr > 0 && ntr % scalarCountTr == 0);
        assert(scalarCountTe > 0 && nte % scalarCountTe == 0);
        const size_t vdim = ntr / scalarCountTr;
        assert(vdim == nte / scalarCountTe);

        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& rc = m_qf->getPoint(qp);

          const ScalarType csv = coeff(p);

          for (size_t ib = 0; ib < scalarCountTe; ++ib)
          {
            const ScalarType phi_te = testScalarFE.getBasis(ib)(rc);
            for (size_t ia = 0; ia < scalarCountTr; ++ia)
            {
              const ScalarType kij =
                wdet * csv * phi_te * trialScalarFE.getBasis(ia)(rc);
              for (size_t c = 0; c < vdim; ++c)
                m_matrix(ib * vdim + c, ia * vdim + c) += kij;
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((1543, 1733, spec7))

# ============================================================
# Spec 8 (lines 1784-1953): ∫ (∇·u)q dx - Div-pressure coupling
# ============================================================
spec8 = """\
    public:
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Math::Vector<Real>, LHSMesh>;
      using TestFESType  = P1<Real, RHSMesh>;

      using LHSType = ShapeFunctionBase<
        Div<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
        TrialFESType, TrialSpace>;

      using RHSType = ShapeFunctionBase<
        ShapeFunction<RHSDerived, TestFESType, TestSpace>,
        TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_refGrad(std::move(other.m_refGrad)),
          m_testBasis(std::move(other.m_testBasis)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const auto& rc = m_qf->getPoint(0);

          const P1Element<Math::Vector<Real>> trialVecFE(geometry, trialfes.getVectorDimension());
          const P1Element<Real> testScalarFE(geometry);

          m_testBasis.resize(testScalarFE.getCount());
          for (size_t i = 0; i < testScalarFE.getCount(); ++i)
            m_testBasis[i] = testScalarFE.getBasis(i)(rc);

          m_refGrad.resize(trialVecFE.getCount());
          for (size_t i = 0; i < trialVecFE.getCount(); ++i)
          {
            auto& refg = m_refGrad[i];
            refg.resize(trialfes.getVectorDimension());
            for (size_t comp = 0; comp < trialfes.getVectorDimension(); ++comp)
            {
              refg[comp].resize(d);
              const auto& basis = trialVecFE.getBasis(i);
              for (size_t j = 0; j < d; ++j)
                refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
            }
          }

          m_matrix.resize(m_testBasis.size(), m_refGrad.size());
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const size_t vdim = m_refGrad.empty() ? 0 : m_refGrad.front().size();

          for (size_t i = 0; i < m_refGrad.size(); ++i)
          {
            ScalarType div = 0;
            for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
            {
              ScalarType physComp = 0;
              for (size_t j = 0; j < d; ++j)
                physComp += Jinv(comp, j) * m_refGrad[i][comp](j);
              div += physComp;
            }

            const size_t nte = m_testBasis.size();
            for (size_t te = 0; te < nte; ++te)
              m_matrix(te, i) += wdet * m_testBasis[te] * div;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;
      std::vector<ScalarType> m_testBasis;
      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((1784, 1953, spec8))

# ============================================================
# Spec 9 (lines 2004-2172): ∫ p(∇·v) dx - Pressure-divergence
# ============================================================
spec9 = """\
    public:
      using ScalarType = typename FormLanguage::Traits<P1<Real, LHSMesh>>::ScalarType;
      using TrialFESType = P1<Real, LHSMesh>;
      using TestFESType  = P1<Math::Vector<Real>, RHSMesh>;

      using LHSType = ShapeFunctionBase<
        ShapeFunction<LHSDerived, TrialFESType, TrialSpace>,
        TrialFESType, TrialSpace>;

      using RHSType = ShapeFunctionBase<
        Div<ShapeFunction<RHSDerived, TestFESType, TestSpace>>,
        TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_refGrad(std::move(other.m_refGrad)),
          m_trialBasis(std::move(other.m_trialBasis)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes  = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const auto& rc = m_qf->getPoint(0);

          const P1Element<Real> trialScalarFE(geometry);
          const P1Element<Math::Vector<Real>> testVecFE(geometry, testfes.getVectorDimension());

          m_trialBasis.resize(trialScalarFE.getCount());
          for (size_t i = 0; i < trialScalarFE.getCount(); ++i)
            m_trialBasis[i] = trialScalarFE.getBasis(i)(rc);

          m_refGrad.resize(testVecFE.getCount());
          for (size_t i = 0; i < testVecFE.getCount(); ++i)
          {
            auto& refg = m_refGrad[i];
            refg.resize(testfes.getVectorDimension());
            for (size_t comp = 0; comp < testfes.getVectorDimension(); ++comp)
            {
              refg[comp].resize(d);
              const auto& basis = testVecFE.getBasis(i);
              for (size_t j = 0; j < d; ++j)
                refg[comp](j) = basis.template getDerivative<1>(comp, j)(rc);
            }
          }

          m_matrix.resize(m_refGrad.size(), m_trialBasis.size());
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const size_t vdim = m_refGrad.empty() ? 0 : m_refGrad.front().size();

          for (size_t i = 0; i < m_refGrad.size(); ++i)
          {
            ScalarType div = 0;
            for (size_t comp = 0; comp < std::min(vdim, d); ++comp)
            {
              ScalarType physComp = 0;
              for (size_t j = 0; j < d; ++j)
                physComp += Jinv(comp, j) * m_refGrad[i][comp](j);
              div += physComp;
            }

            const size_t ntr = m_trialBasis.size();
            for (size_t tr = 0; tr < ntr; ++tr)
              m_matrix(i, tr) += wdet * m_trialBasis[tr] * div;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<std::vector<Math::SpatialVector<ScalarType>>> m_refGrad;
      std::vector<ScalarType> m_trialBasis;
      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((2004, 2172, spec9))

# ============================================================
# Spec 10 (lines 2232-2446): ∫ Ju:Jv dx - Isotropic Jacobian form
# ============================================================
spec10 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using LHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, RHSFESType, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, RHSFESType, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);
      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_trialRefJac(std::move(other.m_trialRefJac)),
          m_testRefJac(std::move(other.m_testRefJac)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const P1Element<Math::Vector<ScalarType>> trialVecFE(
            geometry, trialfes.getVectorDimension());
          const P1Element<Math::Vector<ScalarType>> testVecFE(
            geometry, testfes.getVectorDimension());

          const auto& rc = m_qf->getPoint(0);

          m_trialRefJac.resize(trialVecFE.getCount());
          for (size_t local = 0; local < trialVecFE.getCount(); ++local)
          {
            auto& J = m_trialRefJac[local];
            J.resize(trialfes.getVectorDimension(), d);
            const auto& basis = trialVecFE.getBasis(local);
            for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
              for (size_t j = 0; j < d; ++j)
                J(i, j) = basis.template getDerivative<1>(i, j)(rc);
          }

          if (trialfes == testfes)
          {
            m_testRefJac = m_trialRefJac;
          }
          else
          {
            m_testRefJac.resize(testVecFE.getCount());
            for (size_t local = 0; local < testVecFE.getCount(); ++local)
            {
              auto& J = m_testRefJac[local];
              J.resize(testfes.getVectorDimension(), d);
              const auto& basis = testVecFE.getBasis(local);
              for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }
          }

          m_matrix.resize(m_testRefJac.size(), m_trialRefJac.size());
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();

          if (trialfes == testfes)
          {
            const size_t n = m_trialRefJac.size();
            for (size_t i = 0; i < n; ++i)
            {
              const auto Ji = m_trialRefJac[i] * Jinv;
              m_matrix(i, i) += wdet * Ji.squaredNorm();

              for (size_t j = 0; j < i; ++j)
                m_matrix(i, j) += wdet * Math::dot(m_trialRefJac[j] * Jinv, Ji);
            }
          }
          else
          {
            const size_t ntr = m_trialRefJac.size();
            const size_t nte = m_testRefJac.size();

            for (size_t te = 0; te < nte; ++te)
            {
              const auto Jte = m_testRefJac[te] * Jinv;
              for (size_t tr = 0; tr < ntr; ++tr)
                m_matrix(te, tr) += wdet * Math::dot(m_trialRefJac[tr] * Jinv, Jte);
            }
          }
        }

        if (trialfes == testfes)
        {
          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.template triangularView<Eigen::Lower>().transpose();
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<Math::SpatialMatrix<ScalarType>> m_trialRefJac, m_testRefJac;

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((2232, 2446, spec10))

# ============================================================
# Spec 11 (lines 2523-2791): ∫ (AJu):Jv dx - Anisotropic Jacobian form
# ============================================================
spec11 = """\
    public:
      using LHSFESType = P1<LHSRange, LHSMesh>;

      using RHSFESType = P1<RHSRange, RHSMesh>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using CoefficientRangeType =
        typename FormLanguage::Traits<CoefficientType>::RangeType;

      using MultiplicandType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, LHSFESType, TrialSpace>>>;

      using LHSType =
        ShapeFunctionBase<Mult<CoefficientType, MultiplicandType>>;

      using LHSOperandType =
        ShapeFunction<LHSDerived, LHSFESType, TrialSpace>;

      using LHSOperandRangeType =
        typename FormLanguage::Traits<LHSOperandType>::RangeType;

      using RHSType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>>>;

      using RHSOperandType =
        ShapeFunction<RHSDerived, P1<RHSRange, RHSMesh>, TestSpace>;

      using RHSOperandRangeType =
        typename FormLanguage::Traits<RHSOperandType>::RangeType;

      using IntegrandType = Dot<LHSType, RHSType>;

      using IntegrandRangeType =
        typename FormLanguage::Traits<IntegrandType>::RangeType;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

      static_assert(std::is_same_v<LHSOperandRangeType, Math::Vector<ScalarType>>);
      static_assert(std::is_same_v<RHSOperandRangeType, Math::Vector<ScalarType>>);

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_trialRefJac(std::move(other.m_trialRefJac)),
          m_testRefJac(std::move(other.m_testRefJac)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& coeff = lhs.getDerived().getLHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const P1Element<Math::Vector<ScalarType>> trialVecFE(
            geometry, trialfes.getVectorDimension());
          const P1Element<Math::Vector<ScalarType>> testVecFE(
            geometry, testfes.getVectorDimension());

          const auto& rc = m_qf->getPoint(0);

          m_trialRefJac.resize(trialVecFE.getCount());
          for (size_t local = 0; local < trialVecFE.getCount(); ++local)
          {
            auto& J = m_trialRefJac[local];
            J.resize(trialfes.getVectorDimension(), d);
            const auto& basis = trialVecFE.getBasis(local);
            for (size_t i = 0; i < trialfes.getVectorDimension(); ++i)
              for (size_t j = 0; j < d; ++j)
                J(i, j) = basis.template getDerivative<1>(i, j)(rc);
          }

          if (trialfes == testfes)
          {
            m_testRefJac = m_trialRefJac;
          }
          else
          {
            m_testRefJac.resize(testVecFE.getCount());
            for (size_t local = 0; local < testVecFE.getCount(); ++local)
            {
              auto& J = m_testRefJac[local];
              J.resize(testfes.getVectorDimension(), d);
              const auto& basis = testVecFE.getBasis(local);
              for (size_t i = 0; i < testfes.getVectorDimension(); ++i)
                for (size_t j = 0; j < d; ++j)
                  J(i, j) = basis.template getDerivative<1>(i, j)(rc);
            }
          }

          m_matrix.resize(m_testRefJac.size(), m_trialRefJac.size());
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();

          if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
          {
            const ScalarType csv = coeff.getValue(p);

            if (trialfes == testfes)
            {
              const size_t n = m_trialRefJac.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ji = m_trialRefJac[i] * Jinv;
                m_matrix(i, i) += wdet * csv * Ji.squaredNorm();

                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * csv * Math::dot(m_trialRefJac[j] * Jinv, Ji);
              }
            }
            else
            {
              const size_t ntr = m_trialRefJac.size();
              const size_t nte = m_testRefJac.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Jte = m_testRefJac[te] * Jinv;
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * csv * Math::dot(m_trialRefJac[tr] * Jinv, Jte);
              }
            }
          }
          else if constexpr (std::is_same_v<CoefficientRangeType, Math::Matrix<ScalarType>>)
          {
            static thread_local Math::Matrix<ScalarType> s_cmv;
            coeff.getValue(s_cmv, p);

            if (trialfes == testfes)
            {
              const size_t n = m_trialRefJac.size();
              for (size_t i = 0; i < n; ++i)
              {
                const auto Ji = m_trialRefJac[i] * Jinv;
                m_matrix(i, i) += wdet * Math::dot(s_cmv * Ji, Ji);

                for (size_t j = 0; j < i; ++j)
                  m_matrix(i, j) += wdet * Math::dot(s_cmv * (m_trialRefJac[j] * Jinv), Ji);
              }

              for (size_t i = 0; i < n; ++i)
                for (size_t j = i + 1; j < n; ++j)
                  m_matrix(i, j) += wdet * Math::dot(
                    s_cmv * (m_trialRefJac[j] * Jinv), m_trialRefJac[i] * Jinv);
            }
            else
            {
              const size_t ntr = m_trialRefJac.size();
              const size_t nte = m_testRefJac.size();

              for (size_t te = 0; te < nte; ++te)
              {
                const auto Jte = m_testRefJac[te] * Jinv;
                for (size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(te, tr) += wdet * Math::dot(
                    s_cmv * (m_trialRefJac[tr] * Jinv), Jte);
              }
            }
          }
          else
          {
            assert(false);
          }
        }

        if constexpr (std::is_same_v<CoefficientRangeType, ScalarType>)
        {
          if (trialfes == testfes)
          {
            m_matrix.template triangularView<Eigen::Upper>() =
              m_matrix.template triangularView<Eigen::Lower>().transpose();
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      std::vector<Math::SpatialMatrix<ScalarType>> m_trialRefJac, m_testRefJac;

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((2523, 2791, spec11))

# ============================================================
# Spec 12 (lines 2871-3057): ∫ (Ju·f)·v dx - Linearized convection
# ============================================================
spec12 = """\
    public:
      using TrialFESType = P1<LHSRange, LHSMesh>;
      using TestFESType  = P1<RHSRange, RHSMesh>;

      using TrialSFType =
        ShapeFunctionBase<
          Jacobian<ShapeFunction<LHSDerived, TrialFESType, TrialSpace>>,
          TrialFESType, TrialSpace>;

      using CoefficientType = FunctionBase<CoefficientDerived>;

      using LHSType =
        ShapeFunctionBase<
          Mult<TrialSFType, CoefficientType>,
          TrialFESType, TrialSpace>;

      using RHSType =
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, TestFESType, TestSpace>,
          TestFESType, TestSpace>;

      using IntegrandType = Dot<LHSType, RHSType>;

      using ScalarType = typename FormLanguage::Traits<IntegrandType>::ScalarType;
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_qf(nullptr), m_polytope(nullptr),
          m_set(false), m_order(0),
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
          m_refGrad(std::move(other.m_refGrad)),
          m_basis(std::move(other.m_basis)),
          m_matrix(std::move(other.m_matrix))
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

        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = rhs.getFiniteElementSpace().getFiniteElement(d, idx);

        const size_t k_tr = trialfe.getOrder();
        const size_t k_te = testfe.getOrder();
        const size_t order = (k_tr == 0 || k_te == 0) ? 0 : (k_tr + k_te - 2);

        const auto geometry = polytope.getGeometry();
        const bool recompute = !m_set || m_order != order || m_geometry != geometry;

        // The coefficient is the RHS of the Mult node
        const auto& coeff = lhs.getDerived().getRHS();

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));

          const P1Element<ScalarType> scalarFE(geometry);
          const size_t n = scalarFE.getCount();

          const auto& rc = m_qf->getPoint(0);

          m_refGrad.resize(n);
          for (size_t a = 0; a < n; ++a)
          {
            m_refGrad[a].resize(d);
            const auto& basisFn = scalarFE.getBasis(a);
            for (size_t j = 0; j < d; ++j)
              m_refGrad[a](j) = basisFn.template getDerivative<1>(j)(rc);
          }

          m_basis.resize(n);
          for (size_t b = 0; b < n; ++b)
            m_basis[b] = scalarFE.getBasis(b)(rc);
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        const size_t n = m_refGrad.size();
        const size_t vdim = trialfes.getVectorDimension();

        const size_t ntr = lhs.getDOFs(polytope);
        const size_t nte = rhs.getDOFs(polytope);
        m_matrix.resize(static_cast<Eigen::Index>(nte), static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        for (size_t qp = 0; qp < m_ps.size(); ++qp)
        {
          const auto& p = m_ps[qp];
          const ScalarType wdet =
            static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

          const auto& Jinv = p.getJacobianInverse();
          const auto fval = coeff.getValue(p);

          for (size_t a = 0; a < n; ++a)
          {
            const auto physGrad = Jinv.transpose() * m_refGrad[a];
            const ScalarType gradDotF = Math::dot(physGrad, fval);

            for (size_t b = 0; b < n; ++b)
            {
              const ScalarType val = wdet * gradDotF * m_basis[b];

              for (size_t c = 0; c < vdim; ++c)
              {
                const size_t row = b * vdim + c;
                const size_t col = a * vdim + c;
                m_matrix(
                  static_cast<Eigen::Index>(row),
                  static_cast<Eigen::Index>(col)) += val;
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

      std::vector<Math::SpatialVector<ScalarType>> m_refGrad;
      std::vector<ScalarType> m_basis;

      Math::Matrix<ScalarType> m_matrix;
"""
replacements.append((2871, 3057, spec12))

# Sort replacements by start line descending to apply bottom-up
replacements.sort(key=lambda x: x[0], reverse=True)

for start, end, new_content in replacements:
    # Convert to 0-indexed
    s = start - 1
    e = end  # end is inclusive in 1-indexed, so lines[s:e] gives lines start..end
    new_lines = new_content.split('\n')
    # Ensure each line ends with newline
    new_lines = [line + '\n' for line in new_lines]
    # Remove last empty line if content ends with newline
    if new_lines and new_lines[-1].strip() == '':
        new_lines = new_lines[:-1]
    lines[s:e] = new_lines

with open(FILE, 'w') as f:
    f.writelines(lines)

print(f"Applied {len(replacements)} replacements successfully.")
print(f"New file has {len(lines)} lines.")
