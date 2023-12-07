/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_POTENTIAL_H
#define RODIN_VARIATIONAL_POTENTIAL_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/QF/QuadratureFormula.h"

#include "ForwardDecls.h"
#include "Mult.h"
#include "Dot.h"
#include "Function.h"
#include "Integral.h"
#include "ShapeFunction.h"
#include "BilinearForm.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHSType, class RHSDerived>
  struct Traits<
    Variational::Potential<
      LHSType,
      Variational::FunctionBase<RHSDerived>>>
  {
    using LHS = LHSType;

    using RHS = Variational::FunctionBase<Variational::FunctionBase<RHSDerived>>;

    using Kernel = LHS;

    using Operand = RHS;

    using RHSRange =
      typename FormLanguage::Traits<RHS>::RangeType;

    using LHSRange =
      std::conditional_t<
        // If
        std::is_same_v<RHSRange, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Math::Vector>,
          // Then
          Math::Matrix,
          // Else
          void>>;

    using RangeType = RHSRange;
  };

  /**
   * @ingroup TraitsSpecializations
   * @brief Traits for Potential.
   */
  template <class LHSType, class RHSDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Potential<
      LHSType,
      Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FESType, SpaceType>>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;

    using LHS = LHSType;

    using RHS = Variational::ShapeFunctionBase<Variational::ShapeFunction<RHSDerived, FES, Space>, FES, Space>;

    using Kernel = LHS;

    using Operand = RHS;

    using RHSRange =
      typename FormLanguage::Traits<RHS>::RangeType;

    using LHSRange =
      std::conditional_t<
        // If
        std::is_same_v<RHSRange, Scalar>,
        // Then
        Scalar,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Math::Vector>,
          // Then
          Math::Matrix,
          // Else
          void>>;

    using RangeType = RHSRange;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup PotentialSpecializations Potential Template Specializations
   * @brief Template specializations of the Potential class.
   * @see Potential
   */

  /**
   * @ingroup PotentialSpecializations
   */
  template <class LHSType, class RHSDerived>
  class Potential<LHSType, FunctionBase<RHSDerived>> final
    : public FunctionBase<Potential<LHSType, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = LHSType;

      using Kernel = LHS;

      using RHS = FunctionBase<RHSDerived>;

      using Operand = RHS;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      using LHSRange =
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Scalar>,
          // Then
          Scalar,
          // Else
          std::conditional_t<
            // If
            std::is_same_v<RHSRange, Math::Vector>,
            // Then
            Math::Matrix,
            // Else
            void>>;

      using Parent = FunctionBase<Potential<LHS, RHS>>;

      Potential(const Kernel& kernel, const Operand& u)
        : m_kernel(kernel), m_u(u.copy())
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel),
          m_u(other.m_u->copy())
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)),
          m_u(std::move(other.m_u))
      {}

      inline
      const auto& getKernel() const
      {
        return m_kernel.get();
      }

      inline
      const auto& getOperand() const
      {
        assert(m_u);
        return *m_u;
      }

      inline
      auto getValue(const Geometry::Point& p) const
      {
        const auto& kernel = getKernel();
        const auto& operand = getOperand();
        const auto& mesh = p.getPolytope().getMesh();
        if (m_qf.has_value())
        {
          if constexpr (std::is_same_v<RHSRange, Scalar>)
          {
            Scalar res = 0;
            for (auto it = mesh.getCell(); it; ++it)
            {
              const auto& polytope = *it;
              const auto& qf = m_qf.value()(polytope);
              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Geometry::Point y(polytope, trans, std::ref(qf.getPoint(i)));
                res += qf.getWeight(i) * y.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRange, Math::Vector>)
          {
            Math::Vector res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
        else
        {
          if constexpr (std::is_same_v<RHSRange, Scalar>)
          {
            Scalar res = 0;
            for (auto it = mesh.getCell(); it; ++it)
            {
              const auto& polytope = *it;
              const QF::GenericPolytopeQuadrature qf(polytope.getGeometry());

              const auto& trans = polytope.getTransformation();
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Geometry::Point y(polytope, trans, std::cref(qf.getPoint(i)));
                res += qf.getWeight(i) * y.getDistortion() * kernel(p, y) * operand(y);
              }
            }
            return res;
          }
          else if constexpr (std::is_same_v<RHSRange, Math::Vector>)
          {
            Math::Vector res;
            getValue(res, p);
            return res;
          }
          else
          {
            assert(false);
            return void();
          }
        }
      }

      inline
      constexpr
      void getValue(Math::Vector& res, const Geometry::Point& p) const
      {
        assert(false);
        res.setConstant(NAN);
      }

      Potential& setQuadratureFormula(
          const std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>& qf)
      {
        m_qf.emplace(qf);
        return *this;
      }

      inline
      const auto& getQuadratureFormula() const
      {
        return m_qf;
      }

    private:
      std::reference_wrapper<const Kernel> m_kernel;
      std::unique_ptr<Operand> m_u;
      std::optional<
        std::function<const QF::QuadratureFormulaBase&(const Geometry::Polytope&)>> m_qf;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived>
  Potential(const LHSType&, const FunctionBase<RHSDerived>&)
    -> Potential<LHSType, FunctionBase<RHSDerived>>;

  /**
   * @ingroup PotentialSpecializations
   *
   * Represents the expression:
   * @f[
   *  (K u)(x) = \sum_{\tau \in \mathcal{T}_h} (K_\tau u)(x)
   *  = \sum_{\tau \in \mathcal{T}_h} \sum^{n(\tau)}_{\ell = 1} w_{\tau, \ell}
   *  (K_\tau \phi_{\tau, \ell})(x), \quad w_{\tau, \ell} \in \mathbb{R}
   *  @f]
   *  where:
   *  @f[
   *    (K_\tau \phi_{\tau, \ell})(x) := \ \mathrm{p.v.} \int_{\tau} k(x, y) \phi_{\tau, \ell}(y) \ dy
   *  @f]
   *  are the basis functions.
   */
  template <class LHSType, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Potential<
    LHSType,
    ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>> final
    : public LinearFormBase<Math::Vector>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using LHS = LHSType;

      using Kernel = LHS;

      using RHS = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, SpaceType>, FES, Space>;

      using Operand = RHS;

      using Parent = LinearFormBase<Math::Vector>;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      using LHSRange =
        std::conditional_t<
          // If
          std::is_same_v<RHSRange, Scalar>,
          // Then
          Scalar,
          // Else
          std::conditional_t<
            // If
            std::is_same_v<RHSRange, Math::Vector>,
            // Then
            Math::Matrix,
            // Else
            void>>;

      static_assert(
          (std::is_same_v<LHSRange, Scalar> && std::is_same_v<RHSRange, Scalar>) ||
          (std::is_same_v<LHSRange, Math::Matrix> || std::is_same_v<RHSRange, Math::Vector>));

      Potential(const Kernel& kernel, const Operand& u)
        : m_kernel(kernel), m_u(u)
      {}

      Potential(const Potential& other)
        : Parent(other),
          m_kernel(other.m_kernel), m_u(other.m_u)
      {}

      Potential(Potential&& other)
        : Parent(std::move(other)),
          m_kernel(std::move(other.m_kernel)), m_u(std::move(other.m_u))
      {}

      inline
      const LHS& getLHS() const
      {
        return getKernel();
      }

      inline
      const RHS& getRHS() const
      {
        return getOperand();
      }

      inline
      const Kernel& getKernel() const
      {
        return m_kernel;
      }

      inline
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          static_assert(std::is_same_v<RHSRange, Scalar>);
          return { 1, 1 };
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Matrix>)
        {
          static_assert(std::is_same_v<RHSRange, Math::Vector>);
          return getRHS().getRangeShape()[0];
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      void assemble() final override
      {
        assert(false);
      }

      Math::Vector& getVector() override
      {
        return m_vec;
      }

      const Math::Vector& getVector() const override
      {
        return m_vec;
      }

      inline
      const Operand& getTestFunction() const override
      {
        return getOperand();
      }

      inline Potential* copy() const noexcept override
      {
        return new Potential(*this);
      }

    private:
      std::reference_wrapper<Kernel> m_kernel;
      std::reference_wrapper<const Operand> m_u;
      Math::Vector m_vec;
  };

  /**
   * @brief CTAD for Potential.
   */
  template <class LHSType, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  Potential(const LHSType&, const ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>&)
    -> Potential<LHSType, ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, SpaceType>, FESType, SpaceType>>;

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Integral<
    Dot<
      Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
        : public BilinearFormBase<Math::Matrix>
  {
    public:
      static void add(
          Math::Matrix& out, const Math::Matrix& in,
          const IndexArray& rows, const IndexArray& cols)
      {
        assert(rows.size() >= 0);
        assert(cols.size() >= 0);
        assert(in.rows() == rows.size());
        assert(in.cols() == cols.size());
        for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
        {
          for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
            out(rows(i), cols(j)) += in(i, j);
        }
      }

      using Kernel = KernelType;

      using LHS = Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>;

      using RHS =
        ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      using Integrand = Dot<LHS, RHS>;

      using Parent = BilinearFormBase<Math::Matrix>;

      static_assert(std::is_same_v<LHSRange, RHSRange>);

      Integral(const LHS& lhs, const RHS& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      Integral(const Integrand& prod)
        : m_integrand(prod.copy())
      {}

      Integral(const Integral& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      Integral(Integral&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      void assemble() final override
      {
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& mesh = trialfes.getMesh();
        const auto& kernel = lhs.getKernel();
        auto& res = m_mat;
        res.resize(testfes.getSize(), trialfes.getSize());
        res.setZero();
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          for (auto itt = mesh.getCell(); itt; ++itt)
          {
            const auto& t = *itt;
            for (auto itau = mesh.getCell(); itau; ++itau)
            {
              const auto& tau = *itau;
              integrate(res, t, tau);
            }
          }
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Vector>)
        {
          assert(false);
        }
        else
        {
          assert(false);
        }
      }

      virtual void integrate(Math::Matrix& out, const Geometry::Polytope& trgeom, const Geometry::Polytope& tegeom)
      {
        const auto& tau = trgeom;
        const auto& t = tegeom;
        const auto& tautrans = tau.getTransformation();
        const auto& ttrans = t.getTransformation();
        const auto& integrand = getIntegrand();
        const auto& lhs = integrand.getLHS();
        const auto& rhs = integrand.getRHS();
        const auto& trialfes = lhs.getFiniteElementSpace();
        const auto& testfes = rhs.getFiniteElementSpace();
        const auto& testfe = testfes.getFiniteElement(t.getDimension(), t.getIndex());
        const auto& trialfe = trialfes.getFiniteElement(tau.getDimension(), tau.getIndex());
        const auto& kernel = lhs.getKernel();
        if (t == tau)
        {
          switch (t.getGeometry())
          {
            case Geometry::Polytope::Type::Point:
            case Geometry::Polytope::Type::Segment:
            case Geometry::Polytope::Type::Tetrahedron:
            {
              assert(false);
              break;
            }
            case Geometry::Polytope::Type::Quadrilateral:
            {
              assert(false);
              break;
            }
            case Geometry::Polytope::Type::Triangle:
            {
              const size_t order = std::max(trialfe.getOrder(), testfe.getOrder());
              const QF::GenericPolytopeQuadrature qf(order, Geometry::Polytope::Type::Segment);
              Math::SpatialVector rx1(2), rz1(2),
                                  rx2(2), rz2(2),
                                  rx3(2), rz3(2),
                                  rx4(2), rz4(2),
                                  rx5(2), rz5(2),
                                  rx6(2), rz6(2);
              for (size_t i = 0; i < qf.getSize(); i++)
              {
                const Scalar xi = 1 - qf.getPoint(i).value();
                for (size_t j = 0; j < qf.getSize(); j++)
                {
                  const Scalar eta1 = qf.getPoint(j).value();
                  for (size_t h = 0; h < qf.getSize(); h++)
                  {
                    const Scalar eta2 = qf.getPoint(h).value();
                    rx1 << xi, xi * (1 - eta1 + eta1 * eta2);
                    rz2 << xi, xi * (1 - eta1 + eta1 * eta2);
                    rz3 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    rx4 << xi * (1 - eta1 * eta2), xi * eta1 * (1 - eta2);
                    rz5 << xi, xi * eta1 * (1 - eta2);
                    rx6 << xi, xi * eta1 * (1 - eta2);
                    const Geometry::Point x1(t, ttrans, std::cref(rx1));
                    const Geometry::Point z2(t, ttrans, std::cref(rz2));
                    const Geometry::Point z3(t, ttrans, std::cref(rz3));
                    const Geometry::Point x4(t, ttrans, std::cref(rx4));
                    const Geometry::Point z5(t, ttrans, std::cref(rz5));
                    const Geometry::Point x6(t, ttrans, std::cref(rx6));
                    const Scalar d = xi * xi * xi * eta1 * eta1 * eta2;
                    for (size_t k = 0; k < qf.getSize(); k++)
                    {
                      const Scalar w = qf.getWeight(i) * qf.getWeight(j) * qf.getWeight(h) * qf.getWeight(k);
                      const Scalar eta3 = qf.getPoint(k).value();
                      rz1 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      rx2 << xi * (1 - eta1 * eta2 * eta3), xi * (1 - eta1);
                      rx3 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      rz4 << xi, xi * eta1 * (1 - eta2 + eta2 * eta3);
                      rx5 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      rz6 << xi * (1 - eta1 * eta2 * eta3), xi * eta1 * (1 - eta2 * eta3);
                      const Geometry::Point z1(t, ttrans, std::cref(rz1));
                      const Geometry::Point x2(t, ttrans, std::cref(rx2));
                      const Geometry::Point x3(t, ttrans, std::cref(rx3));
                      const Geometry::Point z4(t, ttrans, std::cref(rz4));
                      const Geometry::Point x5(t, ttrans, std::cref(rx5));
                      const Geometry::Point z6(t, ttrans, std::cref(rz6));
                      const Scalar s1 = kernel(x1, z1) * x1.getDistortion() * z1.getDistortion();
                      const Scalar s2 = kernel(x2, z2) * x2.getDistortion() * z2.getDistortion();
                      const Scalar s3 = kernel(x3, z3) * x3.getDistortion() * z3.getDistortion();
                      const Scalar s4 = kernel(x4, z4) * x4.getDistortion() * z4.getDistortion();
                      const Scalar s5 = kernel(x5, z5) * x5.getDistortion() * z5.getDistortion();
                      const Scalar s6 = kernel(x6, z6) * x6.getDistortion() * z6.getDistortion();
                      for (size_t l = 0; l < testfe.getCount(); l++)
                      {
                        const auto& teb = testfe.getBasis(l);
                        const Index teglobal = testfes.getGlobalIndex({ t.getDimension(), t.getIndex() }, l);
                        for (size_t m = 0; m < trialfe.getCount(); m++)
                        {
                          const auto& trb = trialfe.getBasis(m);
                          const Index trglobal = trialfes.getGlobalIndex({ tau.getDimension(), tau.getIndex() }, m);
                          out(teglobal, trglobal) += d * w * s1 * trb(rx1) * teb(rz1);
                          assert(std::isfinite(s1));
                          out(teglobal, trglobal) += d * w * s2 * trb(rx2) * teb(rz2);
                          assert(std::isfinite(s2));
                          out(teglobal, trglobal) += d * w * s3 * trb(rx3) * teb(rz3);
                          assert(std::isfinite(s3));
                          out(teglobal, trglobal) += d * w * s4 * trb(rx4) * teb(rz4);
                          assert(std::isfinite(s4));
                          out(teglobal, trglobal) += d * w * s5 * trb(rx5) * teb(rz5);
                          assert(std::isfinite(s5));
                          out(teglobal, trglobal) += d * w * s6 * trb(rx6) * teb(rz6);
                          assert(std::isfinite(s6));
                        }
                      }
                    }
                  }
                }
              }
              break;
            }
          }
        }
        else
        {
          const QF::GenericPolytopeQuadrature qftr(trialfe.getOrder(), tau.getGeometry());
          const QF::GenericPolytopeQuadrature qfte(testfe.getOrder(), t.getGeometry());
          for (size_t i = 0; i < qfte.getSize(); i++)
          {
            const Geometry::Point x(t, ttrans, std::cref(qfte.getPoint(i)));
            for (size_t j = 0; j < qftr.getSize(); j++)
            {
              const Geometry::Point y(tau, tautrans, std::cref(qftr.getPoint(j)));
              const Scalar d = x.getDistortion() * y.getDistortion();
              const Scalar kxy = kernel(x, y);
              const Scalar w = qfte.getWeight(i) * qftr.getWeight(j);
              for (size_t l = 0; l < testfe.getCount(); l++)
              {
                const auto& teb = testfe.getBasis(l);
                const Index teglobal = testfes.getGlobalIndex({ t.getDimension(), t.getIndex() }, l);
                for (size_t m = 0; m < trialfe.getCount(); m++)
                {
                  const auto& trb = trialfe.getBasis(m);
                  const Index trglobal = testfes.getGlobalIndex({ tau.getDimension(), tau.getIndex() }, m);
                  out(teglobal, trglobal) += w * d * kxy * teb(qfte.getPoint(i)) * trb(qftr.getPoint(j));
                }
              }
            }
          }
        }
      }

      Math::Matrix& getOperator() override
      {
        return m_mat;
      }

      const Math::Matrix& getOperator() const override
      {
        return m_mat;
      }

      const TrialFunction<TrialFES>& getTrialFunction() const override
      {
        return getIntegrand().getLHS().getLeaf();
      }

      const TestFunction<TestFES>& getTestFunction() const override
      {
        return getIntegrand().getRHS().getLeaf();
      }

      inline Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      std::unique_ptr<Integrand> m_integrand;
      Math::Matrix m_mat;
  };

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const
      Dot<
        Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
        ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>&,
        const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> Integral<
          Dot<
            Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
            ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;
}

#endif

