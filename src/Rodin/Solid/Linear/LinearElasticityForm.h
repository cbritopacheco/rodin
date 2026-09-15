/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearElasticityForm.h
 * @brief Named bilinear form of isotropic linear elasticity.
 */
#ifndef RODIN_SOLID_LINEAR_LINEARELASTICITYFORM_H
#define RODIN_SOLID_LINEAR_LINEARELASTICITYFORM_H

#include <algorithm>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/Traits.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/RealFunction.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /// @brief Traits of a Variational::LinearElasticityForm.
  template <class Solution, class FES, class Operator, class LambdaDerived,
    class MuDerived>
  struct Traits<
    Variational::LinearElasticityForm<Solution, FES, Operator, LambdaDerived, MuDerived>>
  {
      /// @brief Solution type of the trial function.
      using SolutionType = Solution;

      /// @brief Trial finite element space type.
      using TrialFESType = FES;

      /// @brief Test finite element space type.
      using TestFESType = FES;

      /// @brief Assembled operator type.
      using OperatorType = Operator;
  };

  /// @brief Marks Variational::LinearElasticityForm as a named form.
  template <class Solution, class FES, class Operator, class LambdaDerived,
    class MuDerived>
  struct IsNamedForm<
    Variational::LinearElasticityForm<Solution, FES, Operator, LambdaDerived, MuDerived>>
    : std::true_type
  {
      /// @brief True for linear elasticity forms.
      static constexpr bool Value = true;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Discrete isotropic linear elasticity form
   * @f$ a(u,v) = \int_\Omega \lambda\, (\nabla \cdot u)(\nabla \cdot v)
   *   + 2 \mu\, \varepsilon(u) : \varepsilon(v)\,dx @f$.
   *
   * Computes the same operator as
   * <tt>LinearElasticityIntegral(u, v)(lambda, mu)</tt>, including the
   * quadrature order that integrator picks: the highest order either Lamé
   * parameter reports, counting an unknown order as zero, plus the orders of
   * the two bases. Vector-valued spaces whose dimension matches the mesh's.
   *
   * @tparam LambdaDerived Derived type of the first Lamé parameter.
   * @tparam MuDerived Derived type of the shear modulus.
   */
  template <class Solution, class FES, class Operator, class LambdaDerived,
    class MuDerived>
  class LinearElasticityForm final : public BilinearFormBase<Operator>
  {
      using MeshType = typename FormLanguage::Traits<FES>::MeshType;
      using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;

      static_assert(FormLanguage::IsVectorRange<
        typename FormLanguage::Traits<FES>::RangeType>::Value);

    public:
      /// @brief Solution type of the trial function.
      using SolutionType = Solution;

      /// @brief Finite element space type of both functions.
      using FESType = FES;

      /// @brief Assembled operator type.
      using OperatorType = Operator;

      /// @brief Scalar type of the operator entries.
      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      /// @brief First Lamé parameter type.
      using LambdaType = FunctionBase<LambdaDerived>;

      /// @brief Shear modulus type.
      using MuType = FunctionBase<MuDerived>;

      /**
       * @brief Local kernel of the linear elasticity form.
       *
       * Evaluates the reference derivatives of every basis function once per
       * polytope geometry and quadrature order and reuses them across
       * polytopes; at each polytope they are pulled back through @f$ J^{-1} @f$
       * and combined into divergences and symmetric gradients. The kernel
       * owns copies of the Lamé parameters, so that each thread assembling
       * with its own copy of the kernel evaluates its own copies.
       */
      class Kernel
      {
        public:
          /// @brief Local matrix type.
          using MatrixType = Math::Matrix<ScalarType>;

          /**
           * @brief Constructs the kernel.
           * @param[in] trialFES Trial space.
           * @param[in] testFES Test space.
           * @param[in] lambda First Lamé parameter, copied.
           * @param[in] mu Shear modulus, copied.
           */
          Kernel(const FES& trialFES, const FES& testFES, const LambdaType& lambda,
            const MuType& mu)
            : m_trialFES(trialFES),
              m_testFES(testFES),
              m_lambda(lambda.copy()),
              m_mu(mu.copy())
          {}

          /**
           * @brief Copies the kernel, including its own copies of the parameters.
           * @param[in] other Kernel to copy.
           */
          Kernel(const Kernel& other)
            : m_trialFES(other.m_trialFES),
              m_testFES(other.m_testFES),
              m_lambda(other.m_lambda->copy()),
              m_mu(other.m_mu->copy())
          {}

          /**
           * @brief Computes the local stiffness matrix of @p polytope.
           * @param[out] out Local matrix, test rows by trial columns.
           * @param[in] polytope Polytope to integrate over.
           */
          void compute(MatrixType& out, const Geometry::Polytope& polytope) const
          {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto geometry = polytope.getGeometry();
            const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
            const auto& testFE = m_testFES.get().getFiniteElement(d, i);
            const size_t trialCount = trialFE.getCount();
            const size_t testCount = testFE.getCount();
            const size_t order =
              std::max(m_lambda->getOrder(polytope).value_or(size_t(0)),
                m_mu->getOrder(polytope).value_or(size_t(0))) +
              trialFE.getOrder() + testFE.getOrder();
            const auto& qf = QF::PolytopeQuadratureFormula::get(order, geometry);
            const auto& quadrature = polytope.getQuadrature(qf);

            const bool rebuild = !m_cached || m_geometry != geometry ||
              m_order != order || m_trialCount != trialCount || m_testCount != testCount;
            if (rebuild)
            {
              m_cached = true;
              m_geometry = geometry;
              m_order = order;
              m_trialCount = trialCount;
              m_testCount = testCount;
              m_trialDerivatives.assign(qf.getSize(), {});
              m_testDerivatives.assign(qf.getSize(), {});
              for (size_t qp = 0; qp < qf.getSize(); ++qp)
              {
                const auto& reference = qf.getPoint(qp);
                m_trialDerivatives[qp] = getReferenceDerivatives(trialFE, reference, d);
                m_testDerivatives[qp] = getReferenceDerivatives(testFE, reference, d);
              }
            }

            out.resize(testCount, trialCount);
            out.setZero();
            m_trialSymmetric.resize(trialCount);
            m_trialDivergence.resize(trialCount);
            m_testSymmetric.resize(testCount);
            m_testDivergence.resize(testCount);
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const IntegrationPoint ip(point, &qf, qp);
              const ScalarType weight =
                static_cast<ScalarType>(qf.getWeight(qp) * point.getDistortion());
              const auto lambda = m_lambda->getValue(ip);
              const auto mu = m_mu->getValue(ip);
              const auto& jacobianInverse = point.getJacobianInverse();
              for (size_t tr = 0; tr < trialCount; ++tr)
              {
                const Math::SpatialMatrix<ScalarType> jacobian =
                  m_trialDerivatives[qp][tr] * jacobianInverse;
                m_trialSymmetric[tr] = jacobian + jacobian.adjoint();
                m_trialDivergence[tr] = jacobian.trace();
              }
              for (size_t te = 0; te < testCount; ++te)
              {
                const Math::SpatialMatrix<ScalarType> jacobian =
                  m_testDerivatives[qp][te] * jacobianInverse;
                m_testSymmetric[te] = jacobian + jacobian.adjoint();
                m_testDivergence[te] = jacobian.trace();
              }
              for (size_t te = 0; te < testCount; ++te)
              {
                for (size_t tr = 0; tr < trialCount; ++tr)
                {
                  out(te, tr) += weight *
                    (Math::dot(lambda * m_trialDivergence[tr], m_testDivergence[te]) +
                      static_cast<ScalarType>(0.5) *
                        Math::dot(mu * m_trialSymmetric[tr], m_testSymmetric[te]));
                }
              }
            }
          }

        private:
          /**
           * @brief Evaluates the reference Jacobian of every basis function of
           * @p fe at @p reference.
           */
          template <class FiniteElement, class Point>
          static std::vector<Math::SpatialMatrix<ScalarType>> getReferenceDerivatives(
            const FiniteElement& fe, const Point& reference, size_t d)
          {
            std::vector<Math::SpatialMatrix<ScalarType>> derivatives(fe.getCount());
            for (size_t b = 0; b < fe.getCount(); ++b)
            {
              const auto& basis = fe.getBasis(b);
              auto& derivative = derivatives[b];
              derivative.resize(d, d);
              for (size_t r = 0; r < d; ++r)
                for (size_t c = 0; c < d; ++c)
                  derivative(r, c) = basis.template getDerivative<1>(r, c)(reference);
            }
            return derivatives;
          }

          std::reference_wrapper<const FES> m_trialFES;
          std::reference_wrapper<const FES> m_testFES;
          std::unique_ptr<LambdaType> m_lambda;
          std::unique_ptr<MuType> m_mu;
          mutable bool m_cached = false;
          mutable Geometry::Polytope::Type m_geometry = Geometry::Polytope::Type::Point;
          mutable size_t m_order = 0;
          mutable size_t m_trialCount = 0;
          mutable size_t m_testCount = 0;
          mutable std::vector<std::vector<Math::SpatialMatrix<ScalarType>>>
            m_trialDerivatives;
          mutable std::vector<std::vector<Math::SpatialMatrix<ScalarType>>>
            m_testDerivatives;
          mutable std::vector<Math::SpatialMatrix<ScalarType>> m_trialSymmetric;
          mutable std::vector<ScalarType> m_trialDivergence;
          mutable std::vector<Math::SpatialMatrix<ScalarType>> m_testSymmetric;
          mutable std::vector<ScalarType> m_testDivergence;
      };

      /// @brief Local kernel type, as the assembly asks for it.
      using KernelType = Kernel;

      /// @brief Parent class type.
      using Parent = BilinearFormBase<OperatorType>;

      /// @brief Assembly type selected for the space's context.
      using AssemblyType = typename Assembly::Default<ContextType,
        ContextType>::template Type<OperatorType, LinearElasticityForm>;

      /**
       * @brief Constructs and assembles the linear elasticity form.
       *
       * Each Lamé parameter is either a function, copied, or a value lifted to
       * a RealFunction.
       *
       * @param[in] lambda First Lamé parameter.
       * @param[in] mu Shear modulus.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      template <class L, class M>
      LinearElasticityForm(const L& lambda, const M& mu,
        const TrialFunction<SolutionType, FESType>& u, const TestFunction<FESType>& v)
        : m_u(u),
          m_v(v),
          m_lambda(lift<LambdaDerived>(lambda)),
          m_mu(lift<MuDerived>(mu))
      {
        assemble();
      }

      /**
       * @brief Copy constructor, copying the Lamé parameters.
       * @param[in] other Form to copy.
       */
      LinearElasticityForm(const LinearElasticityForm& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v),
          m_lambda(other.m_lambda->copy()),
          m_mu(other.m_mu->copy()),
          m_attributes(other.m_attributes),
          m_operator(other.m_operator),
          m_assembly(other.m_assembly)
      {}

      /// @brief Move constructor.
      LinearElasticityForm(LinearElasticityForm&&) = default;

      OperatorType& getOperator() override
      {
        return m_operator;
      }

      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      void assemble() override
      {
        m_assembly.execute(m_operator, *this);
      }

      const TrialFunction<SolutionType, FESType>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<FESType>& getTestFunction() const override
      {
        return m_v.get();
      }

      /// @brief Gets the first Lamé parameter.
      const LambdaType& getLameFirstParameter() const
      {
        return *m_lambda;
      }

      /// @brief Gets the shear modulus.
      const MuType& getShearModulus() const
      {
        return *m_mu;
      }

      /// @brief Gets the region the form integrates over.
      Geometry::Region getRegion() const
      {
        return Geometry::Region::Cells;
      }

      /// @brief Gets the attributes the form is restricted to, empty for all.
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_attributes;
      }

      /**
       * @brief Restricts the form to cells carrying one of @p attributes and
       * reassembles it.
       * @param[in] attributes Non-empty set of attributes.
       * @returns Reference to this form.
       */
      LinearElasticityForm& over(const FlatSet<Geometry::Attribute>& attributes)
      {
        assert(attributes.size() > 0);
        m_attributes = attributes;
        assemble();
        return *this;
      }

      /**
       * @brief Restricts the form to cells carrying the given attributes and
       * reassembles it.
       * @param[in] a1 First attribute.
       * @param[in] as Further attributes.
       * @returns Reference to this form.
       */
      template <class A1, class... As>
      LinearElasticityForm& over(const A1& a1, const As&... as)
      {
        return over(FlatSet<Geometry::Attribute>{a1, as...});
      }

      /**
       * @brief Builds the kernel the assembly runs on each polytope.
       * @returns Kernel of this form.
       */
      Kernel getKernel() const
      {
        return Kernel(getTrialFunction().getFiniteElementSpace(),
          getTestFunction().getFiniteElementSpace(), *m_lambda, *m_mu);
      }

      LinearElasticityForm* copy() const noexcept override
      {
        return new LinearElasticityForm(*this);
      }

    private:
      /**
       * @brief Copies @p value as a Lamé parameter, lifting a non-function to
       * a RealFunction.
       */
      template <class Derived, class Value>
      static std::unique_ptr<FunctionBase<Derived>> lift(const Value& value)
      {
        if constexpr (std::is_base_of_v<FunctionBase<Derived>, Value>)
        {
          return std::unique_ptr<FunctionBase<Derived>>(
            static_cast<const FunctionBase<Derived>&>(value).copy());
        }
        else
        {
          static_assert(std::is_same_v<Derived,
            typename FormLanguage::FunctionDerived<RealFunction<Value>>::Type>);
          return std::unique_ptr<FunctionBase<Derived>>(
            static_cast<FunctionBase<Derived>*>(RealFunction<Value>(value).copy()));
        }
      }

      std::reference_wrapper<const TrialFunction<SolutionType, FESType>> m_u;
      std::reference_wrapper<const TestFunction<FESType>> m_v;
      std::unique_ptr<LambdaType> m_lambda;
      std::unique_ptr<MuType> m_mu;
      FlatSet<Geometry::Attribute> m_attributes;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  /// @brief Deduction guide for two function Lamé parameters.
  template <class LambdaDerived, class MuDerived, class Solution, class FES>
  LinearElasticityForm(const FunctionBase<LambdaDerived>&, const FunctionBase<MuDerived>&,
    const TrialFunction<Solution, FES>&, const TestFunction<FES>&)
    -> LinearElasticityForm<Solution, FES,
      Math::SparseMatrix<typename FormLanguage::Traits<FES>::ScalarType>, LambdaDerived,
      MuDerived>;

  /// @brief Deduction guide for a function first parameter and a lifted shear modulus.
  template <class LambdaDerived, class M, class Solution, class FES>
    requires(!std::is_base_of_v<FormLanguage::Base, M>)
  LinearElasticityForm(const FunctionBase<LambdaDerived>&, const M&,
    const TrialFunction<Solution, FES>&, const TestFunction<FES>&)
    -> LinearElasticityForm<Solution, FES,
      Math::SparseMatrix<typename FormLanguage::Traits<FES>::ScalarType>, LambdaDerived,
      typename FormLanguage::FunctionDerived<RealFunction<M>>::Type>;

  /// @brief Deduction guide for a lifted first parameter and a function shear modulus.
  template <class L, class MuDerived, class Solution, class FES>
    requires(!std::is_base_of_v<FormLanguage::Base, L>)
  LinearElasticityForm(const L&, const FunctionBase<MuDerived>&,
    const TrialFunction<Solution, FES>&, const TestFunction<FES>&)
    -> LinearElasticityForm<Solution, FES,
      Math::SparseMatrix<typename FormLanguage::Traits<FES>::ScalarType>,
      typename FormLanguage::FunctionDerived<RealFunction<L>>::Type, MuDerived>;

  /// @brief Deduction guide for two lifted Lamé parameters.
  template <class L, class M, class Solution, class FES>
    requires(!std::is_base_of_v<FormLanguage::Base, L> &&
              !std::is_base_of_v<FormLanguage::Base, M>)
  LinearElasticityForm(
    const L&, const M&, const TrialFunction<Solution, FES>&, const TestFunction<FES>&)
    -> LinearElasticityForm<Solution, FES,
      Math::SparseMatrix<typename FormLanguage::Traits<FES>::ScalarType>,
      typename FormLanguage::FunctionDerived<RealFunction<L>>::Type,
      typename FormLanguage::FunctionDerived<RealFunction<M>>::Type>;
}

#endif
