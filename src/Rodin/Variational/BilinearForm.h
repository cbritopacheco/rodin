/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BilinearForm.h
 * @brief Bilinear form classes for finite element assembly.
 *
 * This file defines the BilinearForm classes which represent bilinear forms
 * @f$ a(u,v) @f$ in variational formulations. Bilinear forms are assembled into
 * matrices and form the left-hand side of finite element systems.
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include "Rodin/Math/Traits.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Assembly/ForwardDecls.h"

#include "Exceptions/TrialFunctionMismatchException.h"
#include "Exceptions/TestFunctionMismatchException.h"

#include "ForwardDecls.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  /// @brief Traits for bilinear form bases.
  template <class Operator>
  struct Traits<Variational::BilinearFormBase<Operator>>
  {
    /// @brief Assembled operator type.
      using OperatorType = Operator;
  };

  /// @brief Traits for bilinear forms.
  template <class Solution, class TrialFES, class TestFES, class Operator>
  struct Traits<Variational::BilinearForm<Solution, TrialFES, TestFES, Operator>>
  {
    /// @brief Solution vector type.
      using SolutionType = Solution;
    /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;
    /// @brief Test finite element space type.
      using TestFESType = TestFES;
    /// @brief Assembled operator type.
      using OperatorType = Operator;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup BilinearFormSpecializations BilinearForm Template Specializations
   * @brief Template specializations of the BilinearForm class.
   * @see BilinearForm
   */

  /**
   * @ingroup RodinVariational
   * @brief Base class for bilinear form representations.
   *
   * BilinearFormBase provides the foundation for representing bilinear forms
   * @f$ a(u,v) : V \times V \to \mathbb{R} @f$ in finite element computations.
   * A bilinear form is a function that is linear in both arguments and forms
   * the basis for defining variational problems.
   *
   * @tparam Operator Matrix type for the discrete representation
   *
   * ## Mathematical Foundation
   * A bilinear form @f$ a(u,v) @f$ satisfies:
   * - **Linearity in first argument**: @f$ a(\alpha u_1 + \beta u_2, v) = \alpha a(u_1,v) + \beta a(u_2,v) @f$
   * - **Linearity in second argument**: @f$ a(u, \alpha v_1 + \beta v_2) = \alpha a(u,v_1) + \beta a(u,v_2) @f$
   *
   * ## Discrete Representation  
   * The discrete matrix representation satisfies @f$ A_{ij} = a(\phi_j, \psi_i) @f$
   * where @f$ \phi_j @f$ are trial basis functions and @f$ \psi_i @f$ are test basis functions.
   */
  template <class Operator>
  class BilinearFormBase : public FormLanguage::Base
  {
    public:
      /// @brief Matrix operator type for the discrete representation
      using OperatorType = Operator;

      /// @brief Scalar type for matrix entries
      using ScalarType =
        typename FormLanguage::Traits<OperatorType>::ScalarType;

      /// @brief Parent class type.
      using Parent =
        FormLanguage::Base;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType =
        LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType =
        GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief List of local bilinear form integrators.
      using LocalBilinearFormIntegratorBaseListType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief List of global bilinear form integrators.
      using GlobalBilinearFormIntegratorBaseListType =
        FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /**
       * @brief Default constructor creating an empty bilinear form.
       *
       * Constructs a bilinear form with no integrators. Integrators must be
       * added separately to define the bilinear form @f$ a(u,v) @f$.
       */
      constexpr
      BilinearFormBase() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other Bilinear form to copy
       *
       * Creates a copy including all local and global integrators.
       */
      constexpr
      BilinearFormBase(const BilinearFormBase& other)
        : Parent(other),
          m_lbfis(other.m_lbfis),
          m_gbfis(other.m_gbfis)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Bilinear form to move from
       */
      constexpr
      BilinearFormBase(BilinearFormBase&& other)
        : Parent(std::move(other)),
          m_lbfis(std::move(other.m_lbfis)),
          m_gbfis(std::move(other.m_gbfis))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Bilinear form to copy
       * @returns Reference to this bilinear form
       */
      constexpr
      BilinearFormBase& operator=(const BilinearFormBase& other)
      {
        if (this != &other)
        {
          m_lbfis = other.m_lbfis;
          m_gbfis = other.m_gbfis;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Bilinear form to move from
       * @returns Reference to this bilinear form
       */
      constexpr
      BilinearFormBase& operator=(BilinearFormBase&& other) noexcept
      {
        if (this != &other)
        {
          m_lbfis = std::move(other.m_lbfis);
          m_gbfis = std::move(other.m_gbfis);
        }
        return *this;
      }

      /**
       * @brief Gets the list of local integrators.
       * @return Reference to the list of local element integrators
       *
       * Local integrators operate element-by-element, contributing to matrix
       * entries through local assembly:
       * @f$ A^K_{ij} = \int_K \text{integrand}(\phi_j, \psi_i) \, dx @f$
       */
      constexpr
      LocalBilinearFormIntegratorBaseListType& getLocalIntegrators()
      {
        return m_lbfis;
      }

      /**
       * @brief Gets the list of local integrators (const version).
       * @return Const reference to the list of local element integrators
       */
      constexpr
      const LocalBilinearFormIntegratorBaseListType& getLocalIntegrators() const
      {
        return m_lbfis;
      }

      /**
       * @brief Gets the list of global integrators.
       * @return Reference to the list of global integrators
       *
       * Global integrators operate on the entire global matrix, typically used
       * for non-local operators or constraint enforcement.
       */
      constexpr
      GlobalBilinearFormIntegratorBaseListType& getGlobalIntegrators()
      {
        return m_gbfis;
      }

      /**
       * @brief Gets the list of global integrators (const version).
       * @return Const reference to the list of global integrators
       */
      constexpr
      const GlobalBilinearFormIntegratorBaseListType& getGlobalIntegrators() const
      {
        return m_gbfis;
      }

      /**
       * @brief Builds the bilinear form the given bilinear integrator
       * @param[in] bfi Bilinear integrator which will be used to
       * build the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearFormBase& operator=(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        m_lbfis.clear();
        m_lbfis.add(bfi);
        return *this;
      }

      /**
       * @brief Replaces the local integrator list.
       * @param[in] bfi List of local bilinear integrators
       * @returns Reference to this bilinear form
       */
      BilinearFormBase& operator=(const LocalBilinearFormIntegratorBaseListType& bfi)
      {
        m_lbfis.clear();
        m_lbfis.add(bfi);
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearFormBase& operator+=(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_lbfis.add(bfi);
        return *this;
      }

      /**
       * @brief Adds local bilinear integrators to the bilinear form.
       * @param[in] bfis Local bilinear integrators to add
       * @returns Reference to this bilinear form
       */
      BilinearFormBase& operator+=(const LocalBilinearFormIntegratorBaseListType& bfis)
      {
        m_lbfis.add(bfis);
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearFormBase& operator+=(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_gbfis.add(bfi);
        return *this;
      }

      /**
       * @brief Adds global bilinear integrators to the bilinear form.
       * @param[in] bfis Global bilinear integrators to add
       * @returns Reference to this bilinear form
       */
      BilinearFormBase& operator+=(const GlobalBilinearFormIntegratorBaseListType& bfis)
      {
        m_gbfis.add(bfis);
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearFormBase& operator-=(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_lbfis.add(UnaryMinus(bfi));
        return *this;
      }

      /**
       * @brief Subtracts local bilinear integrators from the bilinear form.
       * @param[in] bfis Local bilinear integrators to subtract
       * @returns Reference to this bilinear form
       */
      BilinearFormBase& operator-=(const LocalBilinearFormIntegratorBaseListType& bfis)
      {
        m_lbfis.add(UnaryMinus(bfis));
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearFormBase& operator-=(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_gbfis.add(UnaryMinus(bfi));
        return *this;
      }

      /**
       * @brief Subtracts global bilinear integrators from the bilinear form.
       * @param[in] bfis Global bilinear integrators to subtract
       * @returns Reference to this bilinear form
       */
      BilinearFormBase& operator-=(const GlobalBilinearFormIntegratorBaseListType& bfis)
      {
        m_gbfis.add(UnaryMinus(bfis));
        return *this;
      }

      /**
       * @brief Gets the reference to the associated operator of the bilinear
       * form.
       */
      virtual OperatorType& getOperator() = 0;

      /** @brief Gets a constant reference to the associated operator of the
       * bilinear form.
       */
      virtual const OperatorType& getOperator() const = 0;

      /**
       * @brief Assembles the bilinear form.
       *
       * This method will assemble the underlying sparse matrix associated
       * the bilinear form.
       *
       * @see getMatrix()
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the reference to the associated TrialFunction object.
       * @returns Reference to this (for method chaining)
       */
      virtual const FormLanguage::Base& getTrialFunction() const = 0;

      /**
       * @brief Gets the reference to the associated TestFunction object.
       * @returns Reference to this (for method chaining)
       */
      virtual const FormLanguage::Base& getTestFunction() const = 0;

      virtual BilinearFormBase* copy() const noexcept override = 0;

    private:
      LocalBilinearFormIntegratorBaseListType m_lbfis;
      GlobalBilinearFormIntegratorBaseListType m_gbfis;
  };

  /**
   * @ingroup BilinearFormSpecializations
   * @brief Speciallization of BilinearForm for a matrix type.
   *
   * This specialization aids in the construction of a @f$ m \times n @f$
   * matrix @f$ A @f$, which is associated to the bilinear form. Here, @f$ n
   * @f$ represents the size (total number of degrees-of-freedom) of the trial
   * space, and @f$ m @f$ represents the size of the test space.
   */
  template <class Solution, class TrialFES, class TestFES, class Scalar>
  class BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Scalar>> final
    : public BilinearFormBase<Math::SparseMatrix<Scalar>>
  {
    using TrialFESMeshType =
      typename FormLanguage::Traits<TrialFES>::MeshType;

    using TestFESMeshType =
      typename FormLanguage::Traits<TestFES>::MeshType;

    using TrialFESContextType =
      typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

    using TestFESContextType =
      typename FormLanguage::Traits<TestFESMeshType>::ContextType;

    public:
      /// @brief Solution vector type.
      using SolutionType =
        Solution;

      /// @brief Scalar value type.
      using ScalarType =
        Scalar;

      /// @brief Type of operator associated to the bilinear form.
      using OperatorType =
        Math::SparseMatrix<ScalarType>;

      /// @brief Default assembly backend selected from the trial/test contexts.
      using DefaultAssemblyType =
        typename Assembly::Default<TrialFESContextType, TestFESContextType>
          ::template Type<OperatorType, BilinearForm>;

      /// @brief Assembly backend used by this bilinear form.
      using AssemblyType =
        DefaultAssemblyType;

      /// @brief Parent class type.
      using Parent = BilinearFormBase<OperatorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a sparse bilinear form from trial and test functions.
       * @param[in] u Trial function
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Bilinear form to copy
       */
      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_operator(other.m_operator),
          m_assembly(other.m_assembly)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Bilinear form to move from
       */
      constexpr
      BilinearForm(BilinearForm&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_operator(std::move(other.m_operator)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Bilinear form to copy
       * @returns Reference to this bilinear form
       */
      BilinearForm& operator=(const BilinearForm& other)
      {
        if (this != &other)
        {
          m_u = other.m_u;
          m_v = other.m_v;
          m_operator = other.m_operator;
          m_assembly = other.m_assembly;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Bilinear form to move from
       * @returns Reference to this bilinear form
       */
      BilinearForm& operator=(BilinearForm&& other) noexcept
      {
        if (this != &other)
        {
          m_u = std::move(other.m_u);
          m_v = std::move(other.m_v);
          m_operator = std::move(other.m_operator);
          m_assembly = std::move(other.m_assembly);
        }
        return *this;
      }

      /**
       * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
       * v @f$.
       *
       * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
       * compute the action of the bilinear mapping @f$ a(u, v) @f$.
       *
       * @returns The action @f$ a(u, v) @f$ which the bilinear form takes
       * at @f$ ( u, v ) @f$.
       */
      template <class UData, class VData>
      constexpr
      ScalarType operator()(
        const GridFunction<TrialFES, UData>& u, const GridFunction<TestFES, VData>& v) const
      {
        return (this->getOperator() * v.getData()).dot(u.getData());
      }

      /**
       * @brief Gets the assembled sparse operator.
       * @returns Mutable reference to the assembled matrix
       */
      OperatorType& getOperator() override
      {
        return m_operator;
      }

      /**
       * @brief Gets the assembled sparse operator.
       * @returns Const reference to the assembled matrix
       */
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      void assemble() override
      {
        const auto& trialFES = getTrialFunction().getFiniteElementSpace();
        const auto& testFES = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(m_operator, {
          trialFES, testFES, this->getLocalIntegrators(), this->getGlobalIntegrators() });
      }

      const TrialFunction<SolutionType, TrialFES>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<TestFES>& getTestFunction() const override
      {
        return m_v.get();
      }

      BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<Solution, TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>> m_v;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  /// @brief Deduces the default sparse bilinear form type.
  template <class Solution, class TrialFES, class TestFES>
  BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
    -> BilinearForm<
        Solution, TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>
          ::Type>>;

  /**
   * @ingroup BilinearFormSpecializations
   * @brief Specialization of BilinearForm for a dense matrix operator.
   */
  template <class Solution, class TrialFES, class TestFES, class Scalar>
  class BilinearForm<Solution, TrialFES, TestFES, Math::Matrix<Scalar>> final
    : public BilinearFormBase<Math::Matrix<Scalar>>
  {
    using TrialFESContextType = typename FormLanguage::Traits<TrialFES>::ContextType;

    using TestFESContextType = typename FormLanguage::Traits<TestFES>::ContextType;

    public:
      /// @brief Solution vector type.
      using SolutionType = Solution;

      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// @brief Type of operator associated to the bilinear form.
      using OperatorType =
        Math::Matrix<ScalarType>;

      /// @brief Parent class type.
      using Parent =
        BilinearFormBase<OperatorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /// @brief Default assembly backend selected from the trial/test contexts.
      using DefaultAssemblyType =
        typename Assembly::Default<TrialFESContextType, TestFESContextType>
          ::template Type<OperatorType, BilinearForm>;

      /// @brief Assembly backend used by this bilinear form.
      using AssemblyType = DefaultAssemblyType;

      /**
       * @brief Constructs a dense bilinear form from trial and test functions.
       * @param[in] u Trial function
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Bilinear form to copy
       */
      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_operator(other.m_operator),
          m_assembly(other.m_assembly)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Bilinear form to move from
       */
      constexpr
      BilinearForm(BilinearForm&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_operator(std::move(other.m_operator)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Bilinear form to copy
       * @returns Reference to this bilinear form
       */
      BilinearForm& operator=(const BilinearForm& other)
      {
        if (this != &other)
        {
          m_u = other.m_u;
          m_v = other.m_v;
          m_operator = other.m_operator;
          m_assembly = other.m_assembly;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Bilinear form to move from
       * @returns Reference to this bilinear form
       */
      BilinearForm& operator=(BilinearForm&& other) noexcept
      {
        if (this != &other)
        {
          m_u = std::move(other.m_u);
          m_v = std::move(other.m_v);
          m_operator = std::move(other.m_operator);
          m_assembly = std::move(other.m_assembly);
        }
        return *this;
      }

      /**
       * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
       * v @f$.
       *
       * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
       * compute the action of the bilinear mapping @f$ a(u, v) @f$.
       *
       * @returns The action @f$ a(u, v) @f$ which the bilinear form takes
       * at @f$ ( u, v ) @f$.
       */
      template <class UData, class VData>
      constexpr
      ScalarType operator()(
        const GridFunction<TrialFES, UData>& u, const GridFunction<TestFES, VData>& v) const
      {
        return (this->getOperator() * v.getData()).dot(u.getData());
      }

      /**
       * @brief Gets the assembled dense operator.
       * @returns Mutable reference to the assembled matrix
       */
      OperatorType& getOperator() override
      {
        return m_operator;
      }

      /**
       * @brief Gets the assembled dense operator.
       * @returns Const reference to the assembled matrix
       */
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      void assemble() override
      {
        const auto& trialFES = getTrialFunction().getFiniteElementSpace();
        const auto& testFES = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(this->getOperator(), {
          trialFES, testFES, this->getLocalIntegrators(), this->getGlobalIntegrators() });
      }

      const TrialFunction<SolutionType, TrialFES>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<TestFES>& getTestFunction() const override
      {
        return m_v.get();
      }

      BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<Solution, TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>> m_v;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };
}

#endif
