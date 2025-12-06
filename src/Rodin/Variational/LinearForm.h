/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearForm.h
 * @brief Linear form classes for finite element assembly.
 *
 * This file defines the LinearForm classes which represent linear functionals
 * @f$ l(v) @f$ in variational formulations. Linear forms are assembled into
 * vectors and form the right-hand side of finite element systems.
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_H
#define RODIN_VARIATIONAL_LINEARFORM_H

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Assembly/ForwardDecls.h"

#include "Rodin/Math/Vector.h"

#include "Exceptions/TestFunctionMismatchException.h"

#include "ForwardDecls.h"
#include "TestFunction.h"
#include "LinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class Vector>
  struct Traits<Variational::LinearFormBase<Vector>>
  {
    using VectorType = Vector;
  };

  template <class FES, class Vector>
  struct Traits<Variational::LinearForm<FES, Vector>>
  {
    using FESType = FES;
    using VectorType = Vector;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Base class for linear form representations.
   *
   * LinearFormBase provides the foundation for representing linear forms
   * @f$ l(v) : V \to \mathbb{R} @f$ in finite element computations. A linear
   * form is a linear functional that maps functions from the function space
   * @f$ V @f$ to real numbers, typically representing loads, sources, or
   * boundary data in variational formulations.
   *
   * @tparam Vector Vector type for the discrete representation
   *
   * ## Mathematical Foundation
   * A linear form @f$ l(v) @f$ satisfies:
   * - **Linearity**: @f$ l(\alpha v_1 + \beta v_2) = \alpha l(v_1) + \beta l(v_2) @f$
   * - **Boundedness**: @f$ |l(v)| \leq C \|v\|_V @f$ for some constant @f$ C @f$
   *
   * ## Discrete Representation
   * The discrete vector representation satisfies @f$ b_i = l(\psi_i) @f$ where
   * @f$ \psi_i @f$ are test basis functions.
   *
   * ## Common Examples
   * - **Load vector**: @f$ l(v) = \int_\Omega f \cdot v \, dx @f$
   * - **Neumann boundary**: @f$ l(v) = \int_{\partial\Omega} g \cdot v \, ds @f$
   * - **Point sources**: @f$ l(v) = \sum_k f_k v(x_k) @f$
   */
  template <class Vector>
  class LinearFormBase : public FormLanguage::Base
  {
    public:
      /// @brief Vector type for the discrete representation
      using VectorType = Vector;

      /// @brief Scalar type for vector entries
      using ScalarType = typename FormLanguage::Traits<VectorType>::ScalarType;

      /// @brief Base type for linear form integrators
      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      /// @brief List type for managing integrators
      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      /// @brief Parent class type
      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs a linear form with a default constructed vector.
       *
       * Creates a linear form instance with an owned vector that will store
       * the discrete representation of the linear form after assembly.
       */
      constexpr
      LinearFormBase() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other Linear form to copy
       *
       * Creates a copy including all integrators defining the linear form.
       */
      constexpr
      LinearFormBase(const LinearFormBase& other)
        : Parent(other),
          m_lfis(other.m_lfis)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Linear form to move from
       */
      constexpr
      LinearFormBase(LinearFormBase&& other)
        : Parent(std::move(other)),
          m_lfis(std::move(other.m_lfis))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Linear form to copy
       * @return Reference to this linear form
       */
      LinearFormBase& operator=(const LinearFormBase& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_lfis = other.m_lfis;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Linear form to move from
       * @return Reference to this linear form
       */
      LinearFormBase& operator=(LinearFormBase&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_lfis = std::move(other.m_lfis);
        }
        return *this;
      }

      /**
       * @brief Gets the list of integrators.
       * @return Reference to the list of linear form integrators
       *
       * Linear form integrators define contributions to the right-hand side:
       * @f$ b_i = \sum_{\text{integrators}} \int \text{integrand}(\psi_i) \, dx @f$
       */
      constexpr
      LinearFormIntegratorBaseListType& getIntegrators()
      {
        return m_lfis;
      }

      /**
       * @brief Gets the list of integrators (const version).
       * @return Const reference to the list of linear form integrators
       */
      constexpr
      const LinearFormIntegratorBaseListType& getIntegrators() const
      {
        return m_lfis;
      }

      /**
       * @brief Adds an integrator to this linear form.
       * @param[in] lfi Linear form integrator to add
       * @return Reference to this linear form for method chaining
       *
       * Adds a contribution @f$ l_i(v) @f$ to the linear form:
       * @f$ l(v) = l_1(v) + l_2(v) + \ldots @f$
       */
      constexpr
      LinearFormBase& operator+=(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(lfi);
        return *this;
      }

      /**
       * @brief Adds a list of integrators to this linear form.
       * @param[in] lfis List of linear form integrators to add
       * @return Reference to this linear form for method chaining
       */
      constexpr
      LinearFormBase& operator+=(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(lfis);
        return *this;
      }

      /**
       * @brief Subtracts an integrator from this linear form.
       * @param[in] lfi Linear form integrator to subtract
       * @return Reference to this linear form for method chaining
       *
       * Subtracts a contribution: @f$ l(v) = l(v) - l_i(v) @f$
       */
      constexpr
      LinearFormBase& operator-=(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(UnaryMinus(lfi));
        return *this;
      }

      /**
       * @brief Subtracts a list of integrators from this linear form.
       * @param[in] lfis List of linear form integrators to subtract
       * @return Reference to this linear form for method chaining
       */
      constexpr
      LinearFormBase& operator-=(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(UnaryMinus(lfis));
        return *this;
      }

      /**
       * @brief Assigns an integrator to this linear form.
       * @param[in] lfi Linear form integrator to assign
       * @return Reference to this linear form
       *
       * Replaces all existing integrators with a single new integrator.
       */
      constexpr
      LinearFormBase& operator=(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.clear();
        m_lfis.add(lfi);
        return *this;
      }

      constexpr
      LinearFormBase& operator=(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.clear();
        m_lfis.add(lfis);
        return *this;
      }

      /**
       * @brief Assembles the linear form.
       *
       * This method will assemble the underlying vector associated
       * the linear form.
       *
       * @see getVector()
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      virtual VectorType& getVector() = 0;

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      virtual const VectorType& getVector() const = 0;

      /**
       * @brief Gets the test function argument associated to this linear
       * form.
       */
      virtual const FormLanguage::Base& getTestFunction() const = 0;

      virtual LinearFormBase* copy() const noexcept override = 0;

    private:
      LinearFormIntegratorBaseListType m_lfis;
  };

  /**
   * @brief Represents a linear form defined over some finite element space
   *
   * An object of type LinearForm represents a linear map
   * @f[
   * \begin{aligned}
   *   L : V &\rightarrow \mathbb{R}\\
   *      v &\mapsto L(v)
   * \end{aligned}
   * @f]
   * where @f$ V @f$ is a finite element space.
   *
   * A linear form can be specified by from one or more
   * LinearFormIntegratorBase instances.
   */
  template <class FES>
  class LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>> final
    : public LinearFormBase<Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      using FESType =
        FES;

      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      using VectorType =
        Math::Vector<ScalarType>;

      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      using DefaultAssemblyType =
        typename Assembly::Default<FESMeshContextType>::template Type<VectorType, LinearForm>;

      using AssemblyType =
        DefaultAssemblyType;

      using Parent =
        LinearFormBase<VectorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      LinearForm(const TestFunction<FES>& v)
        : m_v(v)
      {}

      constexpr
      LinearForm(const LinearForm& other)
        : Parent(other),
          m_v(other.m_v),
          m_vector(other.m_vector),
          m_assembly(other.m_assembly)
      {}

      constexpr
      LinearForm(LinearForm&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_vector(std::move(other.m_vector)),
          m_assembly(std::move(other.m_assembly))
      {}

      LinearForm& operator=(const LinearForm& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_v = other.m_v;
          m_vector = other.m_vector;
          m_assembly = other.m_assembly;
        }
        return *this;
      }

      LinearForm& operator=(LinearForm&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_v = std::move(other.m_v);
          m_vector = std::move(other.m_vector);
          m_assembly = std::move(other.m_assembly);
        }
        return *this;
      }

      /**
       * @brief Evaluates the linear form at the function @f$ u @f$.
       *
       * Given a grid function @f$ u @f$, this function will compute the
       * action of the linear mapping @f$ L(u) @f$.
       *
       * @returns The value which the linear form takes at @f$ u @f$.
       */
      template <class Data>
      constexpr
      ScalarType operator()(const GridFunction<FES, Data>& u) const
      {
        return this->getVector().dot(u.getData());
      }

      void assemble() override
      {
        const auto& fes = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(this->getVector(), { fes, this->getIntegrators() });
      }

      VectorType& getVector() override
      {
        return m_vector;
      }

      const VectorType& getVector() const override
      {
        return m_vector;
      }

      const TestFunction<FES>& getTestFunction() const override
      {
        return m_v.get();
      }

      LinearForm* copy() const noexcept override
      {
        return new LinearForm(*this);
      }

    private:
      std::reference_wrapper<const TestFunction<FES>> m_v;
      VectorType m_vector;
      AssemblyType m_assembly;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LinearForm.
   * @param[in] v Reference to a TestFunction
   *
   * The constructor taking a single TestFunction reference deduces a
   * LinearForm with a default-constructed Math::Vector owned by the LinearForm
   * instance.
   */
  template <class FES>
  LinearForm(const TestFunction<FES>& v)
    -> LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>;
}

#endif
