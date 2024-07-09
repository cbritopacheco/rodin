/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_H
#define RODIN_VARIATIONAL_LINEARFORM_H

#include "Rodin/Configure.h"

#include "Rodin/FormLanguage/List.h"

#include "Rodin/Assembly/ForwardDecls.h"
#include "Rodin/Assembly/Multithreaded.h"

#include "Rodin/Alert/MemberFunctionException.h"
#include "Exceptions/TestFunctionMismatchException.h"

#include "ForwardDecls.h"
#include "TestFunction.h"
#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
  template <class Vector>
  class LinearFormBase : public FormLanguage::Base
  {
    public:
      using VectorType = Vector;

      LinearFormBase() = default;

      LinearFormBase(const LinearFormBase& other)
        : FormLanguage::Base(other)
      {}

      LinearFormBase(LinearFormBase&& other)
        : FormLanguage::Base(std::move(other))
      {}

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
  template <class FES, class Vector>
  class LinearForm final
    : public LinearFormBase<Vector>
  {
    public:
      using FESType = FES;

      using VectorType = Vector;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using SequentialAssembly = Assembly::Sequential<VectorType, LinearForm>;

      using MultithreadedAssembly = Assembly::Multithreaded<VectorType, LinearForm>;

      using Parent = LinearFormBase<VectorType>;

      /**
       * @brief Constructs a linear form defined on some finite element
       * space
       * @param[in] fes Reference to the finite element space
       */
      constexpr
      LinearForm(std::reference_wrapper<const TestFunction<FES>> v)
        : m_v(v)
      {
#ifdef RODIN_MULTITHREADED
        m_assembly.reset(new MultithreadedAssembly);
#else
        m_assembly.reset(new SequentialAssembly);
#endif
      }

      constexpr
      LinearForm(const LinearForm& other)
        : Parent(other),
          m_v(other.m_v),
          m_assembly(other.m_assembly->copy()),
          m_lfis(other.m_lfis)
      {}

      constexpr
      LinearForm(LinearForm&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_lfis(std::move(other.m_lfis))
      {}

      /**
       * @brief Evaluates the linear form at the function @f$ u @f$.
       *
       * Given a grid function @f$ u @f$, this function will compute the
       * action of the linear mapping @f$ L(u) @f$.
       *
       * @returns The value which the linear form takes at @f$ u @f$.
       */
      constexpr
      ScalarType operator()(const GridFunction<FES>& u) const
      {
        const auto& weights = u.getWeights();
        if (!weights.has_value())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GridFunction weights have not been calculated. "
            << "Call " << Alert::Identifier::Function("setWeights()")
            << " on the GridFunction object."
            << Alert::Raise;
        }
        assert(weights.has_value());
        return getVector().dot(weights.value());
      }

      constexpr
      LinearFormIntegratorBaseListType& getIntegrators()
      {
        return m_lfis;
      }

      constexpr
      const LinearFormIntegratorBaseListType& getIntegrators() const
      {
        return m_lfis;
      }

      LinearForm& setAssembly(const Assembly::AssemblyBase<VectorType, LinearForm>& assembly)
      {
        m_assembly.reset(assembly.copy());
        return *this;
      }

      const Assembly::AssemblyBase<VectorType, LinearForm>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

      void assemble() override
      {
        const auto& fes = getTestFunction().getFiniteElementSpace();
        m_vector = getAssembly().execute({ fes, getIntegrators() });
      }

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      inline
      VectorType& getVector() override
      {
        return m_vector;
      }

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      inline
      const VectorType& getVector() const override
      {
        return m_vector;
      }

      inline
      const TestFunction<FES>& getTestFunction() const override
      {
        return m_v.get();
      }
      /**
       * @brief Builds the linear form the given LinearFormIntegratorBase
       * instance
       * @param[in] lfi Integrator which will be used to build the linear form.
       * @returns Reference to this (for method chaining)
       */
      virtual LinearForm& from(const LinearFormIntegratorBaseType& lfi)
      {
        m_lfis.clear();
        add(lfi).assemble();
        return *this;
      }

      virtual LinearForm& from(const LinearFormIntegratorBaseListType& lfi)
      {
        m_lfis.clear();
        add(lfi).assemble();
        return *this;
      }

      /**
       * @brief Builds the linear form the given LinearFormIntegratorBase
       * instance
       * @param[in] lfi Integrator which will be used to build the linear form.
       * @returns Reference to this (for method chaining)
       */
      virtual LinearForm& add(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(lfi);
        return *this;
      }

      virtual LinearForm& add(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(lfis);
        return *this;
      }

      virtual LinearForm& operator=(const LinearFormIntegratorBaseType& lfi)
      {
        from(lfi).assemble();
        return *this;
      }

      virtual LinearForm& operator=(const LinearFormIntegratorBaseListType& lfis)
      {
        from(lfis).assemble();
        return *this;
      }

      LinearForm* copy() const noexcept override
      {
        return new LinearForm(*this);
      }

    private:
      std::reference_wrapper<const TestFunction<FES>> m_v;
      std::unique_ptr<Assembly::AssemblyBase<VectorType, LinearForm>> m_assembly;
      LinearFormIntegratorBaseListType m_lfis;
      VectorType m_vector;
  };

  template <class FES>
  LinearForm(TestFunction<FES>&)
    -> LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>;
}

#include "LinearForm.hpp"

#endif
