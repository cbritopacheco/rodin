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

      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs a linear form with a default constructed vector
       * which is owned by the LinearFormBase instance.
       */
      LinearFormBase()
        : m_vector(VectorType())
      {}

      /**
       * @brief Constructs a linear form with reference to vector which is not
       * owned by the LinearFormBase instance.
       */
      LinearFormBase(VectorType& vec)
        : m_vector(std::ref(vec))
      {}

      /**
       * @brief Constructs a linear form with a vector which is owned by the
       * LinearFormBase instance.
       */
      LinearFormBase(VectorType&& vec)
        : m_vector(std::move(vec))
      {}

      LinearFormBase(const LinearFormBase& other)
        : FormLanguage::Base(other),
          m_vector(other.m_vector)
      {}

      LinearFormBase(LinearFormBase&& other)
        : FormLanguage::Base(std::move(other)),
          m_vector(std::move(other.m_vector))
      {}

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      VectorType& getVector()
      {
        auto& ref = std::visit([](auto& m) -> VectorType& { return m; }, m_vector);
        return ref;
      }

      /**
       * @brief Gets the reference to the (local) associated vector
       * to the LinearForm.
       */
      const VectorType& getVector() const
      {
        const auto& ref = std::visit([](const auto& m) -> const VectorType& { return m; }, m_vector);
        return ref;
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
       * @brief Gets the test function argument associated to this linear
       * form.
       */
      virtual const FormLanguage::Base& getTestFunction() const = 0;

      virtual LinearFormBase* copy() const noexcept override = 0;

    private:
      std::variant<std::reference_wrapper<VectorType>, VectorType> m_vector;
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
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      LinearForm(const TestFunction<FES>& v)
        : LinearForm(v, Vector())
      {}

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * an non-owned vector.
       * @param[in] v Reference to a TestFunction
       * @param[in] vec Reference to a vector
       */
      constexpr
      LinearForm(const TestFunction<FES>& v, Vector& vec)
        : Parent(vec),
          m_v(v)
      {
#ifdef RODIN_MULTITHREADED
        m_assembly.reset(new MultithreadedAssembly);
#else
        m_assembly.reset(new SequentialAssembly);
#endif
      }

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * an owned vector.
       * @param[in] v Reference to a TestFunction
       * @param[in] vec Vector which will be owned by the LinearForm
       */
      constexpr
      LinearForm(const TestFunction<FES>& v, Vector&& vec)
        : Parent(std::move(vec)),
          m_v(v)
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
        return this->getVector().dot(weights.value());
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
        getAssembly().execute(this->getVector(), { fes, getIntegrators() });
      }

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
      LinearForm& from(const LinearFormIntegratorBaseType& lfi)
      {
        m_lfis.clear();
        add(lfi);
        return *this;
      }

      LinearForm& from(const LinearFormIntegratorBaseListType& lfi)
      {
        m_lfis.clear();
        add(lfi);
        return *this;
      }

      /**
       * @brief Builds the linear form the given LinearFormIntegratorBase
       * instance
       * @param[in] lfi Integrator which will be used to build the linear form.
       * @returns Reference to this (for method chaining)
       */
      LinearForm& add(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(lfi);
        return *this;
      }

      LinearForm& add(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(lfis);
        return *this;
      }

      LinearForm& operator=(const LinearFormIntegratorBaseType& lfi)
      {
        from(lfi);
        return *this;
      }

      LinearForm& operator=(const LinearFormIntegratorBaseListType& lfis)
      {
        from(lfis);
        return *this;
      }

      LinearForm* copy() const noexcept override
      {
        return new LinearForm(*this);
      }

      LinearForm& clear()
      {
        m_lfis.clear();
        return *this;
      }

    private:
      std::reference_wrapper<const TestFunction<FES>> m_v;
      std::unique_ptr<Assembly::AssemblyBase<VectorType, LinearForm>> m_assembly;
      LinearFormIntegratorBaseListType m_lfis;
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

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LinearForm.
   * @param[in] v Reference to a TestFunction
   * @param[in] vec Vector reference, not owned by the LinearForm
   *
   * The constructor taking a TestFunction reference and a vector reference
   * deduces a LinearForm with a Vector reference which is not owned by
   * the LinearForm instance.
   */
  template <class FES, class Vector>
  LinearForm(const TestFunction<FES>& v, Vector& vec)
    -> LinearForm<FES, Vector>;

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LinearForm.
   * @param[in] v Reference to a TestFunction
   * @param[in] vec Vector which will be owned by the LinearForm
   *
   * The constructor taking a TestFunction reference and a vector r-value
   * reference deduces a LinearForm with a Vector value which is owned by
   * the LinearForm instance.
   */
  template <class FES, class Vector>
  LinearForm(const TestFunction<FES>& v, Vector&& vec)
    -> LinearForm<FES, Vector>;
}

#include "LinearForm.hpp"

#endif
