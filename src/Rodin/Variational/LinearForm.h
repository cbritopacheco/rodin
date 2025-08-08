/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
  template <class Vector>
  class LinearFormBase : public FormLanguage::Base
  {
    public:
      using VectorType = Vector;

      using ScalarType = typename FormLanguage::Traits<VectorType>::ScalarType;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseListType = FormLanguage::List<LinearFormIntegratorBaseType>;

      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs a linear form with a default constructed vector
       * which is owned by the LinearFormBase instance.
       */
      constexpr
      LinearFormBase() = default;

      constexpr
      LinearFormBase(const LinearFormBase& other)
        : Parent(other),
          m_lfis(other.m_lfis)
      {}

      constexpr
      LinearFormBase(LinearFormBase&& other)
        : Parent(std::move(other)),
          m_lfis(std::move(other.m_lfis))
      {}

      LinearFormBase& operator=(const LinearFormBase& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_lfis = other.m_lfis;
        }
        return *this;
      }

      LinearFormBase& operator=(LinearFormBase&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_lfis = std::move(other.m_lfis);
        }
        return *this;
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

      constexpr
      LinearFormBase& operator+=(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(lfi);
        return *this;
      }

      constexpr
      LinearFormBase& operator+=(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(lfis);
        return *this;
      }

      constexpr
      LinearFormBase& operator-=(const LinearFormIntegratorBaseType& lfi)
      {
        if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
        m_lfis.add(UnaryMinus(lfi));
        return *this;
      }

      constexpr
      LinearFormBase& operator-=(const LinearFormIntegratorBaseListType& lfis)
      {
        m_lfis.add(UnaryMinus(lfis));
        return *this;
      }

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
