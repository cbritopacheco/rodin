#ifndef RODIN_PETSC_VARIATIONAL_LINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_LINEARFORM_H

#include "Rodin/PETSc/Variational/GridFunction.h"
#include "Rodin/PETSc/Math/Vector.h"

#include "Rodin/Variational/LinearForm.h"

namespace Rodin::Variational
{
  template <class FES>
  class LinearForm<FES, PETSc::Vector> final
    : public LinearFormBase<PETSc::Vector>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using VectorType = PETSc::Vector;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using DefaultAssembly =
        typename Assembly::Default<ContextType>::template Type<VectorType, LinearForm>;

      using Parent = LinearFormBase<VectorType>;

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
        : m_v(v),
          m_assembly(new DefaultAssembly)
      {}

      constexpr
      LinearForm(const LinearForm& other)
        : Parent(other),
          m_v(other.m_v),
          m_assembly(other.m_assembly->copy()),
          m_vector(other.m_vector)
      {}

      constexpr
      LinearForm(LinearForm&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_vector(std::move(other.m_vector))
      {}

      LinearForm& operator=(const LinearForm& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_v = other.m_v;
          m_assembly.reset(other.m_assembly->copy());
          m_vector = other.m_vector;
        }
        return *this;
      }

      LinearForm& operator=(LinearForm&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_v = std::move(other.m_v);
          m_assembly = std::move(other.m_assembly);
          m_vector = std::move(other.m_vector);
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
      ScalarType operator()(const GridFunction<FES, PETSc::Vector>& u) const
      {
        ScalarType result;
        this->getVector().dot(u.getData(), &result);
        return result;
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
      DefaultAssembly m_assembly;
      VectorType m_vector;
  };
}

namespace Rodin::PETSc::Variational
{
  template <class FES>
  using LinearForm = Rodin::Variational::LinearForm<FES, PETSc::Vector>;
}

#endif

