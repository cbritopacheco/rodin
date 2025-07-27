#ifndef RODIN_PETSC_VARIATIONAL_BILINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_BILINEARFORM_H

#include "Rodin/PETSc/Variational/GridFunction.h"
#include "Rodin/PETSc/Math/Matrix.h"

#include "Rodin/Variational/BilinearForm.h"

namespace Rodin::Variational
{
  template <class Solution, class TrialFES, class TestFES>
  class BilinearForm<Solution, TrialFES, TestFES, PETSc::Matrix> final
    : public BilinearFormBase<PETSc::Matrix>
  {
    using TrialFESContextType = typename FormLanguage::Traits<TrialFES>::ContextType;

    using TestFESContextType = typename FormLanguage::Traits<TestFES>::ContextType;

    public:
      using SolutionType = Solution;

      using ScalarType = PetscScalar;

      /// Type of operator associated to the bilinear form.
      using OperatorType = PETSc::Matrix;

      /// Parent class.
      using Parent = BilinearFormBase<OperatorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      using DefaultAssembly =
        typename Assembly::Default<TrialFESContextType, TestFESContextType>
          ::template Type<OperatorType, BilinearForm>;

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v),
          m_assembly(new DefaultAssembly)
      {}

      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_assembly(other.m_assembly->copy()),
          m_operator(other.m_operator)
      {}

      constexpr
      BilinearForm(BilinearForm&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_operator(other.m_operator)
      {}

      BilinearForm& operator=(const BilinearForm& bf)
      {
        if (this != &bf)
        {
          m_u = bf.m_u;
          m_v = bf.m_v;
          m_assembly.reset(bf.m_assembly->copy());
          m_operator = bf.m_operator;
        }
        return *this;
      }

      BilinearForm& operator=(BilinearForm&& bf) noexcept
      {
        if (this != &bf)
        {
          m_u = std::move(bf.m_u);
          m_v = std::move(bf.m_v);
          m_assembly = std::move(bf.m_assembly);
          m_operator = std::move(bf.m_operator);
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
      ScalarType operator()(
        const PETSc::GridFunction<TrialFES>& u, const PETSc::GridFunction<TestFES>& v) const
      {
        PETSc::Vector tmp;
        u.duplicate(tmp);
        this->getOperator().mult(v.getData(), tmp);
        ScalarType result;
        tmp.dot(u.getData(), &result);
        return result;
      }

      const Assembly::AssemblyBase<OperatorType, BilinearForm>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

      BilinearForm& setAssembly(const Assembly::AssemblyBase<OperatorType, BilinearForm>& assembly)
      {
        m_assembly.reset(assembly.copy());
        return *this;
      }

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
        const auto& trialFES = getTrialFunction().getFiniteElementSpace();
        const auto& testFES = getTestFunction().getFiniteElementSpace();
        this->getAssembly().execute(m_operator, {
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
      std::reference_wrapper<const TrialFunction<Solution, TrialFES>>   m_u;
      std::reference_wrapper<const TestFunction<TestFES>>               m_v;
      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearForm>> m_assembly;
      OperatorType m_operator;
  };
}

#endif
