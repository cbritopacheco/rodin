#ifndef RODIN_PETSC_VARIATIONAL_BILINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_BILINEARFORM_H

#include "Rodin/PETSc/Variational/GridFunction.h"
#include "Rodin/PETSc/Math/Matrix.h"

#include "Rodin/Variational/BilinearForm.h"
#include <cstdint>
#include <petscsystypes.h>

namespace Rodin::Variational
{
  template <class Solution, class TrialFES, class TestFES>
  class BilinearForm<Solution, TrialFES, TestFES, ::Mat> final
    : public BilinearFormBase<::Mat>
  {
    using TrialFESContextType = typename FormLanguage::Traits<TrialFES>::ContextType;

    using TestFESContextType = typename FormLanguage::Traits<TestFES>::ContextType;

    public:
      using SolutionType =
        Solution;

      using ScalarType =
        PetscScalar;

      template <class FES>
      using GridFunctionType =
        PETSc::Variational::GridFunction<FES>;

      /// Type of operator associated to the bilinear form.
      using OperatorType = ::Mat;

      using DefaultAssembly =
        typename Assembly::Default<TrialFESContextType, TestFESContextType>
          ::template Type<OperatorType, BilinearForm>;

      /// Parent class.
      using Parent = BilinearFormBase<OperatorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      BilinearForm(const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v)
      {}

      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_assembly(other.m_assembly),
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
          const GridFunctionType<TrialFES>& u, const GridFunctionType<TestFES>& v) const
      {
        ScalarType result;
        PetscErrorCode ierr;
        ::Vec tmp;
        ierr = VecDuplicate(tmp, u.getData());
        assert(ierr == PETSC_SUCCESS);
        ierr = MatMult(this->getOperator(), v.getData(), tmp);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecDot(tmp, u.getData(), &result);
        assert(ierr == PETSC_SUCCESS);
        return result;
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
      DefaultAssembly m_assembly;
      OperatorType m_operator;
  };
}

namespace Rodin::PETSc::Math::Variational
{
  template <class Solution, class TrialFES, class TestFES>
  using BilinearForm =
    Rodin::Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>;
}

#endif
