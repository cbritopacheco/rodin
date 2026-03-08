#ifndef RODIN_PETSC_VARIATIONAL_BILINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_BILINEARFORM_H

#include <petscmacros.h>
#include <petscsystypes.h>

#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/PETSc/Variational/TrialFunction.h"

#include "Rodin/Variational/BilinearForm.h"

namespace Rodin::Variational
{
  template <class Solution, class TrialFES, class TestFES>
  class BilinearForm<Solution, TrialFES, TestFES, ::Mat> final
    : public BilinearFormBase<::Mat>
  {
    using TrialFESMeshType =
      typename FormLanguage::Traits<TrialFES>::MeshType;

    using TestFESMeshType =
      typename FormLanguage::Traits<TestFES>::MeshType;

    using TrialFESMeshContextType =
      typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

    using TestFESMeshContextType =
      typename FormLanguage::Traits<TestFESMeshType>::ContextType;

    public:
      using ScalarType =
        PetscScalar;

      template <class FES>
      using GridFunctionType =
        PETSc::Variational::GridFunction<FES>;

      using SolutionType = Solution;

      /// Type of operator associated to the bilinear form.
      using OperatorType = ::Mat;

      using DefaultAssembly =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
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
      BilinearForm(
        const TrialFunction<SolutionType, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u),
          m_v(v),
          m_operator(PETSC_NULLPTR)
      {
        PetscErrorCode ierr;
        const MPI_Comm comm = getPETScComm(u, v);
        ierr = MatCreate(comm, &m_operator);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v),
          m_assembly(other.m_assembly),
          m_operator(PETSC_NULLPTR)
      {
        PetscErrorCode ierr;
        const MPI_Comm comm = getPETScComm(m_u.get(), m_v.get());
        ierr = MatCreate(comm, &m_operator);
        assert(ierr == PETSC_SUCCESS);

        if (other.m_operator)
        {
          ierr = MatDuplicate(other.m_operator, MAT_COPY_VALUES, &m_operator);
          assert(ierr == PETSC_SUCCESS);
        }
        (void) ierr;
      }

      BilinearForm(BilinearForm&& other) noexcept
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_operator(other.m_operator)
      {
        other.m_operator = PETSC_NULLPTR;
      }

      ~BilinearForm() override
      {
        this->destroy();
      }

      BilinearForm& operator=(const BilinearForm& other)
      {
        if (this != &other)
        {
          this->destroy();

          static_cast<Parent&>(*this) = static_cast<const Parent&>(other);
          m_u = other.m_u;
          m_v = other.m_v;
          m_assembly = other.m_assembly;

          PetscErrorCode ierr;
          const MPI_Comm comm = getPETScComm(m_u.get(), m_v.get());
          ierr = MatCreate(comm, &m_operator);
          assert(ierr == PETSC_SUCCESS);

          if (other.m_operator)
          {
            ierr = MatDuplicate(other.m_operator, MAT_COPY_VALUES, &m_operator);
            assert(ierr == PETSC_SUCCESS);
          }
          (void) ierr;
        }
        return *this;
      }

      BilinearForm& operator=(BilinearForm&& other) noexcept
      {
        if (this != &other)
        {
          destroy();

          static_cast<Parent&>(*this) = std::move(static_cast<Parent&>(other));
          m_u = std::move(other.m_u);
          m_v = std::move(other.m_v);
          m_assembly = std::move(other.m_assembly);
          m_operator = other.m_operator;
          other.m_operator = PETSC_NULLPTR;
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
          const GridFunctionType<TrialFES>& u,
          const GridFunctionType<TestFES>& v) const
      {
        ScalarType result = 0;
        PetscErrorCode ierr;
        ::Vec tmp = PETSC_NULLPTR;

        ierr = VecDuplicate(u.getData(), &tmp);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatMult(this->getOperator(), u.getData(), tmp);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecDot(v.getData(), tmp, &result);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecDestroy(&tmp);
        assert(ierr == PETSC_SUCCESS);

        (void) ierr;
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

      void destroy() noexcept
      {
        if (m_operator)
        {
          PetscErrorCode ierr = MatDestroy(&m_operator);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          m_operator = PETSC_NULLPTR;
        }
      }

    private:
      static MPI_Comm getPETScComm(
          const TrialFunction<SolutionType, TrialFES>& u, const TestFunction<TestFES>& v)
      {
        const auto& trialMesh = u.getFiniteElementSpace().getMesh();
        const auto& testMesh  = v.getFiniteElementSpace().getMesh();
        const auto& trialCtx = trialMesh.getContext();
        const auto& testCtx  = testMesh.getContext();
        (void) testCtx;
        if constexpr (
            std::is_same_v<TrialFESMeshContextType, Context::Local> &&
            std::is_same_v<TestFESMeshContextType, Context::Local>)
        {
          return PETSC_COMM_SELF;
        }
        else
        {
          return trialCtx.getCommunicator();
        }
      }

      std::reference_wrapper<const TrialFunction<SolutionType, TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>> m_v;
      DefaultAssembly m_assembly;
      OperatorType m_operator;
  };

  template <class Solution, class TrialFES, class TestFES>
  BilinearForm(
      const PETSc::Variational::TrialFunction<Solution, TrialFES>& u,
      const PETSc::Variational::TestFunction<TestFES>& v)
    -> BilinearForm<Solution, TrialFES, TestFES, ::Mat>;
}

namespace Rodin::PETSc::Variational
{
  template <class Solution, class TrialFES, class TestFES>
  using BilinearForm =
    Rodin::Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>;
}

#endif
