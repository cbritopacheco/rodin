#ifndef RODIN_PETSC_VARIATIONAL_BILINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_BILINEARFORM_H

/**
 * @file
 * @brief PETSc specialization of variational bilinear forms.
 */

#include <petscmacros.h>
#include <petscsystypes.h>

#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/PETSc/Variational/TrialFunction.h"

#include "Rodin/Variational/BilinearForm.h"

namespace Rodin::Variational
{
  /**
   * @brief Bilinear form specialization assembled into a PETSc matrix.
   */
  template <class Solution, class TrialFES, class TestFES>
  class BilinearForm<Solution, TrialFES, TestFES, ::Mat> final
    : public BilinearFormBase<::Mat>
  {
    /// @brief Mesh type for the trial finite element space.
    using TrialFESMeshType =
      typename FormLanguage::Traits<TrialFES>::MeshType;

    /// @brief Mesh type for the test finite element space.
    using TestFESMeshType =
      typename FormLanguage::Traits<TestFES>::MeshType;

    /// @brief Context type (Local or MPI) for the trial mesh.
    using TrialFESMeshContextType =
      typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

    /// @brief Context type (Local or MPI) for the test mesh.
    using TestFESMeshContextType =
      typename FormLanguage::Traits<TestFESMeshType>::ContextType;

    public:
      /// @brief Scalar type (PETSc scalar).
      using ScalarType =
        PetscScalar;

      /// @brief Grid function type template for a given FES.
      template <class FES>
      using GridFunctionType =
        PETSc::Variational::GridFunction<FES>;

      /// @brief Solution type for the trial function.
      using SolutionType = Solution;

      /// @brief Type of operator associated to the bilinear form.
      using OperatorType = ::Mat;

      /// @brief Default assembly type for this bilinear form.
      using DefaultAssembly =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<OperatorType, BilinearForm>;

      /// @brief Parent class type.
      using Parent = BilinearFormBase<OperatorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a PETSc bilinear form for trial and test functions.
       * @param[in] u Trial function.
       * @param[in] v Test function.
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

      /// @brief Copy constructor (deep-copies the PETSc matrix).
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

      /// @brief Move constructor.
      BilinearForm(BilinearForm&& other) noexcept
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_operator(other.m_operator)
      {
        other.m_operator = PETSC_NULLPTR;
      }

      /// @brief Destructor; destroys the owned PETSc matrix.
      ~BilinearForm() override
      {
        this->destroy();
      }

      /**
       * @brief Copy assignment operator.
       * @param[in] other Bilinear form to copy.
       * @return Reference to this bilinear form.
       */
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

      /**
       * @brief Move assignment operator.
       * @param[in] other Bilinear form to move from.
       * @return Reference to this bilinear form.
       */
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

      /// @brief Returns a mutable reference to the assembled PETSc matrix.
      OperatorType& getOperator() override
      {
        return m_operator;
      }

      /// @brief Returns a read-only reference to the assembled PETSc matrix.
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      /// @brief Assembles the bilinear form into the PETSc matrix.
      void assemble() override
      {
        const auto& trialFES = getTrialFunction().getFiniteElementSpace();
        const auto& testFES = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(m_operator, {
          trialFES, testFES, this->getLocalIntegrators(), this->getGlobalIntegrators() });
      }

      /// @brief Returns a reference to the associated trial function.
      const TrialFunction<SolutionType, TrialFES>& getTrialFunction() const override
      {
        return m_u.get();
      }

      /// @brief Returns a reference to the associated test function.
      const TestFunction<TestFES>& getTestFunction() const override
      {
        return m_v.get();
      }

      /// @brief Creates a heap-allocated copy of this bilinear form.
      BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

      /// @brief Destroys the owned PETSc matrix, releasing resources.
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
      /**
       * @brief Deduces the MPI communicator from the trial and test meshes.
       * @param[in] u Trial function reference.
       * @param[in] v Test function reference.
       * @returns MPI communicator for PETSc object creation.
       */
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

      std::reference_wrapper<const TrialFunction<SolutionType, TrialFES>> m_u; ///< Reference to the trial function.
      std::reference_wrapper<const TestFunction<TestFES>> m_v;                 ///< Reference to the test function.
      DefaultAssembly m_assembly;                                              ///< Assembly strategy.
      OperatorType m_operator;                                                 ///< Owned PETSc matrix.
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for PETSc-backed BilinearForm.
   */
  template <class Solution, class TrialFES, class TestFES>
  BilinearForm(
      const PETSc::Variational::TrialFunction<Solution, TrialFES>& u,
      const PETSc::Variational::TestFunction<TestFES>& v)
    -> BilinearForm<Solution, TrialFES, TestFES, ::Mat>;
}

namespace Rodin::PETSc::Variational
{
  /**
   * @brief Convenient PETSc alias for Rodin::Variational::BilinearForm.
   */
  template <class Solution, class TrialFES, class TestFES>
  using BilinearForm =
    Rodin::Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>;
}

#endif
