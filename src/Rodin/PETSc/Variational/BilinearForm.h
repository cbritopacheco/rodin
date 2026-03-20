#ifndef RODIN_PETSC_VARIATIONAL_BILINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_BILINEARFORM_H

/**
 * @file BilinearForm.h
 * @brief PETSc specialization of variational bilinear forms.
 *
 * Provides the partial specialization of
 * @ref Rodin::Variational::BilinearForm that assembles bilinear-form
 * contributions into a PETSc `Mat`.  The resulting matrix represents
 * the stiffness (or system) matrix @f$ A @f$ of a discrete finite
 * element system @f$ A\mathbf{x} = \mathbf{b} @f$.
 *
 * ## Mathematical Background
 *
 * A bilinear form @f$ a : U_h \times V_h \to \mathbb{R} @f$ acts on
 * a trial function @f$ u \in U_h @f$ and a test function @f$ v \in V_h @f$.
 * After discretisation, its entries are:
 * @f[
 *   A_{ij} = a(\psi_j, \phi_i), \quad
 *   i = 1,\ldots,M, \;\; j = 1,\ldots,N
 * @f]
 * where @f$ \psi_j @f$ and @f$ \phi_i @f$ are trial and test basis
 * functions, respectively.
 *
 * @see Rodin::PETSc::Variational::LinearForm,
 *      Rodin::PETSc::Variational::Problem,
 *      Rodin::PETSc::Variational::TrialFunction,
 *      Rodin::PETSc::Variational::TestFunction
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
   * @brief Bilinear form specialization that assembles into a PETSc matrix.
   *
   * Owns a PETSc `Mat` representing the system matrix @f$ A @f$.  The
   * matrix is created in the constructor (using an MPI communicator
   * deduced from the trial and test meshes) and destroyed in the
   * destructor.
   *
   * @tparam Solution  Solution grid function type for the trial function.
   * @tparam TrialFES  Finite element space type of the trial function.
   * @tparam TestFES   Finite element space type of the test function.
   *
   * @see Rodin::Variational::BilinearFormBase,
   *      Rodin::PETSc::Variational::BilinearForm
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

    /// @brief Context type (either @ref Rodin::Context::Local or
    ///        @ref Rodin::Context::MPI) for the trial mesh.
    using TrialFESMeshContextType =
      typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

    /// @brief Context type (either @ref Rodin::Context::Local or
    ///        @ref Rodin::Context::MPI) for the test mesh.
    using TestFESMeshContextType =
      typename FormLanguage::Traits<TestFESMeshType>::ContextType;

    public:
      /// @brief Scalar type (`PetscScalar`) for matrix entries.
      using ScalarType =
        PetscScalar;

      /// @brief Template alias mapping a finite element space to its
      ///        PETSc grid function type.
      template <class FES>
      using GridFunctionType =
        PETSc::Variational::GridFunction<FES>;

      /// @brief Solution (grid function) type associated with the trial
      ///        function.
      using SolutionType = Solution;

      /// @brief PETSc matrix type (`::Mat`) used to store the system
      ///        matrix @f$ A @f$.
      using OperatorType = ::Mat;

      /// @brief Default assembly strategy deduced from the trial and test
      ///        mesh context types.
      using DefaultAssembly =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<OperatorType, BilinearForm>;

      /// @brief Parent class providing the generic `BilinearFormBase<Mat>`
      ///        interface.
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
       * @brief Evaluates the bilinear form at grid functions @f$ u @f$
       *        and @f$ v @f$.
       *
       * Computes @f$ a(u, v) = \mathbf{v}^\top A \mathbf{u} @f$ by
       * performing a matrix–vector product `MatMult(A, u, tmp)` followed
       * by `VecDot(v, tmp, &result)`.
       *
       * @param[in] u Trial grid function @f$ u @f$.
       * @param[in] v Test grid function @f$ v @f$.
       * @returns The scalar value @f$ a(u, v) @f$.
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

      /// @brief Returns a mutable reference to the assembled PETSc system matrix @f$ A @f$.
      OperatorType& getOperator() override
      {
        return m_operator;
      }

      /// @brief Returns a read-only reference to the assembled PETSc system matrix @f$ A @f$.
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      /// @brief Assembles the bilinear form into the PETSc matrix using the
      ///        @ref DefaultAssembly strategy with both local and global integrators.
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
       *
       * When both meshes use `Context::Local`, returns `PETSC_COMM_SELF`.
       * Otherwise, returns the communicator from the trial mesh context.
       *
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
