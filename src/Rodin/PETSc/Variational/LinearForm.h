#ifndef RODIN_PETSC_VARIATIONAL_LINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_LINEARFORM_H

/**
 * @file LinearForm.h
 * @brief PETSc specialization of variational linear forms.
 *
 * Provides the partial specialization of
 * @ref Rodin::Variational::LinearForm that assembles linear-form
 * contributions into a PETSc `Vec`.  The resulting vector represents
 * the right-hand side @f$ \mathbf{b} @f$ of a discrete finite element
 * system @f$ A\mathbf{x} = \mathbf{b} @f$.
 *
 * ## Mathematical Background
 *
 * A linear form @f$ L : V_h \to \mathbb{R} @f$ is a mapping that takes
 * a test function @f$ v \in V_h @f$ and produces a scalar.  After
 * discretisation the action of @f$ L @f$ on the basis functions
 * yields a load vector:
 * @f[
 *   b_i = L(\phi_i), \quad i = 1,\ldots,N
 * @f]
 *
 * This specialization stores @f$ \mathbf{b} @f$ in a PETSc `Vec` and
 * evaluates @f$ L(u_h) = \mathbf{b}^\top \mathbf{u} @f$ via `VecDot`.
 *
 * @see Rodin::PETSc::Variational::TestFunction,
 *      Rodin::PETSc::Variational::BilinearForm,
 *      Rodin::PETSc::Variational::Problem
 */

#include <petscsystypes.h>

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Variational/TestFunction.h"

#include "Rodin/Variational/LinearForm.h"

namespace Rodin::Variational
{
  /**
   * @brief Linear form specialization that assembles into a PETSc vector.
   *
   * Owns a PETSc `Vec` representing the load vector @f$ \mathbf{b} @f$
   * of the finite element system.  The vector is created in the
   * constructor and destroyed in the destructor.
   *
   * @tparam FES Finite element space type of the associated test function.
   *
   * @see Rodin::Variational::LinearFormBase,
   *      Rodin::PETSc::Variational::LinearForm
   */
  template <class FES>
  class LinearForm<FES, ::Vec> final
    : public LinearFormBase<::Vec>
  {
    public:
      /// @brief Finite element space type of the associated test function.
      using FESType = FES;

      /// @brief Mesh type underlying the finite element space.
      using FESMeshType = typename FormLanguage::Traits<FESType>::MeshType;

      /// @brief Scalar type of the DOF coefficients (`PetscScalar`).
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief PETSc vector type (`::Vec`) used to store the load vector @f$ \mathbf{b} @f$.
      using VectorType = ::Vec;

      /// @brief Context type (either @ref Rodin::Context::Local or
      ///        @ref Rodin::Context::MPI) determined by the mesh.
      using ContextType = typename FormLanguage::Traits<FESMeshType>::ContextType;

      /// @brief Default assembly strategy deduced from the context type.
      ///        Sequential assembly for local contexts, MPI assembly for
      ///        distributed contexts.
      using DefaultAssembly =
        typename Assembly::Default<ContextType>::template Type<VectorType, LinearForm>;

      /// @brief Parent class providing the generic `LinearFormBase<Vec>` interface.
      using Parent = LinearFormBase<VectorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a linear form associated with the given test
       *        function, initialising an empty PETSc vector.
       *
       * The vector is created on `PETSC_COMM_SELF`; its sizes will be
       * set during assembly.
       *
       * @param[in] v Reference to a PETSc test function.
       */
      LinearForm(const TestFunction<FES>& v)
        : m_v(v),
          m_vector(PETSC_NULLPTR)
      {
        PetscErrorCode ierr = VecCreate(PETSC_COMM_SELF, &m_vector);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

      /// @brief Copy constructor (deep-copies the PETSc vector).
      LinearForm(const LinearForm& other)
        : Parent(other),
          m_v(other.m_v),
          m_assembly(other.m_assembly),
          m_vector(PETSC_NULLPTR)
      {
        PetscErrorCode ierr;
        if (other.m_vector)
        {
          ierr = VecDuplicate(other.m_vector, &m_vector);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecCopy(other.m_vector, m_vector);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          ierr = VecCreate(PETSC_COMM_SELF, &m_vector);
          assert(ierr == PETSC_SUCCESS);
        }
        (void) ierr;
      }

      /// @brief Move constructor.
      LinearForm(LinearForm&& other) noexcept
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_vector(other.m_vector)
      {
        other.m_vector = PETSC_NULLPTR;
      }

      /// @brief Destructor; destroys the owned PETSc vector.
      ~LinearForm() override
      {
        destroy();
      }

      /**
       * @brief Copy assignment operator.
       * @param[in] other Linear form to copy.
       * @return Reference to this linear form.
       */
      LinearForm& operator=(const LinearForm& other)
      {
        if (this != &other)
        {
          destroy();

          Parent::operator=(other);
          m_v = other.m_v;
          m_assembly = other.m_assembly;

          PetscErrorCode ierr;
          if (other.m_vector)
          {
            ierr = VecDuplicate(other.m_vector, &m_vector);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecCopy(other.m_vector, m_vector);
            assert(ierr == PETSC_SUCCESS);
          }
          else
          {
            ierr = VecCreate(PETSC_COMM_SELF, &m_vector);
            assert(ierr == PETSC_SUCCESS);
          }
          (void) ierr;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Linear form to move from.
       * @return Reference to this linear form.
       */
      LinearForm& operator=(LinearForm&& other) noexcept
      {
        if (this != &other)
        {
          destroy();

          Parent::operator=(std::move(other));
          m_v = std::move(other.m_v);
          m_assembly = std::move(other.m_assembly);
          m_vector = other.m_vector;
          other.m_vector = PETSC_NULLPTR;
        }
        return *this;
      }

      /**
       * @brief Evaluates the linear form at a grid function @f$ u_h @f$.
       *
       * Computes the action @f$ L(u_h) = \mathbf{b}^\top \mathbf{u} @f$
       * via `VecDot(b, u, &result)`.
       *
       * @param[in] u The grid function @f$ u_h @f$ to evaluate at.
       * @returns The scalar value @f$ L(u_h) @f$.
       */
      ScalarType operator()(const GridFunction<FES, ::Vec>& u) const
      {
        ScalarType result;
        PetscErrorCode ierr;
        ierr = VecDot(this->getVector(), u.getData(), &result);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return result;
      }

      /// @brief Assembles the linear form into the owned PETSc vector
      ///        using the @ref DefaultAssembly strategy.
      void assemble() override
      {
        const auto& fes = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(this->getVector(), { fes, this->getIntegrators() });
      }

      /// @brief Returns a mutable reference to the assembled PETSc load vector @f$ \mathbf{b} @f$.
      VectorType& getVector() override
      {
        return m_vector;
      }

      /// @brief Returns a read-only reference to the assembled PETSc load vector @f$ \mathbf{b} @f$.
      const VectorType& getVector() const override
      {
        return m_vector;
      }

      /// @brief Returns a reference to the associated test function.
      const TestFunction<FES>& getTestFunction() const override
      {
        return m_v.get();
      }

      /// @brief Creates a heap-allocated copy of this linear form.
      LinearForm* copy() const noexcept override
      {
        return new LinearForm(*this);
      }

      /// @brief Destroys the owned PETSc vector, releasing resources.
      void destroy() noexcept
      {
        if (m_vector)
        {
          PetscErrorCode ierr = VecDestroy(&m_vector);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          m_vector = PETSC_NULLPTR;
        }
      }

    private:
      std::reference_wrapper<const TestFunction<FES>> m_v; ///< Reference to the test function.
      DefaultAssembly m_assembly;                          ///< Assembly strategy.
      VectorType m_vector;                                 ///< Owned PETSc vector.
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for PETSc-backed LinearForm.
   */
  template <class FES>
  LinearForm(const PETSc::Variational::TestFunction<FES>&) -> LinearForm<FES, ::Vec>;
}

namespace Rodin::PETSc::Variational
{
  /**
   * @brief Convenient PETSc alias for Rodin::Variational::LinearForm.
   */
  template <class FES>
  using LinearForm = Rodin::Variational::LinearForm<FES, ::Vec>;
}

#endif
