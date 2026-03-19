#ifndef RODIN_PETSC_VARIATIONAL_LINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_LINEARFORM_H

/**
 * @file
 * @brief PETSc specialization of variational linear forms.
 */

#include <petscsystypes.h>

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Variational/TestFunction.h"

#include "Rodin/Variational/LinearForm.h"

namespace Rodin::Variational
{
  /**
   * @brief Linear form specialization assembled into a PETSc vector.
   */
  template <class FES>
  class LinearForm<FES, ::Vec> final
    : public LinearFormBase<::Vec>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Mesh type for the finite element space.
      using FESMeshType = typename FormLanguage::Traits<FESType>::MeshType;

      /// @brief Scalar type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief PETSc vector type.
      using VectorType = ::Vec;

      /// @brief Context type (Local or MPI).
      using ContextType = typename FormLanguage::Traits<FESMeshType>::ContextType;

      /// @brief Default assembly type for this linear form.
      using DefaultAssembly =
        typename Assembly::Default<ContextType>::template Type<VectorType, LinearForm>;

      /// @brief Parent class type.
      using Parent = LinearFormBase<VectorType>;

      using Parent::operator=;

      using Parent::operator+=;

      using Parent::operator-=;

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
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
       * @brief Evaluates the linear form at the function @f$ u @f$.
       *
       * Given a grid function @f$ u @f$, this function will compute the
       * action of the linear mapping @f$ L(u) @f$.
       *
       * @returns The value which the linear form takes at @f$ u @f$.
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

      /// @brief Assembles the linear form into the PETSc vector.
      void assemble() override
      {
        const auto& fes = getTestFunction().getFiniteElementSpace();
        m_assembly.execute(this->getVector(), { fes, this->getIntegrators() });
      }

      /// @brief Returns a mutable reference to the assembled PETSc vector.
      VectorType& getVector() override
      {
        return m_vector;
      }

      /// @brief Returns a read-only reference to the assembled PETSc vector.
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
