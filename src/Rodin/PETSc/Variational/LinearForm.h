#ifndef RODIN_PETSC_VARIATIONAL_LINEARFORM_H
#define RODIN_PETSC_VARIATIONAL_LINEARFORM_H

#include <petscsystypes.h>

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Variational/TestFunction.h"

#include "Rodin/Variational/LinearForm.h"

namespace Rodin::Variational
{
  template <class FES>
  class LinearForm<FES, ::Vec> final
    : public LinearFormBase<::Vec>
  {
    public:
      using FESType = FES;

      using FESMeshType = typename FormLanguage::Traits<FESType>::MeshType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using VectorType = ::Vec;

      using ContextType = typename FormLanguage::Traits<FESMeshType>::ContextType;

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
      LinearForm(const TestFunction<FES>& v)
        : m_v(v),
          m_vector(PETSC_NULLPTR)
      {
        PetscErrorCode ierr = VecCreate(PETSC_COMM_SELF, &m_vector);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

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

      LinearForm(LinearForm&& other) noexcept
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_vector(other.m_vector)
      {
        other.m_vector = PETSC_NULLPTR;
      }

      ~LinearForm() override
      {
        destroy();
      }

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
      std::reference_wrapper<const TestFunction<FES>> m_v;
      DefaultAssembly m_assembly;
      VectorType m_vector;
  };

  template <class FES>
  LinearForm(const PETSc::Variational::TestFunction<FES>&) -> LinearForm<FES, ::Vec>;
}

namespace Rodin::PETSc::Variational
{
  template <class FES>
  using LinearForm = Rodin::Variational::LinearForm<FES, ::Vec>;
}

#endif

