#ifndef RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H
#define RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H

#include <petsc.h>
#include <type_traits>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Assembly/Sequential.h"

namespace Rodin::Assembly
{
  // Sequential assembly for PETSc Vec (linear form)
  template <class FES>
  class Sequential<
    ::Vec&,
    Variational::LinearForm<FES, ::Vec&>
  > final
    : public AssemblyBase<
        ::Vec&,
        Variational::LinearForm<FES, ::Vec&>
      >
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar for PETSc Vec assembly"
      );
      using VectorType      = ::Vec&;
      using LinearFormType  = Variational::LinearForm<FES, VectorType>;
      using Parent          = AssemblyBase<VectorType, LinearFormType>;
      using InputType       = typename Parent::InputType;

      void execute(VectorType& res, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const PetscInt n = PetscInt(input.getFES().getSize());

        ierr = VecSetSizes(res, n, PETSC_DECIDE);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecSetFromOptions(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecZeroEntries(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          Internal::SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              lfi.setPolytope(*it);
              const auto& dofs = input.getFES().getDOFs(it.getDimension(), it->getIndex());
              for (PetscInt l = 0; l < PetscInt(dofs.size()); ++l)
              {
                const PetscScalar v = PetscScalar(lfi.integrate(l));
                ierr = VecSetValue(res, PetscInt(dofs[l]), v, ADD_VALUES);
                PetscCallAbort(PETSC_COMM_SELF, ierr);
              }
            }
          }
        }

        ierr = VecAssemblyBegin(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecAssemblyEnd(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };


  // Sequential assembly for PETSc Mat (bilinear form)
  template <class TrialFES, class TestFES>
  class Sequential<
    ::Mat&,
    Variational::BilinearForm<TrialFES, TestFES, ::Mat&>
  > final
    : public AssemblyBase<
        ::Mat&,
        Variational::BilinearForm<TrialFES, TestFES, ::Mat&>
      >
  {

    public:
      using DotType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      using OperatorType    = ::Mat&;
      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;
      using Parent           = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType        = typename Parent::InputType;

      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly"
      );

      void execute(OperatorType& A, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const PetscInt m = PetscInt(input.getTestFES().getSize());
        const PetscInt n = PetscInt(input.getTrialFES().getSize());

        ierr = MatSetSizes(A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatSetFromOptions(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatSetUp(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatZeroEntries(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        const auto& mesh = input.getTrialFES().getMesh();
        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
              for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
                for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                {
                  const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
                  PetscCallAbort(PETSC_COMM_SELF, ierr);
                }
            }
          }
        }
        // Global contributions
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          Internal::SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          Internal::SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.empty() || testAttrs.count(teIt->getAttribute()))
            {
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.empty() || trialAttrs.count(trIt->getAttribute()))
                {
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
                    for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                    {
                      const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                      ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
                      PetscCallAbort(PETSC_COMM_SELF, ierr);
                    }
                }
              }
            }
          }
        }

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(PETSC_COMM_SELF, ierr);
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif // RODIN_ASSEMBLY_SEQUENTIAL_PETSC_H

