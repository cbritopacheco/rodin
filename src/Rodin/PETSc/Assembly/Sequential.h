#ifndef RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H
#define RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/Sequential.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Assembly
{
  // Sequential assembly for PETSc Vec (linear form)
  template <class FES>
  class Sequential<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      using VectorType = ::Vec;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent = AssemblyBase<VectorType, LinearFormType>;
      using InputType = typename Parent::InputType;

      static_assert(std::is_same_v<ScalarType, PetscScalar>);

      void execute(VectorType& res, const InputType& input) const override
      {
        assert(res);
        const size_t n = input.getFES().getSize();

        PetscErrorCode ierr;

        ierr = VecSetSizes(res, n, PETSC_DECIDE);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(res);
        assert(ierr == PETSC_SUCCESS);

        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              lfi.setPolytope(*it);
              const auto& dofs = input.getFES().getDOFs(it.getDimension(), it->getIndex());
              for (PetscInt l = 0; l < dofs.size(); ++l)
              {
                const PetscScalar v = PetscScalar(lfi.integrate(l));
                ierr = VecSetValue(res, dofs[l], v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          }
        }

        ierr = VecAssemblyBegin(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyEnd(res);
        assert(ierr == PETSC_SUCCESS);
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };


  // Sequential assembly for PETSc Mat (bilinear form)
  template <class Solution, class TrialFES, class TestFES>
  class Sequential<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      using DotType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      using OperatorType = ::Mat;
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType = typename Parent::InputType;

      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly"
      );

      void execute(OperatorType& res, const InputType& input) const override
      {
        assert(res);
        const size_t m = input.getTestFES().getSize();
        const size_t n = input.getTrialFES().getSize();

        PetscErrorCode ierr;
        ierr = MatSetSizes(res, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSeqAIJSetPreallocation(res, PETSC_DETERMINE, PETSC_NULLPTR);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(res);
        assert(ierr == PETSC_SUCCESS);

        const auto& mesh = input.getTrialFES().getMesh();
        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index i = it->getIndex();

            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(d, i);
              const auto& cols = input.getTrialFES().getDOFs(d, i);
              for (PetscInt i = 0; i < rows.size(); ++i)
              {
                for (PetscInt j = 0; j < cols.size(); ++j)
                {
                  const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(res, rows[i], cols[j], v, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
              }
            }
          }
        }

        // Global contributions
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());
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
                  for (PetscInt i = 0; i < rows.size(); ++i)
                    for (PetscInt j = 0; j < cols.size(); ++j)
                    {
                      const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                      ierr = MatSetValue(res, rows[i], cols[j], v, ADD_VALUES);
                      assert(ierr == PETSC_SUCCESS);
                    }
                }
              }
            }
          }
        }

        ierr = MatAssemblyBegin(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatAssemblyEnd(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

namespace Rodin::PETSc
{
  template <class LinearAlgebraType, class Operand>
  using Sequential = Assembly::Sequential<LinearAlgebraType, Operand>;
}

#endif // RODIN_ASSEMBLY_SEQUENTIAL_PETSC_H

