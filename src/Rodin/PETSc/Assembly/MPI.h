#ifndef RODIN_PETSC_ASSEMBLY_MPI_H
#define RODIN_PETSC_ASSEMBLY_MPI_H

#include <petsc.h>
#include <petscmacros.h>
#include <petscmat.h>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/MPI/Assembly.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Math/Matrix.h"

namespace Rodin::Assembly
{
  template <class FES>
  class MPI<PETSc::Math::Vector, Variational::LinearForm<FES, PETSc::Math::Vector>> final
    : public AssemblyBase<
        PETSc::Math::Vector, Variational::LinearForm<FES, PETSc::Math::Vector>>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar"
      );
      using VectorType     = PETSc::Math::Vector;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      using InputType      = typename Parent::InputType;

      void execute(VectorType& res, const InputType& input) const override
      {
        const auto& fes = input.getFES();
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        const size_t localSize = fes.getShard().getSize();
        const size_t globalSize = fes.getSize();
        PetscErrorCode ierr;
        ierr = VecCreate(comm, &res);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetSizes(res, localSize, globalSize);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(res);
        assert(ierr == PETSC_SUCCESS);
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          MPIIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index i = it->getIndex();
            if (shard.isGhost(d, i))
              continue;
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              lfi.setPolytope(*it);
              const auto& dofs = fes.getShard().getDOFs(d, i);
              for (PetscInt i = 0; i < dofs.size(); ++i)
              {
                const PetscScalar v = PetscScalar(lfi.integrate(i));
                VecSetValue(res, fes.getGlobalIndex(dofs[i]), v, ADD_VALUES);
              }
            }
          }
        }
        ierr = VecAssemblyBegin(res);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(res);
        assert(ierr == PETSC_SUCCESS);
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  template <class Solution, class TrialFES, class TestFES>
  class MPI<PETSc::Math::Matrix, Variational::BilinearForm<Solution, TrialFES, TestFES, PETSc::Math::Matrix>> final
    : public AssemblyBase<PETSc::Math::Matrix, Variational::BilinearForm<Solution, TrialFES, TestFES, PETSc::Math::Matrix>>
  {
    public:
      using DotType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "DotType must be PetscScalar"
      );
      using OperatorType     = PETSc::Math::Matrix;
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;
      using Parent           = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType        = typename Parent::InputType;

      void execute(OperatorType& res, const InputType& input) const override
      {
        const auto& testFES  = input.getTestFES();
        const auto& trialFES = input.getTrialFES();
        const auto& mesh     = testFES.getMesh();
        const auto& shard    = mesh.getShard();
        const auto& ctx      = mesh.getContext();
        const auto& comm     = ctx.getCommunicator();
        const size_t localRows = testFES.getShard().getSize();
        const size_t localCols = trialFES.getShard().getSize();
        const size_t globalRows = testFES.getSize();
        const size_t globalCols = trialFES.getSize();
        PetscErrorCode ierr;
        ierr = MatCreate(comm, &res);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetSizes(res, localRows, localCols, globalRows, globalCols);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatMPIAIJSetPreallocation(res, PETSC_DECIDE, PETSC_NULLPTR, PETSC_DECIDE, PETSC_NULLPTR);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetUp(res);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatZeroEntries(res);
        assert(ierr == PETSC_SUCCESS);
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MPIIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index i = it->getIndex();

            if (shard.isGhost(d, i))
              continue;
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = testFES.getShard().getDOFs(d, i);
              const auto& cols = trialFES.getShard().getDOFs(d, i);
              for (PetscInt i = 0; i < rows.size(); ++i)
              {
                for (PetscInt j = 0; j < cols.size(); ++j)
                {
                  PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(
                      res,
                      testFES.getGlobalIndex(rows[i]),
                      trialFES.getGlobalIndex(cols[j]), v, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
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

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

namespace Rodin::PETSc::Math::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using MPI = Rodin::Assembly::MPI<LinearAlgebraType, Operand>;
}

#endif // RODIN_PETSC_ASSEMBLY_MPI_H
