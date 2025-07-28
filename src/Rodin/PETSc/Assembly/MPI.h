#ifndef RODIN_PETSC_ASSEMBLY_MPI_H
#define RODIN_PETSC_ASSEMBLY_MPI_H

#include <petsc.h>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/MPI/Assembly.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Math/Matrix.h"

namespace Rodin::Assembly
{
  template <class FES>
  class MPI<PETSc::Vector, Variational::LinearForm<FES, PETSc::Vector>> final
    : public AssemblyBase<
        PETSc::Vector, Variational::LinearForm<FES, PETSc::Vector>>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar"
      );
      using VectorType     = PETSc::Vector;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      using InputType      = typename Parent::InputType;

      void execute(VectorType& res, const InputType& input) const override
      {
        const auto& fes   = input.getFES();
        const auto& mesh  = fes.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const PetscInt localSize = PetscInt(fes.getShard().getSize());
        const PetscInt globalSize = PetscInt(fes.getSize());
        res.setSizes(localSize, globalSize)
           .setFromOptions()
           .zeroEntries();
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
              for (PetscInt i = 0; i < PetscInt(dofs.size()); ++i)
              {
                const PetscScalar v = PetscScalar(lfi.integrate(i));
                res.setValue(
                    PetscInt(fes.getGlobalIndex(dofs[i])), v, ADD_VALUES);
              }
            }
          }
        }
        res.assemblyBegin();
        res.assemblyEnd();
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  template <class Solution, class TrialFES, class TestFES>
  class MPI<PETSc::Matrix, Variational::BilinearForm<Solution, TrialFES, TestFES, PETSc::Matrix>> final
    : public AssemblyBase<PETSc::Matrix, Variational::BilinearForm<Solution, TrialFES, TestFES, PETSc::Matrix>>
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
      using OperatorType     = PETSc::Matrix;
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
        const PetscInt localRows = PetscInt(testFES.getShard().getSize());
        const PetscInt localCols = PetscInt(trialFES.getShard().getSize());
        const PetscInt globalRows = PetscInt(testFES.getSize());
        const PetscInt globalCols = PetscInt(trialFES.getSize());
        res.setSizes(localRows, localCols, globalRows, globalCols)
           .MPIAIJ.setPreallocation(
               PETSC_DECIDE, nullptr, PETSC_DECIDE, nullptr)
           .setFromOptions()
           .setUp()
           .zeroEntries();
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
              for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
              {
                for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                {
                  PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  res.setValue(
                      testFES.getGlobalIndex(rows[i]), trialFES.getGlobalIndex(cols[j]),
                      v, ADD_VALUES);
                }
              }
            }
          }
        }
        res.assemblyBegin(MAT_FINAL_ASSEMBLY);
        res.assemblyEnd(MAT_FINAL_ASSEMBLY);
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using MPI = Rodin::Assembly::MPI<LinearAlgebraType, Operand>;
}

#endif // RODIN_PETSC_ASSEMBLY_MPI_H
