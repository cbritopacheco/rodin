#ifndef RODIN_PETSC_ASSEMBLY_MPI_H
#define RODIN_PETSC_ASSEMBLY_MPI_H

#include <petsc.h>
#include <type_traits>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/MPI/Geometry.h"
#include "Rodin/MPI/Assembly.h"

namespace Rodin::Assembly
{
  // MPI assembly for PETSc Vec (linear form) -- skipping ghosts via Shard + Boost.MPI
  template <class FES>
  class MPI<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<
        ::Vec,
        Variational::LinearForm<FES, ::Vec>
      >
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar"
      );
      using VectorType     = ::Vec;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      using InputType      = typename Parent::InputType;

      void execute(VectorType& res, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const auto& fes   = input.getFES();
        const auto& mesh  = fes.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();

        // Compute local and global sizes via Boost.MPI
        const PetscInt localSize = PetscInt(fes.getShard().getSize());
        const PetscInt globalSize = PetscInt(fes.getSize());

        ierr = VecSetSizes(res, localSize, globalSize);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSet(res, PetscScalar(0));
        assert(ierr == PETSC_SUCCESS);

        // Local contributions skipping ghosts
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
                ierr = VecSetValue(res, PetscInt(dofs[i]), v, ADD_VALUES);
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

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  // MPI assembly for PETSc Mat (bilinear form) -- skipping ghosts via Shard + Boost.MPI
  template <class TrialFES, class TestFES>
  class MPI<::Mat, Variational::BilinearForm<TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<
        ::Mat,
        Variational::BilinearForm<TrialFES, TestFES, ::Mat>
      >
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
      using OperatorType     = ::Mat;
      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;
      using Parent           = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType        = typename Parent::InputType;

      void execute(OperatorType& A, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const auto& testFES  = input.getTestFES();
        const auto& trialFES = input.getTrialFES();
        const auto& mesh     = testFES.getMesh();
        const auto& shard    = mesh.getShard();
        const auto& ctx      = mesh.getContext();
        const auto& comm     = ctx.getCommunicator();

        // Compute local/global sizes
        const PetscInt localRows = PetscInt(testFES.getShard().getSize());
        const PetscInt localCols = PetscInt(trialFES.getShard().getSize());
        const PetscInt globalRows = PetscInt(testFES.getSize());
        const PetscInt globalCols = PetscInt(trialFES.getSize());

        // Create/init Mat
        ierr = MatSetSizes(A, localRows, localCols, globalRows, globalCols);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetFromOptions(A);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatMPIAIJSetPreallocation(A, PETSC_DECIDE, nullptr, PETSC_DECIDE, nullptr);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // Local element contributions skipping ghosts
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MPIIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            auto d = it->getDimension();
            auto i = it->getIndex();
            if (shard.isGhost(d, i))
              continue;

            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = testFES.getShard().getDOFs(d, i);
              const auto& cols = trialFES.getShard().getDOFs(d, i);
              for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
                for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                {
                  PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
            }
          }
        }

        // Global coupling contributions skipping ghosts
        // for (auto& bfi : input.getGlobalBFIs())
        // {
        //   const auto& testAttrs  = bfi.getTestAttributes();
        //   const auto& trialAttrs = bfi.getTrialAttributes();
        //   MPIIteration testSeq(mesh, bfi.getTestRegion());
        //   SequentialIteration trialSeq(mesh, bfi.getTrialRegion());
        //   for (auto teIt = testSeq.getIterator(); teIt; ++teIt)
        //   {
        //     auto dt = teIt->getDimension();
        //     auto tg = teIt->getIndex();
        //     if (shard.isGhost(dt, tg))
        //       continue;

        //     const auto& tmap = shard.getPolytopeMap(dt).left;
        //     auto tLoc = tmap.find(tg);
        //     if (tLoc == tmap.end())
        //       continue;
        //     auto tLid = tLoc->second;

        //     if (testAttrs.empty() || testAttrs.count(teIt->getAttribute()))
        //     {
        //       for (auto trIt = trialSeq.getIterator(); trIt; ++trIt)
        //       {
        //         auto dr = trIt->getDimension();
        //         auto rg = trIt->getIndex();
        //         if (shard.isGhost(dr, rg))
        //           continue;

        //         const auto& rmap = shard.getPolytopeMap(dr).left;
        //         auto rLoc = rmap.find(rg);
        //         if (rLoc == rmap.end())
        //           continue;
        //         auto rLid = rLoc->second;

        //         if (trialAttrs.empty() || trialAttrs.count(trIt->getAttribute()))
        //         {
        //           bfi.setPolytope(*trIt, *teIt);
        //           const auto& rows = testFES .getDOFs(dt, tLid);
        //           const auto& cols = trialFES.getDOFs(dr, rLid);
        //           for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
        //             for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
        //             {
        //               const PetscScalar v = PetscScalar(bfi.integrate(j, i));
        //               ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
        //               assert(ierr == PETSC_SUCCESS);
        //             }
        //         }
        //       }
        //     }
        //   }
        // }

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

#endif // RODIN_PETSC_ASSEMBLY_MPI_H
