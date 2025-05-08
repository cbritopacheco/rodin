#ifndef RODIN_PETSC_ASSEMBLY_MPI_H
#define RODIN_PETSC_ASSEMBLY_MPI_H

#include <petsc.h>
#include <type_traits>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/MPI.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Geometry/MPI/Mesh.h"

namespace Rodin::Assembly
{
  // MPI assembly for PETSc Vec (linear form) -- skipping ghosts via Shard + Boost.MPI
  template <class FES>
  class MPI<
    ::Vec&,
    Variational::LinearForm<FES, ::Vec&>
  > final
    : public AssemblyBase<
        ::Vec&,
        Variational::LinearForm<FES, ::Vec&>
      >
  {
    public:
      using ScalarType     = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar"
      );
      using VectorType     = ::Vec&;
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
        const PetscInt localSize = PetscInt(fes.getSize());
        const PetscInt globalSize = boost::mpi::all_reduce(comm, localSize, std::plus<PetscInt>());

        ierr = VecSetSizes(res, localSize, globalSize);
        PetscCallAbort(comm, ierr);
        ierr = VecSetFromOptions(res);
        PetscCallAbort(comm, ierr);
        ierr = VecSet(res, PetscScalar(0));
        PetscCallAbort(comm, ierr);

        // Local contributions skipping ghosts
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          Internal::SequentialIteration seq(mesh, lfi.getRegion());
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
                PetscCallAbort(comm, ierr);
              }
            }
          }
        }

        ierr = VecAssemblyBegin(res);
        PetscCallAbort(comm, ierr);
        ierr = VecAssemblyEnd(res);
        PetscCallAbort(comm, ierr);
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  // MPI assembly for PETSc Mat (bilinear form) -- skipping ghosts via Shard + Boost.MPI
  template <class TrialFES, class TestFES>
  class MPI<
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
      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "DotType must be PetscScalar"
      );
      using OperatorType     = ::Mat&;
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
        const PetscInt localRows = PetscInt(testFES.getSize());
        const PetscInt localCols = PetscInt(trialFES.getSize());
        const PetscInt globalRows = boost::mpi::all_reduce(comm, localRows, std::plus<PetscInt>());
        const PetscInt globalCols = boost::mpi::all_reduce(comm, localCols, std::plus<PetscInt>());

        // Create/init Mat
        ierr = MatCreate(comm, &A);
        PetscCallAbort(comm, ierr);
        ierr = MatSetSizes(A, localRows, localCols, globalRows, globalCols);
        PetscCallAbort(comm, ierr);
        ierr = MatSetFromOptions(A);
        PetscCallAbort(comm, ierr);
        ierr = MatMPIAIJSetPreallocation(A,
                                         /*d_nz*/PETSC_DECIDE, nullptr,
                                         /*o_nz*/PETSC_DECIDE, nullptr);
        PetscCallAbort(comm, ierr);
        ierr = MatSetUp(A);
        PetscCallAbort(comm, ierr);
        ierr = MatZeroEntries(A);
        PetscCallAbort(comm, ierr);

        // Local element contributions skipping ghosts
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            auto d = it->getDimension();
            auto i = it->getIndex();
            if (shard.isGhost(d, i))
              continue;

            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = testFES .getDOFs(d, i);
              const auto& cols = trialFES.getDOFs(d, i);
              for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
                for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                {
                  PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
                  PetscCallAbort(comm, ierr);
                }
            }
          }
        }

        // Global coupling contributions skipping ghosts
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& testAttrs  = bfi.getTestAttributes();
          const auto& trialAttrs = bfi.getTrialAttributes();
          Internal::SequentialIteration testSeq (mesh, bfi.getTestRegion());
          Internal::SequentialIteration trialSeq(mesh, bfi.getTrialRegion());
          for (auto teIt = testSeq.getIterator(); teIt; ++teIt)
          {
            auto dt = teIt->getDimension();
            auto tg = teIt->getIndex();
            if (shard.isGhost(dt, tg))
              continue;

            const auto& tmap = shard.getPolytopeMap(dt).left;
            auto tLoc = tmap.find(tg);
            if (tLoc == tmap.end())
              continue;
            auto tLid = tLoc->second;

            if (testAttrs.empty() || testAttrs.count(teIt->getAttribute()))
            {
              for (auto trIt = trialSeq.getIterator(); trIt; ++trIt)
              {
                auto dr = trIt->getDimension();
                auto rg = trIt->getIndex();
                if (shard.isGhost(dr, rg))
                  continue;

                const auto& rmap = shard.getPolytopeMap(dr).left;
                auto rLoc = rmap.find(rg);
                if (rLoc == rmap.end())
                  continue;
                auto rLid = rLoc->second;

                if (trialAttrs.empty() || trialAttrs.count(trIt->getAttribute()))
                {
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = testFES .getDOFs(dt, tLid);
                  const auto& cols = trialFES.getDOFs(dr, rLid);
                  for (PetscInt i = 0; i < PetscInt(rows.size()); ++i)
                    for (PetscInt j = 0; j < PetscInt(cols.size()); ++j)
                    {
                      const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                      ierr = MatSetValue(A, rows[i], cols[j], v, ADD_VALUES);
                      PetscCallAbort(comm, ierr);
                    }
                }
              }
            }
          }
        }

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(comm, ierr);
        ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(comm, ierr);
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

#endif // RODIN_PETSC_ASSEMBLY_MPI_H
