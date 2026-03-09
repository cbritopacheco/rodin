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

#include "Rodin/PETSc/Math/LinearSystem.h"

namespace Rodin::Assembly
{
  template <class FES>
  class MPI<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
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
        assert(res);
        const auto& fes = input.getFES();
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        const size_t globalSize = fes.getSize();

        size_t begin, end;
        fes.getOwnershipRange(begin, end);
        const size_t ownedSize = end - begin;

        PetscErrorCode ierr;
        ierr = VecSetSizes(res, ownedSize, globalSize);
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
            const Index idx = it->getIndex();
            if (shard.isGhost(d, idx))
              continue;
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            lfi.setPolytope(*it);
            const auto& dofs = fes.getShard().getDOFs(d, idx);
            for (PetscInt i = 0; i < dofs.size(); ++i)
            {
              const Index r = fes.getGlobalIndex(dofs[i]);
              const PetscScalar v = lfi.integrate(i);
              ierr = VecSetValue(res, r, v, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
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
  class MPI<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
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
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;
      using Parent           = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType        = typename Parent::InputType;

      void execute(OperatorType& res, const InputType& input) const override
      {
        assert(res);

        const auto& trialFES = input.getTrialFES();
        const auto& testFES  = input.getTestFES();
        assert(trialFES.getMesh() == testFES.getMesh());
        const auto& mesh = trialFES.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();

        size_t rbegin, rend;
        testFES.getOwnershipRange(rbegin, rend);
        const size_t localRows = rend - rbegin;

        size_t cbegin, cend;
        trialFES.getOwnershipRange(cbegin, cend);
        const size_t localCols = cend - cbegin;

        const size_t globalRows = testFES.getSize();

        const size_t globalCols = trialFES.getSize();

        PetscErrorCode ierr;
        ierr = MatSetSizes(res, localRows, localCols, globalRows, globalCols);
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
            const Index idx = it->getIndex();
            if (shard.isGhost(d, idx))
              continue;
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            bfi.setPolytope(*it);
            const auto& rows = testFES.getShard().getDOFs(d, idx);
            const auto& cols = trialFES.getShard().getDOFs(d, idx);
            for (Index i = 0; i < rows.size(); ++i)
            {
              const Index r = testFES.getGlobalIndex(rows[i]);
              for (Index j = 0; j < cols.size(); ++j)
              {
                const Index c = trialFES.getGlobalIndex(cols[j]);
                PetscScalar v = bfi.integrate(j, i);
                ierr = MatSetValue(res, r, c, v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          }
        }

        ierr = MatAssemblyBegin(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatAssemblyEnd(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        (void) ierr;
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  // ------------------------------------------------------------
  // MPI assembly for single-variable Problem (PETSc)
  // ------------------------------------------------------------
  template <class U, class V>
  class MPI<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>>
  {
    public:
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      using ProblemType      = Rodin::Variational::Problem<LinearSystemType, U, V>;
      using Parent           = AssemblyBase<LinearSystemType, ProblemType>;
      using InputType        = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      using TrialFESType  = typename Rodin::FormLanguage::Traits<U>::FESType;
      using TrialMeshType = typename Rodin::FormLanguage::Traits<TrialFESType>::MeshType;
      using MeshContextType =
        typename Rodin::FormLanguage::Traits<TrialMeshType>::ContextType;

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::MPI>,
          "PETSc MPI assembly should only be used with MPI mesh context.");

        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();

        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES  = v.getFiniteElementSpace();
        assert(trialFES.getMesh() == testFES.getMesh());

        const auto& mesh  = trialFES.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();

        const size_t globalCols = trialFES.getSize();
        const size_t globalRows = testFES.getSize();

        // Ownership ranges (global indices)
        size_t rbegin, rend;
        testFES.getOwnershipRange(rbegin, rend);
        const size_t localRows = rend - rbegin;

        size_t cbegin, cend;
        trialFES.getOwnershipRange(cbegin, cend);
        const size_t localCols = cend - cbegin;

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (MPIAIJ)
        // ------------------------
        assert(A);
        ierr = MatSetSizes(
            A,
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(globalRows),
            static_cast<PetscInt>(globalCols));
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetFromOptions(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Allocate / reset b (MPI Vec)
        // ------------------------
        assert(b);
        ierr = VecSetSizes(
            b,
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(globalRows));
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(
            x,
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(globalCols));
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Local BFIs (owned elements only)
        // ------------------------
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MPIIteration seq(mesh, bfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d   = it->getDimension();
            const Index  idx = it->getIndex();

            if (shard.isGhost(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            bfi.setPolytope(*it);

            const auto& rowsDOF = testFES.getShard().getDOFs(d, idx);
            const auto& colsDOF = trialFES.getShard().getDOFs(d, idx);

            for (Index i = 0; i < rowsDOF.size(); ++i)
            {
              const PetscInt I = static_cast<PetscInt>(testFES.getGlobalIndex(rowsDOF[i]));
              for (Index j = 0; j < colsDOF.size(); ++j)
              {
                const PetscInt J = static_cast<PetscInt>(trialFES.getGlobalIndex(colsDOF[j]));
                const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(j, i));
                if (val != PetscScalar(0))
                {
                  ierr = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
              }
            }
          }
        }

        // ------------------------
        // Global BFIs (owned test entities only)
        // - This avoids double assembly across ranks.
        // - We allow trial entities to be ghost; PETSc handles off-proc columns.
        // ------------------------
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          MPIIteration trialseq(mesh, bfi.getTrialRegion());
          MPIIteration testseq(mesh,  bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            const size_t td   = teIt->getDimension();
            const Index  tidx = teIt->getIndex();

            if (shard.isGhost(td, tidx))
              continue;

            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }

            const auto& rowsDOF = testFES.getShard().getDOFs(td, tidx);

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              const size_t rd   = trIt->getDimension();
              const Index  ridx = trIt->getIndex();

              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }

              // NOTE: do NOT skip ghost trial entities here: they still contribute to
              // owned test rows, producing off-process columns in MatSetValue.
              const auto& colsDOF = trialFES.getShard().getDOFs(rd, ridx);

              bfi.setPolytope(*trIt, *teIt);

              for (Index i = 0; i < rowsDOF.size(); ++i)
              {
                const PetscInt I = static_cast<PetscInt>(testFES.getGlobalIndex(rowsDOF[i]));
                for (Index j = 0; j < colsDOF.size(); ++j)
                {
                  const PetscInt J = static_cast<PetscInt>(trialFES.getGlobalIndex(colsDOF[j]));
                  const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(j, i));
                  if (val != PetscScalar(0))
                  {
                    ierr = MatSetValue(A, I, J, val, ADD_VALUES);
                    assert(ierr == PETSC_SUCCESS);
                  }
                }
              }
            }
          }
        }

        // Preassembled bilinear forms
        for (auto& bf : pb.getBFs())
        {
          const auto& op = bf.getOperator();
          ierr = MatAXPY(A, 1.0, op, DIFFERENT_NONZERO_PATTERN);
          assert(ierr == PETSC_SUCCESS);
        }

        // Assemble A
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Linear forms (owned elements only)
        // Convention: matches your sequential single-variable version: b += -LF
        // ------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          MPIIteration seq(mesh, lfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d   = it->getDimension();
            const Index  idx = it->getIndex();

            if (shard.isGhost(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            lfi.setPolytope(*it);

            const auto& dofs = testFES.getShard().getDOFs(d, idx);
            for (Index l = 0; l < dofs.size(); ++l)
            {
              const PetscInt I = static_cast<PetscInt>(testFES.getGlobalIndex(dofs[l]));
              const PetscScalar val = static_cast<PetscScalar>(lfi.integrate(static_cast<PetscInt>(l)));
              if (val != PetscScalar(0))
              {
                ierr = VecSetValue(b, I, -val, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          }
        }

        // Preassembled linear forms: b += -LF
        for (auto& lf : pb.getLFs())
        {
          ierr = VecAXPY(b, -1.0, lf.getVector());
          assert(ierr == PETSC_SUCCESS);
        }

        // Assemble b
        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Dirichlet BCs (only for this single trial u)
        // - Use global indices; safe to provide all ranks' bcIdx (PETSc will route).
        // - If dbc.getDOFs() is already global, drop the getGlobalIndex conversion.
        // ------------------------
        std::vector<PetscInt>    bcIdx;
        std::vector<PetscScalar> bcVals;

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
            // local is assumed in the trial FES numbering (global in the FES sense).
            // In your sequential code local already matched the global row index.
            // Here, keep the same interpretation:
            bcIdx.push_back(static_cast<PetscInt>(local));
            bcVals.push_back(static_cast<PetscScalar>(value));
          }
        }

        if (!bcIdx.empty())
        {
          Vec bcVec;
          ierr = VecDuplicate(b, &bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetValues(
              bcVec,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecAssemblyBegin(bcVec);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecAssemblyEnd(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = MatZeroRowsColumns(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              1.0,
              bcVec,
              b);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecDestroy(&bcVec);
          assert(ierr == PETSC_SUCCESS);
        }
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };


  // ------------------------------------------------------------
  // MPI assembly for multi-variable Problem (PETSc)
  // ------------------------------------------------------------
  template <class U1, class U2, class U3, class ... Us>
  class MPI<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;
      using Parent    = AssemblyBase<LinearSystemType, ProblemType>;
      using InputType = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb            = input.getProblemBody();
        auto& us            = input.getTrialFunctions();
        auto& vs            = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto& trialUUIDMap  = input.getTrialUUIDMap();
        auto& testUUIDMap   = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        using FirstTrialMeshType =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        using MeshContextType =
          typename Rodin::FormLanguage::Traits<FirstTrialMeshType>::ContextType;

        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::MPI>,
          "PETSc MPI assembly should only be used with MPI mesh context.");

        // Mesh + shard
        const auto& mesh = [&]() -> const auto&
        {
          const void* addr = nullptr;
          us.apply([&](const auto& uref)
          {
            if (!addr)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace().getMesh());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>*>(addr);
        }();

        const auto& shard = mesh.getShard();

        // Ownership range for the *global* test vector (assembled rows)
        // Here we assume "ownership range" refers to the assembled global test numbering (with offsets applied),
        // hence we use input.getTotalTestSize() and distribute it by PETSc using Vec/Mat layout already set by caller,
        // OR you keep FES-based ownership and accept that blocks can be owned by different ranks.
        // In your code base you already have FES.getOwnershipRange for each field; we follow that and let PETSc handle.
        //
        // Layout: localRows/localCols computed as sum of owned sizes per block.

        // Helpers (same as your sequential)
        const auto findTrialBlock = [&](const auto& uuid) -> size_t
        {
          auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return it->second;
        };

        const auto findTestBlock = [&](const auto& uuid) -> size_t
        {
          auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return it->second;
        };

        const auto& getTrialFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTrialBlock(uuid);
          const void* addr = nullptr;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        const auto& getTestFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTestBlock(uuid);
          const void* addr = nullptr;
          vs.iapply([&](size_t i, const auto& vref)
          {
            if (i == k)
              addr = static_cast<const void*>(&vref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(vs.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        // Compute localRows/localCols as sum of owned sizes of each block on this rank
        size_t localRows = 0;
        vs.iapply([&](size_t /*i*/, const auto& vref)
        {
          size_t rb, re;
          vref.get().getFiniteElementSpace().getOwnershipRange(rb, re);
          localRows += (re - rb);
        });

        size_t localCols = 0;
        us.iapply([&](size_t /*i*/, const auto& uref)
        {
          size_t cb, ce;
          uref.get().getFiniteElementSpace().getOwnershipRange(cb, ce);
          localCols += (ce - cb);
        });

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (MPIAIJ)
        // ------------------------
        assert(A);
        ierr = MatSetSizes(
            A,
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(nrows),
            static_cast<PetscInt>(ncols));
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetFromOptions(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Allocate / reset b (MPI Vec)
        // ------------------------
        assert(b);
        ierr = VecSetSizes(
            b,
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(nrows));
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(
            x,
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(ncols));
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Assemble bilinear terms into A
        // - For each integrator, only assemble contributions associated to an owned "assembly entity":
        //   - local BFIs: owned cell/facet/etc (same as your MPI Mat/BFI specialization)
        //   - global BFIs: owned test entity (rows owned), trial may be ghost
        // ------------------------

        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& uFES = getTrialFESByUUID(uUUID);
          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = bfi.getAttributes();
          MPIIteration seq(mesh, bfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d   = it->getDimension();
            const Index  idx = it->getIndex();

            if (shard.isGhost(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            bfi.setPolytope(*it);

            const auto& rows = vFES.getShard().getDOFs(d, idx);
            const auto& cols = uFES.getShard().getDOFs(d, idx);

            for (Index i = 0; i < rows.size(); ++i)
            {
              const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(vFES.getGlobalIndex(rows[i])));
              for (Index j = 0; j < cols.size(); ++j)
              {
                const PetscInt J = static_cast<PetscInt>(uOff + static_cast<size_t>(uFES.getGlobalIndex(cols[j])));
                const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(static_cast<PetscInt>(j),
                                                                              static_cast<PetscInt>(i)));
                if (val != PetscScalar(0))
                {
                  ierr = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
              }
            }
          }
        }

        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& uFES = getTrialFESByUUID(uUUID);
          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          MPIIteration trialseq(mesh, bfi.getTrialRegion());
          MPIIteration testseq(mesh,  bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            const size_t td   = teIt->getDimension();
            const Index  tidx = teIt->getIndex();

            if (shard.isGhost(td, tidx))
              continue;

            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }

            const auto& rows = vFES.getShard().getDOFs(td, tidx);

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              const size_t rd   = trIt->getDimension();
              const Index  ridx = trIt->getIndex();

              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }

              // do not skip ghost trial entity: off-proc columns are fine
              const auto& cols = uFES.getShard().getDOFs(rd, ridx);

              bfi.setPolytope(*trIt, *teIt);

              for (Index i = 0; i < rows.size(); ++i)
              {
                const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(vFES.getGlobalIndex(rows[i])));
                for (Index j = 0; j < cols.size(); ++j)
                {
                  const PetscInt J = static_cast<PetscInt>(uOff + static_cast<size_t>(uFES.getGlobalIndex(cols[j])));
                  const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(static_cast<PetscInt>(j),
                                                                                static_cast<PetscInt>(i)));
                  if (val != PetscScalar(0))
                  {
                    ierr = MatSetValue(A, I, J, val, ADD_VALUES);
                    assert(ierr == PETSC_SUCCESS);
                  }
                }
              }
            }
          }
        }

        // Assemble A
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Assemble linear terms into b
        // Convention: matches your sequential multi-variable version (no minus).
        // If you want the single-variable convention (-LF), swap sign here.
        // ------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = lfi.getAttributes();
          MPIIteration seq(mesh, lfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d   = it->getDimension();
            const Index  idx = it->getIndex();

            if (shard.isGhost(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            lfi.setPolytope(*it);

            const auto& dofs = vFES.getShard().getDOFs(d, idx);
            for (Index l = 0; l < dofs.size(); ++l)
            {
              const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(vFES.getGlobalIndex(dofs[l])));
              const PetscScalar val = static_cast<PetscScalar>(lfi.integrate(static_cast<PetscInt>(l)));
              if (val != PetscScalar(0))
              {
                ierr = VecSetValue(b, I, -val, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          }
        }

        // Assemble b
        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Dirichlet BCs via MatZeroRowsColumns (serial on each rank, PETSc routes)
        // ------------------------
        std::vector<PetscInt>    bcIdx;
        std::vector<PetscScalar> bcVals;

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
            bcIdx.push_back(static_cast<PetscInt>(uOff + local));
            bcVals.push_back(static_cast<PetscScalar>(value));
          }
        }

        if (!bcIdx.empty())
        {
          Vec bcVec;
          ierr = VecDuplicate(b, &bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetValues(
              bcVec,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecAssemblyBegin(bcVec);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecAssemblyEnd(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = MatZeroRowsColumns(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              1.0,
              bcVec,
              b);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecDestroy(&bcVec);
          assert(ierr == PETSC_SUCCESS);
        }

        (void) ierr;
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
