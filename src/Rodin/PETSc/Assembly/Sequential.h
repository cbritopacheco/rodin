#ifndef RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H
#define RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <optional>

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
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
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

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
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

        // Global contributions
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }
            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }
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

  // Sequential assembly for single-variable Problem (PETSc)
  template <class U, class V>
  class Sequential<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>>
  {
    public:
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U, V>;
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;
      using InputType = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec
      using ScalarType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::ScalarType;   // PetscScalar
      using TrialFESType = typename Rodin::FormLanguage::Traits<U>::FESType;
      using TrialMeshType = typename Rodin::FormLanguage::Traits<TrialFESType>::MeshType;
      using TrialMeshContextType = typename Rodin::FormLanguage::Traits<TrialMeshType>::ContextType;

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        static_assert(std::is_same_v<TrialMeshContextType, Rodin::Context::Local>,
          "PETSc sequential assembly should only be used with Local mesh context.");

        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();

        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES  = v.getFiniteElementSpace();
        const auto& mesh     = trialFES.getMesh();

        const size_t cols = static_cast<size_t>(trialFES.getSize());
        const size_t rows = static_cast<size_t>(testFES.getSize());

        PetscErrorCode ierr;

        // Matrix setup
        assert(A);
        ierr = MatSetSizes(A,
                           static_cast<PetscInt>(rows),
                           static_cast<PetscInt>(cols),
                           PETSC_DETERMINE,
                           PETSC_DETERMINE);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetType(A, MATSEQAIJ);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // Vector setup
        assert(b);
        ierr = VecSetSizes(b, static_cast<PetscInt>(rows), PETSC_DECIDE);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetType(b, VECSEQ);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(x, static_cast<PetscInt>(cols), PETSC_DECIDE);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetType(x, VECSEQ);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

        // Local BFIs
        using MeshType = std::decay_t<decltype(mesh)>;

        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration<MeshType> seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            const size_t d = it->getDimension();
            const Index  p = it->getIndex();

            bfi.setPolytope(*it);

            const auto& rowsDOF = testFES.getDOFs(d, p);
            const auto& colsDOF = trialFES.getDOFs(d, p);

            for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
            {
              for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
              {
                const PetscScalar val = PetscScalar(bfi.integrate(j, i));
                if (val != PetscScalar(0))
                {
                  ierr = MatSetValue(A, rowsDOF[i], colsDOF[j], val, ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
              }
            }
          }
        }

        // Global BFIs
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration<MeshType> trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration<MeshType> testseq(mesh, bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }

            const auto& rowsDOF = testFES.getDOFs(teIt->getDimension(), teIt->getIndex());

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }

              const auto& colsDOF = trialFES.getDOFs(trIt->getDimension(), trIt->getIndex());

              bfi.setPolytope(*trIt, *teIt);

              for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
              {
                for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
                {
                  const PetscScalar val = PetscScalar(bfi.integrate(j, i));
                  if (val != PetscScalar(0))
                  {
                    ierr = MatSetValue(A, rowsDOF[i], colsDOF[j], val, ADD_VALUES);
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

        // Linear forms
        for (auto& lfi : pb.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration<MeshType> seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            lfi.setPolytope(*it);
            const auto& dofs = testFES.getDOFs(it->getDimension(), it->getIndex());
            for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
            {
              const PetscScalar val = PetscScalar(lfi.integrate(l));
              ierr = VecSetValue(b, dofs[l], -val, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
            }
          }
        }

        // Preassembled linear forms
        for (auto& lf : pb.getLFs())
        {
          ierr = VecAXPY(b, -1.0, lf.getVector());
          assert(ierr == PETSC_SUCCESS);
        }

        // Assemble A and b
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // Dirichlet BCs
        std::vector<PetscInt> bcIdx;
        std::vector<PetscScalar> bcVals;
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;
          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
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

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class U1, class U2, class U3, class ... Us>
  class Sequential<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      using LinearSystemType =
        Rodin::PETSc::Math::LinearSystem;

      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      using Parent =
        AssemblyBase<LinearSystemType, ProblemType>;

      using InputType = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        auto&       us           = input.getTrialFunctions();
        auto&       vs           = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto&       trialUUIDMap = input.getTrialUUIDMap();
        auto&       testUUIDMap  = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        using FirstTrialMeshType =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        using MeshContextType =
          typename Rodin::FormLanguage::Traits<FirstTrialMeshType>::ContextType;

        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::Local>,
          "PETSc sequential assembly supports only Local mesh context.");

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (SeqAIJ)
        // ------------------------
        assert(A);
        ierr = MatSetSizes(A,
                           static_cast<PetscInt>(nrows),
                           static_cast<PetscInt>(ncols),
                           PETSC_DETERMINE,
                           PETSC_DETERMINE);
        assert(ierr == PETSC_SUCCESS);

        // Keep it explicitly AIJ (sequential) for now
        ierr = MatSetType(A, MATSEQAIJ);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

        // ------------------------
        // Allocate / reset b (Seq Vec)
        // ------------------------
        assert(b);
        ierr = VecSetSizes(b, static_cast<PetscInt>(nrows), PETSC_DECIDE);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetType(b, VECSEQ);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(x, static_cast<PetscInt>(ncols), PETSC_DECIDE);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetType(x, VECSEQ);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Helpers
        // ------------------------
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

        const auto withTrialFES = [&](const auto& uuid, auto&& fn)
        {
          const size_t k = findTrialBlock(uuid);
          bool found = false;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
            {
              fn(uref.get().getFiniteElementSpace());
              found = true;
            }
          });
          assert(found);
        };

        const auto withTestFES = [&](const auto& uuid, auto&& fn)
        {
          const size_t k = findTestBlock(uuid);
          bool found = false;
          vs.iapply([&](size_t i, const auto& vref)
          {
            if (i == k)
            {
              fn(vref.get().getFiniteElementSpace());
              found = true;
            }
          });
          assert(found);
        };

        using MeshType0 =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        std::optional<std::reference_wrapper<const MeshType0>> meshRef;
        us.apply([&](const auto& uref)
        {
          using MeshT = std::decay_t<decltype(uref.get().getFiniteElementSpace().getMesh())>;
          static_assert(std::is_same_v<MeshT, MeshType0>,
            "Mixed mesh types are not supported in PETSc multi-field sequential assembly.");
          if (!meshRef)
            meshRef = std::cref(uref.get().getFiniteElementSpace().getMesh());
        });
        assert(meshRef.has_value());
        const MeshType0& mesh = meshRef->get();

        // ------------------------
        // Assemble bilinear terms into A
        // ------------------------
        // Local BFIs
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());

          withTrialFES(uUUID, [&](const auto& uFES)
          {
            withTestFES(vUUID, [&](const auto& vFES)
            {
              // Reusable buffers to avoid realloc each cell
              std::vector<PetscInt>    Iidx;
              std::vector<PetscInt>    Jidx;
              std::vector<PetscScalar> Ke; // row-major (i * nc + j)

              for (auto it = seq.getIterator(); it; ++it)
              {
                if (!attrs.empty())
                {
                  const auto a = it->getAttribute();
                  if (!a || !attrs.count(*a))
                    continue;
                }

                const size_t d = it->getDimension();
                const Index  p = it->getIndex();

                bfi.setPolytope(*it);

                const auto& rows = vFES.getDOFs(d, p);
                const auto& cols = uFES.getDOFs(d, p);

                const PetscInt nr = static_cast<PetscInt>(rows.size());
                const PetscInt nc = static_cast<PetscInt>(cols.size());

                Iidx.resize(static_cast<size_t>(nr));
                Jidx.resize(static_cast<size_t>(nc));
                Ke.resize(static_cast<size_t>(nr) * static_cast<size_t>(nc));

                // Map DOFs -> global indices once
                for (PetscInt i = 0; i < nr; ++i)
                {
                  Iidx[static_cast<size_t>(i)] =
                    static_cast<PetscInt>(vOff + static_cast<size_t>(rows[static_cast<size_t>(i)]));
                }

                for (PetscInt j = 0; j < nc; ++j)
                {
                  Jidx[static_cast<size_t>(j)] =
                    static_cast<PetscInt>(uOff + static_cast<size_t>(cols[static_cast<size_t>(j)]));
                }

                // Build element matrix (row-major)
                for (PetscInt i = 0; i < nr; ++i)
                {
                  for (PetscInt j = 0; j < nc; ++j)
                  {
                    Ke[static_cast<size_t>(i) * static_cast<size_t>(nc) + static_cast<size_t>(j)] =
                      static_cast<PetscScalar>(bfi.integrate(static_cast<size_t>(j), static_cast<size_t>(i)));
                  }
                }

                // One PETSc call per element (instead of nr*nc calls)
                ierr = MatSetValues(
                  A,
                  nr, Iidx.data(),
                  nc, Jidx.data(),
                  Ke.data(),
                  ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            });
          });
        }

        // Global BFIs
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          withTrialFES(uUUID, [&](const auto& uFES)
          {
            withTestFES(vUUID, [&](const auto& vFES)
            {
              for (auto teIt = testseq.getIterator(); teIt; ++teIt)
              {
                if (!testAttrs.empty())
                {
                  const auto a = teIt->getAttribute();
                  if (!a || !testAttrs.count(*a))
                    continue;
                }

                const auto& rows = vFES.getDOFs(teIt.getDimension(), teIt->getIndex());

                for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                {
                  if (!trialAttrs.empty())
                  {
                    const auto a = trIt->getAttribute();
                    if (!a || !trialAttrs.count(*a))
                      continue;
                  }

                  const auto& cols = uFES.getDOFs(trIt.getDimension(), trIt->getIndex());

                  bfi.setPolytope(*trIt, *teIt);

                  for (PetscInt i = 0; i < static_cast<PetscInt>(rows.size()); ++i)
                  {
                    const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(rows[i]));
                    for (PetscInt j = 0; j < static_cast<PetscInt>(cols.size()); ++j)
                    {
                      const PetscInt J = static_cast<PetscInt>(uOff + static_cast<size_t>(cols[j]));
                      const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(j, i));
                      ierr = MatSetValue(A, I, J, val, ADD_VALUES);
                      assert(ierr == PETSC_SUCCESS);
                    }
                  }
                }
              }
            });
          });
        }

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Assemble linear terms into b
        // ------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());

          withTestFES(vUUID, [&](const auto& vFES)
          {
            for (auto it = seq.getIterator(); it; ++it)
            {
              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);

              const auto& dofs = vFES.getDOFs(it.getDimension(), it->getIndex());
              for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
              {
                const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(dofs[l]));
                const PetscScalar val = static_cast<PetscScalar>(lfi.integrate(l));
                ierr = VecSetValue(b, I, -val, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          });
        }

        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Impose Dirichlet BCs via MatZeroRowsColumns
        // ------------------------
        std::vector<PetscInt> bcIdx;
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
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using Sequential = Rodin::Assembly::Sequential<LinearAlgebraType, Operand>;
}

#endif // RODIN_ASSEMBLY_SEQUENTIAL_PETSC_H
