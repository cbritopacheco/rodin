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

  template <class TrialFunction, class TestFunction>
  class Sequential<Rodin::PETSc::Math::LinearSystem, Variational::Problem<Rodin::PETSc::Math::LinearSystem, TrialFunction, TestFunction>>
    : public AssemblyBase<Rodin::PETSc::Math::LinearSystem, Variational::Problem<Rodin::PETSc::Math::LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      using LinearSystemType =
        PETSc::Math::LinearSystem;

      using TrialFESType =
        typename FormLanguage::Traits<TrialFunction>::FESType;

      using TestFESType =
        typename FormLanguage::Traits<TestFunction>::FESType;

      using TrialFESMeshType =
        typename FormLanguage::Traits<TrialFESType>::MeshType;

      using TestFESMeshType =
        typename FormLanguage::Traits<TestFESType>::MeshType;

      using TrialFESMeshContextType =
        typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      using TestFESMeshContextType =
        typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      static_assert(std::is_same_v<TrialFESMeshContextType, TestFESMeshContextType>);

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using SolutionType =
        typename FormLanguage::Traits<TrialFunction>::SolutionType;

      using BilinearFormType =
        Variational::BilinearForm<SolutionType, TrialFESType, TestFESType, OperatorType>;

      using BilinearFormAssemblyType =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<OperatorType, BilinearFormType>;

      using LinearFormType =
        Variational::LinearForm<TestFESType, VectorType>;

      using LinearFormAssemblyType =
        typename Assembly::Default<TestFESMeshContextType>
          ::template Type<VectorType, LinearFormType>;

      using Parent =
        AssemblyBase<LinearSystemType, Variational::Problem<LinearSystemType, TrialFunction, TestFunction>>;

      using InputType =
        typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      void execute(LinearSystemType& out, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES = v.getFiniteElementSpace();
        auto& pb = input.getProblemBody();
        auto& [stiffness, solution, mass] = out;
        const MPI_Comm comm = out.getCommunicator();

        const BilinearFormAssemblyType bfa;
        assert(stiffness);
        bfa.execute(
            stiffness,
            {
              u.getFiniteElementSpace(),
              v.getFiniteElementSpace(),
              pb.getLocalBFIs(),
              pb.getGlobalBFIs()
            });

        for (auto& bf : pb.getBFs())
        {
          ierr = MatAXPY(stiffness, 1.0, bf.getOperator(), UNKNOWN_NONZERO_PATTERN);
          assert(ierr == PETSC_SUCCESS);
        }

        const LinearFormAssemblyType lfa;
        assert(mass);
        lfa.execute(
            mass,
            {
              v.getFiniteElementSpace(),
              pb.getLFIs()
            });

        for (const auto& lf : pb.getLFs())
        {
          ierr = VecAXPY(mass, 1.0, lf.getVector());
          assert(ierr == PETSC_SUCCESS);
        }

        ierr = VecScale(mass, -1.0);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecCreate(comm, &solution);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecDuplicate(mass, &solution);
        assert(ierr == PETSC_SUCCESS);

        assert(solution);
        for (auto& dbc : pb.getDBCs())
        {
          dbc.assemble();
          out.eliminate(dbc.getDOFs());
        }

        for (auto& pbc : pb.getPBCs())
        {
          assert(trialFES == testFES);
          pbc.assemble();
          out.merge(pbc.getDOFs());
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

        // Cheap preallocation (works; not optimal). You can improve with nnz estimates later.
        ierr = MatSeqAIJSetPreallocation(A, PETSC_DETERMINE, PETSC_NULLPTR);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

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

        // Retrieve a reference to the FES by UUID (trial)
        const auto& getTrialFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTrialBlock(uuid);
          // us is a Tuple of reference_wrappers; index access is via iapply/map in your Tuple
          // We implement a scan to avoid relying on tuple operator[].
          const void* addr = nullptr;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        // Retrieve a reference to the FES by UUID (test)
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

        // Mesh: take it from the first trial FES
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

          const auto& uFES = getTrialFESByUUID(uUUID);
          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            const size_t d = it->getDimension();
            const Index  p = it->getIndex();

            bfi.setPolytope(*it);

            const auto& rows = vFES.getDOFs(d, p);
            const auto& cols = uFES.getDOFs(d, p);

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

        // Global BFIs
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

          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty() && !testAttrs.count(teIt->getAttribute()))
              continue;

            const auto& rows = vFES.getDOFs(teIt.getDimension(), teIt->getIndex());

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty() && !trialAttrs.count(trIt->getAttribute()))
                continue;

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

          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            lfi.setPolytope(*it);

            const auto& dofs = vFES.getDOFs(it.getDimension(), it->getIndex());
            for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
            {
              const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(dofs[l]));
              const PetscScalar val = static_cast<PetscScalar>(lfi.integrate(l));
              ierr = VecSetValue(b, I, val, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
            }
          }
        }

        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Impose Dirichlet BCs (delegated to LinearSystem eliminate)
        // ------------------------
        // Assumes:
        //   - dbc.getOperand() is a TrialFunction-like object with UUID
        //   - dbc.assemble() populates dbc.getDOFs()
        //   - axb.eliminate(dofs, offset) exists and acts on A and b
        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          axb.eliminate(dofs, uOff);
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

