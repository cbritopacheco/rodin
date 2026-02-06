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
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
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

        ierr = MatMPIAIJSetPreallocation(
            res, PETSC_DECIDE, PETSC_NULLPTR, PETSC_DECIDE, PETSC_NULLPTR);
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
            if (attrs.empty() || attrs.count(it->getAttribute()))
            {
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

  template <class TrialFunction, class TestFunction>
  class MPI<Rodin::PETSc::Math::LinearSystem, Variational::Problem<Rodin::PETSc::Math::LinearSystem, TrialFunction, TestFunction>>
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

      MPI() = default;

      MPI(const MPI& other)
        : Parent(other)
      {}

      MPI(MPI&& other)
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
