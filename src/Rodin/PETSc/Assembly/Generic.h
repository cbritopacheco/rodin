/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_GENERIC_H
#define RODIN_PETSC_ASSEMBLY_GENERIC_H

#include "Rodin/Assembly/ForwardDecls.h"

#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

namespace Rodin::Assembly
{
  template <class TrialFunction, class TestFunction>
  class Generic<
    PETSc::Math::LinearSystem, Variational::Problem<
      PETSc::Math::LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<
        PETSc::Math::LinearSystem,
        Variational::Problem<PETSc::Math::LinearSystem, TrialFunction, TestFunction>>
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

      Generic() = default;

      Generic(const Generic& other)
        : Parent(other)
      {}

      Generic(Generic&& other)
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

      Generic* copy() const noexcept override
      {
        return new Generic(*this);
      }
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using Generic = Rodin::Assembly::Generic<LinearAlgebraType, Operand>;
}

#endif
