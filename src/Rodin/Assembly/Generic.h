/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_GENERIC_H
#define RODIN_ASSEMBLY_GENERIC_H

#include "Rodin/Variational/ForwardDecls.h"

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Generic<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      using LinearSystemType =
        LinearSystem;

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
        AssemblyBase<LinearSystemType, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>;

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
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES = v.getFiniteElementSpace();
        auto& pb = input.getProblemBody();
        auto& [stiffness, solution, mass] = out;

        const BilinearFormAssemblyType bfa;
        bfa.execute(
            stiffness,
            {
              u.getFiniteElementSpace(),
              v.getFiniteElementSpace(),
              pb.getLocalBFIs(),
              pb.getGlobalBFIs()
            });

        for (auto& bf : pb.getBFs())
          stiffness += bf.getOperator();

        const LinearFormAssemblyType lfa;
        lfa.execute(
            mass,
            {
              v.getFiniteElementSpace(),
              pb.getLFIs()
            });

        for (const auto& lf : pb.getLFs())
          mass += lf.getVector();

        mass *= -1;

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

#endif
