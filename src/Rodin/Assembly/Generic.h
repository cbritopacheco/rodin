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
      using LinearSystemType = LinearSystem;

      using TrialFESType =
        typename FormLanguage::Traits<TrialFunction>::FESType;

      using TestFESType =
        typename FormLanguage::Traits<TestFunction>::FESType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using SolutionType =
        typename FormLanguage::Traits<TrialFunction>::SolutionType;

      using LinearFormType =
        Variational::LinearForm<TestFESType, VectorType>;

      using DirichletBCType =
        Variational::DirichletBC<TrialFunction, VectorType>;

      using BilinearFormType =
        Variational::BilinearForm<SolutionType, TrialFESType, TestFESType, OperatorType>;

      using LinearFormIntegratorBaseListType =
        FormLanguage::List<Variational::LinearFormIntegratorBase<ScalarType>>;

      using Parent =
        AssemblyBase<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>;

      using InputType =
        typename Parent::InputType;

      Generic() = default;

      Generic(const Generic& other)
        : Parent(other)
      {}

      Generic(Generic&& other)
        : Parent(std::move(other))
      {}

      void execute(LinearSystem& out, const InputType& input) const override
      {
        // const auto& u = input.getTrialFunction();
        // const auto& v = input.getTestFunction();
        // const auto& trialFES = u.getFiniteElementSpace();
        // const auto& testFES = v.getFiniteElementSpace();
        // auto& pb = input.getProblemBody();
        // auto& [stiffness, solution, mass] = out;

        // m_bfa.execute(stiffness, { u.getFiniteElementSpace(), v.getFiniteElementSpace(), pb.getLocalBFIs(), pb.getGlobalBFIs() });
        // for (auto& bf : pb.getBFs())
        //   stiffness += bf.getOperator();

        // m_lfa.execute(mass, { v.getFiniteElementSpace(), pb.getLFIs() });
        // mass *= -1;

        // for (auto& dbc : pb.getDBCs())
        // {
        //   dbc.assemble();
        //   const auto& dofs = dbc.getDOFs();
        //   if (dbc.isComponent())
        //   {
        //     assert(false);
        //   }
        //   else
        //   {
        //     out.eliminate(dofs);
        //   }
        // }

        // if (trialFES == testFES)
        // {
        //   for (auto& pbc : pb.getPBCs())
        //   {
        //     pbc.assemble();
        //     const auto& dofs = pbc.getDOFs();

        //     if (pbc.isComponent())
        //     {
        //       assert(false);
        //     }
        //     else
        //     {
        //       out.merge(dofs);
        //     }
        //   }
        // }
        // else
        // {
        //   assert(false); // Not handled yet
        // }
      }

      Generic* copy() const noexcept override
      {
        return new Generic(*this);
      }

    private:
      // Default<Context::Local>::Type<VectorType, LinearFormType> m_lfa;
      // Default<Context::Local>::Type<OperatorType, BilinearFormType> m_bfa;
  };
}

#endif
