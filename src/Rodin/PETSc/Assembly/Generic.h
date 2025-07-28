/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_GENERIC_H
#define RODIN_PETSC_ASSEMBLY_GENERIC_H

#include "Rodin/Assembly/ForwardDecls.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

namespace Rodin::Assembly
{
  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Generic<LinearSystem, Variational::Problem<PETSc::Math::LinearSystem, TrialFunction, TestFunction>> final
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
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using Generic = Rodin::Assembly::Generic<LinearAlgebraType, Operand>;
}

#endif
