/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Generic.h
 * @brief Generic problem assembly implementation for complete variational problems.
 *
 * This file defines the Generic assembly class which handles the complete assembly
 * of variational problems including bilinear forms, linear forms, and boundary
 * conditions into a linear system ready for solution.
 */
#ifndef RODIN_ASSEMBLY_GENERIC_H
#define RODIN_ASSEMBLY_GENERIC_H

#include "Rodin/Variational/ForwardDecls.h"

#include "Rodin/FormLanguage/List.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Generic assembly implementation for complete variational problems.
   *
   * This class provides a complete assembly implementation that handles all
   * components of a variational problem:
   * - Bilinear forms @f$ a(u,v) @f$ assembled into stiffness matrix @f$ A @f$
   * - Linear forms @f$ l(v) @f$ assembled into load vector @f$ b @f$
   * - Dirichlet boundary conditions eliminating constrained DOFs
   * - Periodic boundary conditions merging dependent DOFs
   *
   * The resulting linear system has the form:
   * @f[
   *   A u = b
   * @f]
   * after applying boundary conditions.
   *
   * @tparam LinearSystem Linear system type (contains operator, solution, and RHS)
   * @tparam TrialFunction Trial function type
   * @tparam TestFunction Test function type
   *
   * ## Assembly Process
   * 1. **Bilinear Form Assembly**: Assembles all bilinear form integrators into
   *    the stiffness matrix @f$ A @f$
   * 2. **Linear Form Assembly**: Assembles all linear form integrators into the
   *    load vector @f$ b @f$ (stored as negative mass for compatibility)
   * 3. **Dirichlet BC Application**: Eliminates rows/columns corresponding to
   *    essential boundary conditions
   * 4. **Periodic BC Application**: Merges dependent degrees of freedom for
   *    periodic boundary conditions
   *
   * ## Usage Example
   * ```cpp
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(v) + DirichletBC(u, 0.0);
   * 
   * Generic<LinearSystem, decltype(u), decltype(v)> assembly;
   * LinearSystem ls;
   * assembly.execute(ls, problemInput);
   * ```
   */
  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Generic<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<LinearSystem, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      /// @brief Linear system type
      using LinearSystemType = LinearSystem;

      /// @brief Trial finite element space type
      using TrialFESType = typename FormLanguage::Traits<TrialFunction>::FESType;

      /// @brief Test finite element space type
      using TestFESType = typename FormLanguage::Traits<TestFunction>::FESType;

      /// @brief Trial space mesh type
      using TrialFESMeshType = typename FormLanguage::Traits<TrialFESType>::MeshType;

      /// @brief Test space mesh type
      using TestFESMeshType = typename FormLanguage::Traits<TestFESType>::MeshType;

      /// @brief Trial space mesh context type
      using TrialFESMeshContextType = typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      /// @brief Test space mesh context type
      using TestFESMeshContextType = typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      /// @brief Operator type (matrix)
      using OperatorType = typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      /// @brief Vector type
      using VectorType = typename FormLanguage::Traits<LinearSystemType>::VectorType;

      /// @brief Scalar type
      using ScalarType = typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      /// @brief Solution variable type
      using SolutionType = typename FormLanguage::Traits<TrialFunction>::SolutionType;

      /// @brief Bilinear form type
      using BilinearFormType =
        Variational::BilinearForm<SolutionType, TrialFESType, TestFESType, OperatorType>;

      /// @brief Default bilinear form assembly type
      using BilinearFormAssemblyType =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<OperatorType, BilinearFormType>;

      /// @brief Linear form type
      using LinearFormType = Variational::LinearForm<TestFESType, VectorType>;

      /// @brief Default linear form assembly type
      using LinearFormAssemblyType =
        typename Assembly::Default<TestFESMeshContextType>
          ::template Type<VectorType, LinearFormType>;

      /// @brief Parent class type
      using Parent =
        AssemblyBase<LinearSystemType, Variational::Problem<LinearSystem, TrialFunction, TestFunction>>;

      /// @brief Input data type
      using InputType = typename Parent::InputType;

      /// @brief Default constructor
      Generic() = default;

      /**
       * @brief Copy constructor.
       * @param other Generic assembly to copy
       */
      Generic(const Generic& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Generic assembly to move
       */
      Generic(Generic&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the complete problem assembly.
       *
       * Assembles the complete linear system by processing all bilinear forms,
       * linear forms, and boundary conditions specified in the problem.
       *
       * @param out Output linear system containing matrix, solution, and RHS
       * @param input Problem assembly input data
       *
       * ## Assembly Steps
       * 1. Assemble bilinear forms into stiffness matrix
       * 2. Add pre-assembled bilinear forms
       * 3. Assemble linear forms into load vector
       * 4. Add pre-assembled linear forms
       * 5. Apply Dirichlet boundary conditions (eliminate DOFs)
       * 6. Apply periodic boundary conditions (merge DOFs)
       */
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

      /**
       * @brief Creates a polymorphic copy of this assembly object.
       * @return Generic* Pointer to a new copy
       */
      Generic* copy() const noexcept override
      {
        return new Generic(*this);
      }
  };
}

#endif
