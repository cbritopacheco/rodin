/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_HPP
#define RODIN_VARIATIONAL_PROBLEM_HPP

#include "Rodin/Utility.h"

#include "GridFunction.h"
#include "DirichletBC.h"

#include "Problem.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
   // ------------------------------------------------------------------------
   // ---- Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>
   // ------------------------------------------------------------------------

   template <class TrialFES, class TestFES>
   constexpr
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>
   ::Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
      :  m_bilinearForm(u, v),
         m_linearForm(v),
         m_trialFunction(u),
         m_testFunction(v)
   {}

   template <class TrialFES, class TestFES>
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>&
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>::operator=(ProblemBody&& rhs)
   {
      Parent::operator=(std::move(rhs));

      for (auto& bfi : getProblemBody().getBFIs())
         getBilinearForm().add(bfi);

      for (auto& lfi : getProblemBody().getLFIs())
         getLinearForm().add(UnaryMinus(lfi)); // Negate every linear form

      return *this;
   }

   template <class TrialFES, class TestFES>
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>&
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>::update()
   {
      // Update all components of the problem
      getTrialFunction().getFiniteElementSpace().update();
      getTestFunction().getFiniteElementSpace().update();
      getTrialFunction().getGridFunction().update();
      getLinearForm().update();
      getBilinearForm().update();

      return *this;
   }

   template <class TrialFES, class TestFES>
   void
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>::assemble()
   {
      // Assemble both sides
      getLinearForm().assemble();
      getBilinearForm().assemble();

      // Project values onto the essential boundary and compute essential dofs
      getEssentialTrueDOFs().SetSize(0);

      const auto& uuid = getTrialFunction().getUUID();
      auto tfIt = getEssentialBoundary().getTFMap().find(uuid);
      if (tfIt != getEssentialBoundary().getTFMap().end())
      {
         const auto& tfValue = tfIt->second;
         const auto& bdrAttr = tfValue.attributes;
         getTrialFunction().getGridFunction().projectOnBoundary(*tfValue.value, bdrAttr);
         getEssentialTrueDOFs().Append(getTrialFunction().getFiniteElementSpace().getEssentialTrueDOFs(bdrAttr));
      }

      auto tfCompIt = getEssentialBoundary().getTFCompMap().find(uuid);
      if (tfCompIt != getEssentialBoundary().getTFCompMap().end())
      {
         const auto& compMap = tfCompIt->second;
         for (const auto& [component, compValue] : compMap)
         {
            const auto& bdrAttr = compValue.attributes;
            Component comp(getTrialFunction().getGridFunction(), component);
            comp.projectOnBoundary(*compValue.value, bdrAttr);
            getEssentialTrueDOFs().Append(
               getTrialFunction().getFiniteElementSpace().getEssentialTrueDOFs(bdrAttr, component));
         }
      }

      getEssentialTrueDOFs().Sort();
      getEssentialTrueDOFs().Unique();

      // Form linear system
      if constexpr (std::is_same_v<TrialFES, TestFES>)
      {
         mfem::Operator* ptr = nullptr;
         getBilinearForm().getOperator()
            .FormLinearSystem(
               getEssentialTrueDOFs(),
               getTrialFunction().getGridFunction().getData(),
               getLinearForm().getVector(),
               ptr,
               m_guess,
               m_massVector);
         assert(m_stiffnessOp);
         m_stiffnessOp.reset(ptr);
         ptr = nullptr;
      }
      else
      {
         assert(false); // Not supported yet
      }
   }

   template <class TrialFES, class TestFES>
   void
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>
   ::solve(const Solver::SolverBase<OperatorType, VectorType>& solver)
   {
      // Emplace data
      getTrialFunction().emplaceGridFunction();

      // Assemble the system
      assemble();

      // Solve the system
      solver.solve(getStiffnessOperator(), getMassVector(), m_guess);

      // Recover solution
      getBilinearForm().getOperator()
         .RecoverFEMSolution(
            m_guess,
            getLinearForm().getVector(),
            getTrialFunction().getGridFunction().getData());

   }

   template <class TrialFES, class TestFES>
   const EssentialBoundary&
   Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>
   ::getEssentialBoundary() const
   {
      return getProblemBody().getEssentialBoundary();
   }

   // ------------------------------------------------------------------------
   // ---- Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
   // ------------------------------------------------------------------------

   template <class TrialFES, class TestFES>
   constexpr
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
   ::Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
      :  m_bilinearForm(u, v),
         m_linearForm(v),
         m_trialFunction(u),
         m_testFunction(v)
   {}

   template <class TrialFES, class TestFES>
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>&
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
   ::operator=(ProblemBody&& rhs)
   {
      Parent::operator=(std::move(rhs));

      for (auto& bfi : getProblemBody().getBFIs())
         getBilinearForm().add(bfi);

      for (auto& lfi : getProblemBody().getLFIs())
         getLinearForm().add(UnaryMinus(lfi)); // Negate every linear form

      return *this;
   }

   template <class TrialFES, class TestFES>
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>&
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>::update()
   {
      // Update all components of the problem
      getTrialFunction().getFiniteElementSpace().update();
      getTestFunction().getFiniteElementSpace().update();
      getTrialFunction().getGridFunction().update();
      getLinearForm().update();
      getBilinearForm().update();

      return *this;
   }

   template <class TrialFES, class TestFES>
   void
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>::assemble()
   {
      // Assemble both sides
      getLinearForm().assemble();
      getBilinearForm().assemble();

      // Project values onto the essential boundary and compute essential dofs
      getEssentialTrueDOFs().SetSize(0);

      const auto& uuid = getTrialFunction().getUUID();
      auto tfIt = getEssentialBoundary().getTFMap().find(uuid);
      if (tfIt != getEssentialBoundary().getTFMap().end())
      {
         const auto& tfValue = tfIt->second;
         const auto& bdrAttr = tfValue.attributes;
         getTrialFunction().getGridFunction().projectOnBoundary(*tfValue.value, bdrAttr);
         getEssentialTrueDOFs().Append(getTrialFunction().getFiniteElementSpace().getEssentialTrueDOFs(bdrAttr));
      }

      auto tfCompIt = getEssentialBoundary().getTFCompMap().find(uuid);
      if (tfCompIt != getEssentialBoundary().getTFCompMap().end())
      {
         const auto& compMap = tfCompIt->second;
         for (const auto& [component, compValue] : compMap)
         {
            const auto& bdrAttr = compValue.attributes;
            Component comp(getTrialFunction().getGridFunction(), component);
            comp.projectOnBoundary(*compValue.value, bdrAttr);
            getEssentialTrueDOFs().Append(
               getTrialFunction().getFiniteElementSpace().getEssentialTrueDOFs(bdrAttr, component));
         }
      }

      getEssentialTrueDOFs().Sort();
      getEssentialTrueDOFs().Unique();

      // Form linear system
      getBilinearForm().formLinearSystem(
            getEssentialTrueDOFs(),
            getTrialFunction().getGridFunction().getData(),
            getLinearForm().getVector(),
            m_stiffnessOp,
            m_guess,
            m_massVector);
   }

   template <class TrialFES, class TestFES>
   void
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
   ::solve(const Solver::SolverBase<OperatorType, VectorType>& solver)
   {
      // Emplace data
      getTrialFunction().emplaceGridFunction();

      // Assemble the system
      assemble();

      // Solve the system
      solver.solve(getStiffnessOperator(), getMassVector(), m_guess);

      // Recover solution
      getBilinearForm().getOperator()
         .RecoverFEMSolution(
            m_guess,
            getLinearForm().getVector(),
            getTrialFunction().getGridFunction().getData());

   }

   template <class TrialFES, class TestFES>
   const EssentialBoundary&
   Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
   ::getEssentialBoundary() const
   {
      return getProblemBody().getEssentialBoundary();
   }
}

#endif
