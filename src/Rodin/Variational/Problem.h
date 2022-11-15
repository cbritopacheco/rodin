/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <set>
#include <variant>
#include <functional>

#include <mfem.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/Solver/Solver.h"

#include "ForwardDecls.h"

#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "EssentialBoundary.h"


namespace Rodin::Variational
{
   template <class TrialFES, class TestFES>
   Problem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
      -> Problem<TrialFES, TestFES, typename TrialFES::Context, mfem::SparseMatrix, mfem::Vector>;

   template <class OperatorType, class VectorType>
   class ProblemBase : public FormLanguage::Base
   {
      public:
         ProblemBase() = default;

         ProblemBase(ProblemBase&& other) = default;

         ProblemBase(const ProblemBase& other) = default;

         const ProblemBody& getProblemBody() const
         {
            return m_pb;
         }

         virtual ProblemBase& operator=(ProblemBody&& rhs)
         {
            m_pb = std::move(rhs);
            return *this;
         }

         virtual bool isParallel() const = 0;

         virtual void solve(const Solver::SolverBase<OperatorType, VectorType>& solver) = 0;

         /**
          * @brief Assembles the underlying linear system to solve.
          */
         virtual void assemble() = 0;

         /**
          * @brief Updates the ProblemBase instance after a refinement in the mesh
          */
         virtual ProblemBase& update() = 0;

         /**
          * @returns Reference to the mfem::Operator representing the stiffness
          * matrix.
          *
          * This must be called only after assemble() has been called.
          */
         virtual OperatorType& getStiffnessOperator() = 0;

         /**
          * @returns Constant reference to the mfem::Operator representing the stiffness
          * matrix, i.e. the LHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual const OperatorType& getStiffnessOperator() const = 0;

         /**
          * @returns Reference to the mfem::Vector representing the mass
          * vector, i.e. the RHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual VectorType& getMassVector() = 0;

         /**
          * @returns Constant reference to the mfem::Vector representing the mass
          * vector, i.e. the RHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual const VectorType& getMassVector() const = 0;

         virtual const EssentialBoundary& getEssentialBoundary() const = 0;

         virtual ProblemBase* copy() const noexcept override = 0;

      private:
         ProblemBody m_pb;
   };

   /**
    * @ingroup ProblemSpecializations
    */
   template <class TrialFES, class TestFES>
   class Problem<TrialFES, TestFES, Context::Serial, mfem::Operator, mfem::Vector>
      : public ProblemBase<mfem::Operator, mfem::Vector>
   {
         static_assert(std::is_same_v<typename TrialFES::Context, Context::Serial>);
         static_assert(std::is_same_v<typename TestFES::Context, Context::Serial>);

         using OperatorHandle = std::unique_ptr<mfem::Operator>;

      public:
         using Context = Context::Serial;
         using OperatorType = mfem::Operator;
         using VectorType = mfem::Vector;
         using Parent = ProblemBase<mfem::Operator, mfem::Vector>;

         /**
          * @brief Constructs an empty problem involving the trial function @f$ u @f$
          * and the test function @f$ v @f$.
          *
          * @param[in,out] u Trial function @f$ u @f$
          * @param[in,out] v Test function @f$ v @f$
          */
         explicit
         constexpr
         Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v);

         /**
          * @brief Deleted copy constructor.
          */
         Problem(const Problem& other) = delete;

         /**
          * @brief Deleted copy assignment operator.
          */
         void operator=(const Problem& other) = delete;

         constexpr
         TrialFunction<TrialFES>& getTrialFunction()
         {
            return m_trialFunction;
         }

         constexpr
         TestFunction<TestFES>& getTestFunction()
         {
            return m_testFunction;
         }

         constexpr
         const TrialFunction<TrialFES>& getTrialFunction() const
         {
            return m_trialFunction;
         }

         constexpr
         const TestFunction<TestFES>& getTestFunction() const
         {
            return m_testFunction;
         }

         constexpr
         LinearForm<TestFES, Context, VectorType>& getLinearForm()
         {
            return m_linearForm;
         }

         constexpr
         BilinearForm<TrialFES, TestFES, Context, OperatorType>& getBilinearForm()
         {
            return m_bilinearForm;
         }

         constexpr
         mfem::Array<int>& getEssentialTrueDOFs()
         {
            return m_essTrueDofList;
         }

         constexpr
         const mfem::Array<int>& getEssentialTrueDOFs() const
         {
            return m_essTrueDofList;
         }

         constexpr
         const LinearForm<TestFES, Context, VectorType>& getLinearForm() const
         {
            return m_linearForm;
         }

         constexpr
         const BilinearForm<TrialFES, TestFES, Context, OperatorType>& getBilinearForm() const
         {
            return m_bilinearForm;
         }

         bool isParallel() const override
         {
            return false;
         }

         void assemble() override;

         Problem& update() override;

         void solve(const Solver::SolverBase<OperatorType, VectorType>& solver) override;

         Problem& operator=(ProblemBody&& rhs) override;

         virtual VectorType& getMassVector() override
         {
            return m_massVector;
         }

         virtual const VectorType& getMassVector() const override
         {
            return m_massVector;
         }

         virtual OperatorType& getStiffnessOperator() override
         {
            return *m_stiffnessOp;
         }

         virtual const OperatorType& getStiffnessOperator() const override
         {
            return *m_stiffnessOp;
         }

         virtual const EssentialBoundary& getEssentialBoundary() const override;

         virtual Problem* copy() const noexcept override
         {
            assert(false);
            return nullptr;
         }

      private:
         TrialFunction<TrialFES>& m_trialFunction;
         TestFunction<TestFES>&   m_testFunction;

         LinearForm<TestFES, Context, VectorType> m_linearForm;
         BilinearForm<TrialFES, TestFES, Context, OperatorType> m_bilinearForm;

         OperatorHandle       m_stiffnessOp;
         mfem::Vector         m_massVector;
         mfem::Vector         m_guess;

         mfem::Array<int>     m_essTrueDofList;
   };

   /**
    * @ingroup ProblemSpecializations
    */
   template <class TrialFES, class TestFES>
   class Problem<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix, mfem::Vector>
      : public ProblemBase<mfem::SparseMatrix, mfem::Vector>
   {
         static_assert(std::is_same_v<typename TrialFES::Context, Context::Serial>);
         static_assert(std::is_same_v<typename TestFES::Context, Context::Serial>);

      public:
         using Context = Context::Serial;
         using OperatorType = mfem::SparseMatrix;
         using VectorType = mfem::Vector;
         using Parent = ProblemBase<mfem::SparseMatrix, mfem::Vector>;

         /**
          * @brief Constructs an empty problem involving the trial function @f$ u @f$
          * and the test function @f$ v @f$.
          *
          * @param[in,out] u Trial function @f$ u @f$
          * @param[in,out] v Test function @f$ v @f$
          */
         explicit
         constexpr
         Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v);

         /**
          * @brief Deleted copy constructor.
          */
         Problem(const Problem& other) = delete;

         /**
          * @brief Deleted copy assignment operator.
          */
         void operator=(const Problem& other) = delete;

         constexpr
         TrialFunction<TrialFES>& getTrialFunction()
         {
            return m_trialFunction;
         }

         constexpr
         TestFunction<TestFES>& getTestFunction()
         {
            return m_testFunction;
         }

         constexpr
         const TrialFunction<TrialFES>& getTrialFunction() const
         {
            return m_trialFunction;
         }

         constexpr
         const TestFunction<TestFES>& getTestFunction() const
         {
            return m_testFunction;
         }

         constexpr
         LinearForm<TestFES, Context, VectorType>& getLinearForm()
         {
            return m_linearForm;
         }

         constexpr
         BilinearForm<TrialFES, TestFES, Context, OperatorType>& getBilinearForm()
         {
            return m_bilinearForm;
         }

         constexpr
         mfem::Array<int>& getEssentialTrueDOFs()
         {
            return m_essTrueDofList;
         }

         constexpr
         const mfem::Array<int>& getEssentialTrueDOFs() const
         {
            return m_essTrueDofList;
         }

         constexpr
         const LinearForm<TestFES, Context, VectorType>& getLinearForm() const
         {
            return m_linearForm;
         }

         constexpr
         const BilinearForm<TrialFES, TestFES, Context, OperatorType>& getBilinearForm() const
         {
            return m_bilinearForm;
         }

         bool isParallel() const override
         {
            return false;
         }

         void assemble() override;

         Problem& update() override;

         void solve(const Solver::SolverBase<OperatorType, VectorType>& solver) override;

         Problem& operator=(ProblemBody&& rhs) override;

         virtual VectorType& getMassVector() override
         {
            return m_massVector;
         }

         virtual const VectorType& getMassVector() const override
         {
            return m_massVector;
         }

         virtual OperatorType& getStiffnessOperator() override
         {
            return m_stiffnessOp;
         }

         virtual const OperatorType& getStiffnessOperator() const override
         {
            return m_stiffnessOp;
         }

         virtual const EssentialBoundary& getEssentialBoundary() const override;

         virtual Problem* copy() const noexcept override
         {
            assert(false);
            return nullptr;
         }

      private:
         TrialFunction<TrialFES>& m_trialFunction;
         TestFunction<TestFES>&   m_testFunction;

         LinearForm<TestFES, Context, VectorType> m_linearForm;
         BilinearForm<TrialFES, TestFES, Context, OperatorType> m_bilinearForm;

         mfem::SparseMatrix   m_stiffnessOp;
         mfem::Vector         m_massVector;
         mfem::Vector         m_guess;

         mfem::Array<int>     m_essTrueDofList;
   };
}

#include "Problem.hpp"

#endif
