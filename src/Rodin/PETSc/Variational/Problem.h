#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <type_traits>

#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/PETSc/Variational/TrialFunction.h"
#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

#include "Rodin/PETSc/Assembly/Generic.h"
#include "Rodin/Variational/TestFunction.h"

namespace Rodin::Variational
{
  template <class U, class V>
  class Problem<PETSc::Math::LinearSystem, U, V>
    : public Variational::ProblemUVBase<PETSc::Math::LinearSystem, U, V>
  {
    public:
      using LinearSystemType =
        PETSc::Math::LinearSystem;

      using AssemblyType =
        PETSc::Assembly::Generic<LinearSystemType, Problem>;

      using SolverBaseType =
        Solver::SolverBase<LinearSystemType>;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using ProblemBodyType = Variational::ProblemBody<OperatorType, VectorType, ScalarType>;

      using TrialFESType =
        typename FormLanguage::Traits<U>::FESType;

      using TestFESType =
        typename FormLanguage::Traits<V>::FESType;

      using TrialFESMeshType =
        typename FormLanguage::Traits<TrialFESType>::MeshType;

      using TrialFESMeshContextType =
        typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      using TestFESMeshType =
        typename FormLanguage::Traits<TestFESType>::MeshType;

      using TestFESMeshContextType =
        typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      using Parent =
        Variational::ProblemUVBase<LinearSystemType, U, V>;

      static_assert(
          std::is_same_v<TrialFESMeshContextType, Context::Local> ||
          std::is_same_v<TrialFESMeshContextType, Context::MPI>);

      static_assert(
          std::is_same_v<TrialFESMeshContextType, TestFESMeshContextType>);

      Problem(U& u, V& v)
        : Parent(u, v),
          m_axb(
              [&]() -> MPI_Comm
              {
                if constexpr (std::is_same_v<TrialFESMeshContextType, Context::Local>)
                {
                  return PETSC_COMM_SELF;
                }
                else if constexpr (std::is_same_v<TrialFESMeshContextType, Context::MPI>)
                {
                  const auto& fes = u.getFiniteElementSpace();
                  const auto& mesh = fes.getMesh();
                  const auto& ctx = mesh.getContext();
                  const MPI_Comm comm = ctx.getCommunicator();
                  return comm;
                }
                else
                {
                  assert(false);
                }
              }())
      {}

      constexpr
      Problem(const Problem& other)
        : Parent(other),
          m_axb(other.m_axb)
      {}

      constexpr
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_axb(std::move(other.m_axb))
      {}

      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_axb = other.m_axb;
        }
        return *this;
      }

      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_axb = std::move(other.m_axb);
        }
        return *this;
      }

      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      Problem& assemble() override
      {
        m_assembly.execute(m_axb, { m_pb, this->getTrialFunction(), this->getTestFunction() });
        m_assembled = true;
        return *this;
      }

      void solve(SolverBaseType& solver) override
      {
        auto& axb = this->getLinearSystem();
        if (!m_assembled)
           this->assemble();
        solver.solve(axb);
        this->getTrialFunction().getSolution().setData(axb.getSolution());
      }

      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

    private:
      Boolean m_assembled;
      ProblemBodyType m_pb;
      LinearSystemType m_axb;
      AssemblyType m_assembly;
  };

  template <class Solution, class TrialFES, class TestFES>
  Problem(
      PETSc::Variational::TrialFunction<Solution, TrialFES>&,
      PETSc::Variational::TestFunction<TestFES>&)
    -> Problem<
          PETSc::Math::LinearSystem,
          TrialFunction<Solution, TrialFES>,
          TestFunction<TestFES>>;
}

namespace Rodin::PETSc::Variational
{
  template <class TrialFunction, class TestFunction>
  using Problem =
    Rodin::Variational::Problem<PETSc::Math::LinearSystem, TrialFunction, TestFunction>;
}

#endif
