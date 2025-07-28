#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <type_traits>

#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

#include "Rodin/PETSc/Assembly/Generic.h"

namespace Rodin::Variational
{
  template <class U, class V>
  class Problem<PETSc::Math::LinearSystem, U, V> : public Variational::ProblemUVBase<PETSc::Math::LinearSystem, U, V>
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

      using Parent =
        Variational::ProblemUVBase<LinearSystemType, U, V>;

      Problem(U& u, V& v)
        : Parent(u, v)
      {
        using TrialFESContextType = typename FormLanguage::Traits<TrialFESType>::ContextType;
        if constexpr (std::is_same_v<TrialFESContextType, Context::Local>)
        {
          m_axb.create(PETSC_COMM_SELF);
        }
        else if constexpr (std::is_same_v<TrialFESContextType, Context::MPI>)
        {
          const auto& fes = u.getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          const auto& ctx = mesh.getContext();
          const MPI_Comm comm = ctx.getCommunicator();
          m_axb.create(comm);
        }
        else
        {
          assert(false);
        }
      }

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
        Parent::operator=(rhs);
        return *this;
      }

      Problem& assemble() override
      {
        // m_assembly.execute(m_axb, { m_pb, this->getTrialFunction(), this->getTestFunction() });
        m_assembled = true;
        return *this;
      }

      void solve(SolverBaseType& solver) override
      {
         // auto& axb = this->getLinearSystem();
         // if (!m_assembled)
         //    this->assemble();
         // solver.solve(axb);
         // this->getTrialFunction().getSolution().setData(axb.getSolution());
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

  template <class U, class V>
  Problem(U&, V&) -> Problem<U, V>;
}

namespace Rodin::PETSc::Variational
{
  template <class TrialFunction, class TestFunction>
  using Problem =
    Rodin::Variational::Problem<PETSc::Math::LinearSystem, TrialFunction, TestFunction>;
}

#endif
