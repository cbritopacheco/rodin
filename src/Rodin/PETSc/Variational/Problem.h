#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <type_traits>

#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

namespace Rodin::PETSc
{
  template <class U, class V>
  class Problem : public Variational::ProblemUVBase<LinearSystem, U, V>
  {
    public:
      using Parent = Variational::ProblemUVBase<LinearSystem, U, V>;

      using LinearSystemType = LinearSystem;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      using ProblemBodyType = Variational::ProblemBody<OperatorType, VectorType, ScalarType>;

      using TrialFESType = typename FormLanguage::Traits<U>::FESType;

      using TestFESType = typename FormLanguage::Traits<V>::FESType;

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
      LinearSystemType m_axb;
  };

  template <class U, class V>
  Problem(U&, V&) -> Problem<U, V>;
}

#endif
