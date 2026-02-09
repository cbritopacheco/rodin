/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_LINEARSYSTEM_H
#define RODIN_PETSC_MATH_LINEARSYSTEM_H

#include <petsc.h>
#include <boost/mpi/communicator.hpp>
#include <petscmat.h>
#include <petscvec.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"

namespace Rodin::Math
{
  template <>
  class LinearSystem<::Mat, ::Vec>
    : public LinearSystemBase<::Mat, ::Vec, LinearSystem<::Mat, ::Vec>>
  {
    public:
      using MatrixType =
        ::Mat;

      using VectorType =
        ::Vec;

      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      class FieldSplits
      {
        public:
          struct Split
          {
            std::string name;
            ::IS is = PETSC_NULLPTR;
          };

          FieldSplits() = default;

          explicit FieldSplits(std::vector<Split> splits)
            : m_splits(std::move(splits))
          {}

          FieldSplits(const FieldSplits& other)
          {
            copyFrom(other);
          }

          FieldSplits(FieldSplits&& other) noexcept
            : m_splits(std::move(other.m_splits))
          {}

          ~FieldSplits()
          {
            clear();
          }

          FieldSplits& operator=(const FieldSplits& other)
          {
            if (this != &other)
            {
              clear();
              copyFrom(other);
            }
            return *this;
          }

          FieldSplits& operator=(FieldSplits&& other) noexcept
          {
            if (this != &other)
            {
              clear();
              m_splits = std::move(other.m_splits);
            }
            return *this;
          }

          std::vector<Split>& getSplits() noexcept
          {
            return m_splits;
          }

          const std::vector<Split>& getSplits() const noexcept
          {
            return m_splits;
          }

          void clear() noexcept
          {
            PetscErrorCode ierr;
            for (auto& s : m_splits)
            {
              if (s.is)
              {
                ierr = ISDestroy(&s.is);
                assert(ierr == PETSC_SUCCESS);
                s.is = PETSC_NULLPTR;
              }
            }
            m_splits.clear();
            (void) ierr;
          }

          void copyFrom(const FieldSplits& other)
          {
            PetscErrorCode ierr;

            m_splits.reserve(other.m_splits.size());
            for (const auto& s : other.m_splits)
            {
              Split out;
              out.name = s.name;

              if (s.is)
              {
                ierr = ISDuplicate(s.is, &out.is);
                assert(ierr == PETSC_SUCCESS);
                ierr = ISCopy(s.is, out.is);
                assert(ierr == PETSC_SUCCESS);
              }
              else
              {
                out.is = PETSC_NULLPTR;
              }

              m_splits.push_back(std::move(out));
            }

            (void) ierr;
          }

          Split& operator[](size_t i)
          {
            return m_splits[i];
          }

          const Split& operator[](size_t i) const
          {
            return m_splits[i];
          }

          size_t size() const noexcept
          {
            return m_splits.size();
          }

          bool empty() const noexcept
          {
            return m_splits.empty();
          }

        private:
          std::vector<Split> m_splits;
      };

      LinearSystem(MPI_Comm comm);

      LinearSystem(const LinearSystem& other);

      LinearSystem(LinearSystem&& other) noexcept;

      virtual ~LinearSystem();

      LinearSystem& operator=(const LinearSystem& other);

      LinearSystem& operator=(LinearSystem&& other) noexcept;

      template <class DOFScalar>
      LinearSystem& eliminate(const IndexMap<DOFScalar>& dofs, size_t offset = 0)
      {
        auto& a = this->getOperator();
        auto& b = this->getVector();
        auto& x = this->getSolution();

        PetscErrorCode ierr;

        std::vector<PetscInt> rows;
        rows.reserve(dofs.size());
        for (auto const& kv : dofs)
          rows.push_back(kv.first + offset);

        for (auto const& kv : dofs)
        {
          const PetscInt i = kv.first + offset;
          const auto& ui = kv.second;
          ierr = VecSetValue(x, i, ui, INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);
        }

        ierr = VecAssemblyBegin(x);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyEnd(x);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroRows(a, rows.size(), rows.data(), 1.0, x, b);
        assert(ierr == PETSC_SUCCESS);

        (void) ierr;

        return *this;
      }

      template <class DOFScalar>
      LinearSystemBase& merge(
          const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        throw "Unimplemented.";
        return *this;
      }

      constexpr
      MPI_Comm getCommunicator() const noexcept
      {
        return m_comm;
      }

      constexpr
      MatrixType& getOperator()
      {
        return m_operator;
      }

      constexpr
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      constexpr
      VectorType& getVector()
      {
        return m_vector;
      }

      constexpr
      const VectorType& getVector() const
      {
        return m_vector;
      }

      constexpr
      VectorType& getSolution()
      {
        return m_solution;
      }

      constexpr
      const VectorType& getSolution() const
      {
        return m_solution;
      }

      LinearSystem& setFieldSplits(FieldSplits fields)
      {
        m_fieldSplits = std::move(fields);
        return *this;
      }

      const FieldSplits& getFieldSplits() const noexcept
      {
        return m_fieldSplits;
      }

    private:
      MPI_Comm m_comm; ///< The MPI communicator for the linear system.
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
      VectorType m_vector;   ///< The vector of the linear system.

      FieldSplits m_fieldSplits; ///< The field splits for block preconditioning.
  };
}

namespace Rodin::PETSc::Math
{
  using LinearSystem = Rodin::Math::LinearSystem<Matrix, Vector>;
}

#endif
