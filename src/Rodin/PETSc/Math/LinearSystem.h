/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_LINEARSYSTEM_H
#define RODIN_PETSC_MATH_LINEARSYSTEM_H

/**
 * @file
 * @brief PETSc specialization of Rodin linear system containers.
 */

#include <petsc.h>
#include <boost/mpi/communicator.hpp>
#include <petscmat.h>
#include <petscvec.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"

namespace Rodin::Math
{
  /**
   * @brief Specialization of @ref Rodin::Math::LinearSystem for PETSc objects.
   *
   * Couples a PETSc operator (@f$ A @f$), right-hand side vector (@f$ b @f$),
   * and solution vector (@f$ x @f$) into the common
   * @ref Rodin::Math::LinearSystemBase interface used by Rodin variational
   * problems and solvers.
   *
   * ## Implementation Details
   *
   * The linear system @f$ Ax = b @f$ is stored as three PETSc handles
   * owned by this object and destroyed in the destructor.
   *
   * @see Rodin::PETSc::Math::LinearSystem, Rodin::Solver::KSP
   */
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

      /**
       * @brief Container for named PETSc index sets used with field-split
       * preconditioners (e.g. `PCFIELDSPLIT`).
       *
       * Each @ref Split pairs a human-readable name with a PETSc `IS` that
       * identifies the DOF indices belonging to that field.
       */
      class FieldSplits
      {
        public:
          /**
           * @brief A single named index set for field-split preconditioning.
           */
          struct Split
          {
            std::string name;        ///< Human-readable field name.
            ::IS is = PETSC_NULLPTR; ///< PETSc index set for this field's DOFs.
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

      /**
       * @brief Constructs an empty linear system on the given communicator.
       * @param comm MPI communicator for distributed PETSc objects.
       */
      LinearSystem(MPI_Comm comm);

      /// @brief Copy constructor (deep-copies PETSc handles).
      LinearSystem(const LinearSystem& other);

      /// @brief Move constructor.
      LinearSystem(LinearSystem&& other) noexcept;

      /// @brief Destructor; destroys owned PETSc matrix and vector handles.
      virtual ~LinearSystem();

      /// @brief Copy-assignment operator.
      LinearSystem& operator=(const LinearSystem& other);

      /// @brief Move-assignment operator.
      LinearSystem& operator=(LinearSystem&& other) noexcept;

      /**
       * @brief Eliminates Dirichlet DOFs from the linear system.
       *
       * Zeros the matrix rows corresponding to constrained DOFs and adjusts
       * the right-hand side vector so that the solution satisfies
       * @f$ x_i = u_i @f$ for every constrained index @f$ i @f$.
       *
       * @tparam DOFScalar Scalar type of the prescribed values.
       * @param  dofs   Map from global DOF index to prescribed value.
       * @param  offset Global index offset added to every key in @p dofs.
       * @returns Reference to `*this` for method chaining.
       */
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

      /**
       * @brief Merges periodic DOF constraints into the linear system.
       * @note Currently unimplemented — throws at runtime.
       */
      template <class DOFScalar>
      LinearSystemBase& merge(
          const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs, size_t offset = 0)
      {
        throw "Unimplemented.";
        return *this;
      }

      /// @brief Returns the MPI communicator of this linear system.
      constexpr
      MPI_Comm getCommunicator() const noexcept
      {
        return m_comm;
      }

      /// @brief Returns a mutable reference to the system matrix @f$ A @f$.
      constexpr
      MatrixType& getOperator()
      {
        return m_operator;
      }

      /// @brief Returns a read-only reference to the system matrix @f$ A @f$.
      constexpr
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      /// @brief Returns a mutable reference to the right-hand side vector @f$ b @f$.
      constexpr
      VectorType& getVector()
      {
        return m_vector;
      }

      /// @brief Returns a read-only reference to the right-hand side vector @f$ b @f$.
      constexpr
      const VectorType& getVector() const
      {
        return m_vector;
      }

      /// @brief Returns a mutable reference to the solution vector @f$ x @f$.
      constexpr
      VectorType& getSolution()
      {
        return m_solution;
      }

      /// @brief Returns a read-only reference to the solution vector @f$ x @f$.
      constexpr
      const VectorType& getSolution() const
      {
        return m_solution;
      }

      /**
       * @brief Sets the field-split index sets for block preconditioners.
       * @param fields Field splits to associate with this system.
       * @returns Reference to `*this`.
       */
      LinearSystem& setFieldSplits(FieldSplits fields)
      {
        m_fieldSplits = std::move(fields);
        return *this;
      }

      /// @brief Returns the current field-split configuration.
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
  /**
   * @brief Convenient alias for the PETSc linear system specialization.
   *
   * @see Rodin::Math::LinearSystem
   */
  using LinearSystem = Rodin::Math::LinearSystem<Matrix, Vector>;
}

#endif
