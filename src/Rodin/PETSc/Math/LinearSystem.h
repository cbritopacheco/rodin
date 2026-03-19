/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_LINEARSYSTEM_H
#define RODIN_PETSC_MATH_LINEARSYSTEM_H

/**
 * @file LinearSystem.h
 * @brief PETSc specialization of Rodin linear system containers.
 *
 * Provides the `LinearSystem<::Mat, ::Vec>` specialization that stores
 * the system matrix @f$ A @f$, right-hand side vector @f$ \mathbf{b} @f$,
 * and solution vector @f$ \mathbf{x} @f$ as PETSc handles.  Also defines
 * the @ref Rodin::Math::LinearSystem<::Mat,::Vec>::FieldSplits inner
 * class for block-preconditioner support.
 *
 * ## Usage
 *
 * This class is normally not instantiated directly; instead it is created
 * internally by @ref Rodin::Variational::Problem specializations.  Solvers
 * such as @ref Rodin::Solver::KSP receive a reference to the linear
 * system through the `solve(LinearSystem&)` interface.
 *
 * @see Rodin::PETSc::Math::LinearSystem (convenience alias),
 *      Rodin::Solver::KSP,
 *      Rodin::PETSc::Variational::Problem
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
   * @brief Specialization of @ref Rodin::Math::LinearSystem for PETSc
   *        objects (`::Mat` and `::Vec`).
   *
   * Couples a PETSc matrix (operator @f$ A @f$), right-hand side vector
   * (@f$ \mathbf{b} @f$), and solution vector (@f$ \mathbf{x} @f$) into
   * the common @ref Rodin::Math::LinearSystemBase interface used by Rodin
   * variational problems and solvers.
   *
   * ## Implementation Details
   *
   * The linear system @f$ A\mathbf{x} = \mathbf{b} @f$ is stored as
   * three PETSc handles that are allocated in the constructor and
   * destroyed in the destructor.  The `eliminate()` template method
   * enforces Dirichlet boundary conditions by zeroing rows of @f$ A @f$
   * and adjusting @f$ \mathbf{b} @f$ and @f$ \mathbf{x} @f$.
   *
   * ## Field-Split Support
   *
   * The inner @ref FieldSplits class stores named PETSc index sets
   * (`IS`) that can be passed to `PCFIELDSPLIT` for block
   * preconditioning of multi-field systems.
   *
   * @see Rodin::PETSc::Math::LinearSystem (convenience alias),
   *      Rodin::Solver::KSP
   */
  template <>
  class LinearSystem<::Mat, ::Vec>
    : public LinearSystemBase<::Mat, ::Vec, LinearSystem<::Mat, ::Vec>>
  {
    public:
      /// @brief PETSc matrix type (`::Mat`) for the system operator @f$ A @f$.
      using MatrixType =
        ::Mat;

      /// @brief PETSc vector type (`::Vec`) for @f$ \mathbf{b} @f$ and @f$ \mathbf{x} @f$.
      using VectorType =
        ::Vec;

      /// @brief Parent CRTP base class providing the generic
      ///        `LinearSystemBase` interface.
      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      /**
       * @brief Container for named PETSc index sets used with field-split
       *        preconditioners (e.g. `PCFIELDSPLIT`).
       *
       * Each @ref Split pairs a human-readable field name (e.g.
       * `"velocity"`, `"pressure"`) with a PETSc `IS` that identifies the
       * DOF indices belonging to that field within the global system.
       * The class takes ownership of the `IS` handles and destroys them
       * in its destructor.
       *
       * @see Rodin::PETSc::Variational::Problem::setFieldSplits
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

          /// @brief Default constructor.
          FieldSplits() = default;

          /**
           * @brief Constructs from a vector of named splits.
           * @param[in] splits Vector of Split objects.
           */
          explicit FieldSplits(std::vector<Split> splits)
            : m_splits(std::move(splits))
          {}

          /// @brief Copy constructor (deep-copies PETSc index sets).
          FieldSplits(const FieldSplits& other)
          {
            copyFrom(other);
          }

          /// @brief Move constructor.
          FieldSplits(FieldSplits&& other) noexcept
            : m_splits(std::move(other.m_splits))
          {}

          /// @brief Destructor; destroys owned PETSc index sets.
          ~FieldSplits()
          {
            clear();
          }

          /**
           * @brief Copy assignment operator.
           * @param[in] other FieldSplits to copy.
           * @return Reference to this FieldSplits.
           */
          FieldSplits& operator=(const FieldSplits& other)
          {
            if (this != &other)
            {
              clear();
              copyFrom(other);
            }
            return *this;
          }

          /**
           * @brief Move assignment operator.
           * @param[in] other FieldSplits to move from.
           * @return Reference to this FieldSplits.
           */
          FieldSplits& operator=(FieldSplits&& other) noexcept
          {
            if (this != &other)
            {
              clear();
              m_splits = std::move(other.m_splits);
            }
            return *this;
          }

          /// @brief Returns a mutable reference to the splits vector.
          std::vector<Split>& getSplits() noexcept
          {
            return m_splits;
          }

          /// @brief Returns a read-only reference to the splits vector.
          const std::vector<Split>& getSplits() const noexcept
          {
            return m_splits;
          }

          /// @brief Destroys all owned PETSc index sets and clears the splits.
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

          /**
           * @brief Deep-copies splits from another FieldSplits object.
           * @param[in] other Source FieldSplits to copy from.
           */
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

          /**
           * @brief Returns a mutable reference to the split at index @p i.
           * @param[in] i Split index.
           * @returns Reference to the split.
           */
          Split& operator[](size_t i)
          {
            return m_splits[i];
          }

          /**
           * @brief Returns a read-only reference to the split at index @p i.
           * @param[in] i Split index.
           * @returns Const reference to the split.
           */
          const Split& operator[](size_t i) const
          {
            return m_splits[i];
          }

          /// @brief Returns the number of splits.
          size_t size() const noexcept
          {
            return m_splits.size();
          }

          /// @brief Returns true if there are no splits.
          bool empty() const noexcept
          {
            return m_splits.empty();
          }

        private:
          std::vector<Split> m_splits; ///< The stored field splits.
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
       * @brief Enforces Dirichlet boundary conditions by eliminating
       *        constrained DOFs from the linear system.
       *
       * For every entry @f$ (i, u_i) @f$ in @p dofs the method:
       * 1. Sets @f$ x_i = u_i @f$ in the solution vector.
       * 2. Zeros row @f$ i @f$ of the matrix and places a @f$ 1 @f$ on
       *    the diagonal (`MatZeroRows` with `diag = 1.0`).
       * 3. Adjusts @f$ \mathbf{b} @f$ so that @f$ x_i = u_i @f$ is
       *    enforced after the solve.
       *
       * All indices in @p dofs are shifted by @p offset, which is useful
       * for block systems where each field starts at a non-zero position.
       *
       * @tparam DOFScalar Scalar type of the prescribed boundary values.
       * @param[in] dofs   Map from global DOF index to prescribed value
       *                   @f$ u_i @f$.
       * @param[in] offset Global index offset added to every key in
       *                   @p dofs (default 0).
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
