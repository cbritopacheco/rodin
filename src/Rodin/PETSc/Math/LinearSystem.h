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
 * Provides the PETSc-backed linear system specialization that stores
 * the system matrix @f$ A @f$, right-hand side vector @f$ \mathbf{b} @f$,
 * and solution vector @f$ \mathbf{x} @f$ as PETSc handles.  Also defines
 * the @c FieldSplits inner
 * class for block-preconditioner support.
 *
 * # Usage
 *
 * This class is normally not instantiated directly; instead it is created
 * internally by @ref Rodin::Variational::Problem specializations.  Solvers
 * such as @ref Rodin::Solver::KSP receive a reference to the linear
 * system through the `solve(LinearSystem&)` interface.
 *
 * @see <a href="class_rodin_1_1_math_1_1_linear_system_3_1_1_mat_00_01_1_1_vec_01_4.html">Rodin::PETSc::Math::LinearSystem</a> convenience alias
 * @see <a href="class_rodin_1_1_solver_1_1_k_s_p.html">Rodin::Solver::KSP</a>
 * @see <a href="class_rodin_1_1_variational_1_1_problem_3_01_p_e_t_sc_1_1_math_1_1_linear_system_00_01_u_00_01_v_01_4.html">Rodin::PETSc::Variational::Problem</a>
 */

#include <petsc.h>
#include <boost/mpi/communicator.hpp>
#include <petscmat.h>
#include <petscvec.h>
#include <string>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/Vector.h"

namespace Rodin::Math
{
  /**
   * @brief Specialization of @ref Rodin::Math::LinearSystem for PETSc
   *        objects (@c Mat and @c Vec).
   *
   * Couples a PETSc matrix (operator @f$ A @f$), right-hand side vector
   * (@f$ \mathbf{b} @f$), and solution vector (@f$ \mathbf{x} @f$) into
   * the common @ref Rodin::Math::LinearSystemBase interface used by Rodin
   * variational problems and solvers.
   *
   * # Implementation Details
   *
   * The linear system @f$ A\mathbf{x} = \mathbf{b} @f$ is stored as
   * three PETSc handles that are allocated in the constructor and
   * destroyed in the destructor.  The `eliminate()` template method
   * enforces Dirichlet boundary conditions by zeroing rows of @f$ A @f$
   * and adjusting @f$ \mathbf{b} @f$ and @f$ \mathbf{x} @f$.
   *
   * # Field-Split Support
   *
   * The inner @ref FieldSplits class stores named PETSc index sets
   * (`IS`) that can be passed to `PCFIELDSPLIT` for block
   * preconditioning of multi-field systems.
   *
   * @see <a href="class_rodin_1_1_math_1_1_linear_system_3_1_1_mat_00_01_1_1_vec_01_4.html">Rodin::PETSc::Math::LinearSystem</a> convenience alias
   * @see <a href="class_rodin_1_1_solver_1_1_k_s_p.html">Rodin::Solver::KSP</a>
   */
  template <>
  class LinearSystem<::Mat, ::Vec>
    : public LinearSystemBase<::Mat, ::Vec, LinearSystem<::Mat, ::Vec>>
  {
    public:
      /// @brief PETSc matrix type (@c Mat) for the system operator @f$ A @f$.
      using MatrixType = ::Mat;

      /// @brief PETSc vector type (@c Vec) for @f$ \mathbf{b} @f$ and @f$ \mathbf{x} @f$.
      using VectorType = ::Vec;

      /**
       * @brief Parent CRTP base class providing the generic
       *        `LinearSystemBase` interface.
       */
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
       * @see <a href="class_rodin_1_1_variational_1_1_problem_3_01_p_e_t_sc_1_1_math_1_1_linear_system_00_01_u1_00_01_u2_00_01_u3_00_01_us_8_8_8_01_4.html">Rodin::PETSc::Variational::Problem::setFieldSplits</a>
       */
      class FieldSplits
      {
        public:
          /// @brief A single named index set for field-split preconditioning.
          struct Split
          {
            std::string name; ///< Human-readable field name.
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

          /**
           * @brief Copy constructor (deep-copies PETSc index sets).
           * @param other Object to copy from.
           */
          FieldSplits(const FieldSplits& other)
          {
            copyFrom(other);
          }

          /**
           * @brief Move constructor.
           * @param other Object to move from.
           */
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

          /**
           * @brief Returns a mutable reference to the splits vector.
           * @returns A mutable reference to the splits vector.
           */
          std::vector<Split>& getSplits() noexcept
          {
            return m_splits;
          }

          /**
           * @brief Returns a read-only reference to the splits vector.
           * @returns A read-only reference to the splits vector.
           */
          const std::vector<Split>& getSplits() const noexcept
          {
            return m_splits;
          }

          /// @brief Destroys all owned PETSc index sets and clears the splits.
          void clear() noexcept
          {
            [[maybe_unused]] PetscErrorCode ierr;
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
          }

          /**
           * @brief Deep-copies splits from another FieldSplits object.
           * @param[in] other Source FieldSplits to copy from.
           */
          void copyFrom(const FieldSplits& other)
          {
            [[maybe_unused]] PetscErrorCode ierr;

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

          /**
           * @brief Returns the number of splits.
           * @returns The number of splits.
           */
          size_t size() const noexcept
          {
            return m_splits.size();
          }

          /**
           * @brief Returns true if there are no splits.
           * @returns True if there are no splits.
           */
          bool empty() const noexcept
          {
            return m_splits.empty();
          }

          /**
           * @brief Returns the name under which split @p i is installed.
           * @param[in] i Split index.
           * @returns The split name, or its index when it has none.
           */
          std::string getName(size_t i) const
          {
            return m_splits[i].name.empty() ? std::to_string(i) : m_splits[i].name;
          }

          /**
           * @brief Installs the splits on a field-split preconditioner.
           *
           * PETSc appends every split set in code, even one it already
           * holds, so a preconditioner that holds these splits is left as
           * is. A change of names, order, count or indices resets the blocks
           * before installing the new configuration. Empty stored splits
           * leave runtime-defined splits under PETSc control.
           * @param[in] pc Preconditioner of type `PCFIELDSPLIT`.
           * @returns `PETSC_SUCCESS`, or the first PETSc error.
           */
          PetscErrorCode install(::PC pc) const
          {
            if (m_splits.empty())
              return PETSC_SUCCESS;
            // PETSc has no pre-setup split-count query. The indexed query
            // reports an out-of-range error when no extra split exists.
            PetscErrorCode ierr = PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
            if (ierr != PETSC_SUCCESS)
              return ierr;
            ::IS extra = PETSC_NULLPTR;
            const PetscErrorCode probe = PCFieldSplitGetISByIndex(
              pc, static_cast<PetscInt>(m_splits.size()), &extra);
            ierr = PetscPopErrorHandler();
            if (ierr != PETSC_SUCCESS)
              return ierr;
            if (probe != PETSC_SUCCESS && probe != PETSC_ERR_ARG_OUTOFRANGE)
              return probe;
            bool installed = probe == PETSC_ERR_ARG_OUTOFRANGE;
            for (size_t k = 0; installed && k < m_splits.size(); ++k)
            {
              ::IS named = PETSC_NULLPTR, indexed = PETSC_NULLPTR;
              ierr = PCFieldSplitGetIS(pc, getName(k).c_str(), &named);
              if (ierr != PETSC_SUCCESS)
                return ierr;
              if (!named)
              {
                installed = false;
                break;
              }
              ierr = PCFieldSplitGetISByIndex(pc, static_cast<PetscInt>(k), &indexed);
              if (ierr != PETSC_SUCCESS)
                return ierr;
              PetscBool equal = PETSC_FALSE;
              ierr = ISEqual(named, m_splits[k].is, &equal);
              if (ierr != PETSC_SUCCESS)
                return ierr;
              installed = named == indexed && equal;
            }
            if (installed)
              return PETSC_SUCCESS;
            {
              // PCReset also releases the operators, which are set again.
              PetscBool hasOperator = PETSC_FALSE, hasPreconditioner = PETSC_FALSE;
              ierr = PCGetOperatorsSet(pc, &hasOperator, &hasPreconditioner);
              if (ierr != PETSC_SUCCESS)
                return ierr;
              ::Mat operators[2] = {PETSC_NULLPTR, PETSC_NULLPTR};
              if (hasOperator || hasPreconditioner)
              {
                ierr = PCGetOperators(pc, &operators[0], &operators[1]);
                if (ierr != PETSC_SUCCESS)
                  return ierr;
                for (::Mat m : operators)
                {
                  ierr = PetscObjectReference(reinterpret_cast<PetscObject>(m));
                  if (ierr != PETSC_SUCCESS)
                    return ierr;
                }
              }
              ierr = PCReset(pc);
              if (ierr == PETSC_SUCCESS && (hasOperator || hasPreconditioner))
                ierr = PCSetOperators(pc, operators[0], operators[1]);
              for (::Mat& m : operators)
              {
                const PetscErrorCode destroyed = MatDestroy(&m);
                if (ierr == PETSC_SUCCESS)
                  ierr = destroyed;
              }
              if (ierr != PETSC_SUCCESS)
                return ierr;
            }
            for (size_t k = 0; k < m_splits.size(); ++k)
            {
              assert(m_splits[k].is);
              ierr = PCFieldSplitSetIS(pc, getName(k).c_str(), m_splits[k].is);
              if (ierr != PETSC_SUCCESS)
                return ierr;
            }
            return PETSC_SUCCESS;
          }

        private:
          std::vector<Split> m_splits; ///< The stored field splits.
      };

      /**
       * @brief Constructs an empty linear system on the given communicator.
       * @param comm MPI communicator for distributed PETSc objects.
       */
      LinearSystem(MPI_Comm comm);

      /**
       * @brief Copy constructor (increments PETSc handle reference counts).
       * @param other Object to copy from.
       */
      LinearSystem(const LinearSystem& other);

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      LinearSystem(LinearSystem&& other) noexcept;

      /// @brief Destructor; destroys owned PETSc matrix and vector handles.
      virtual ~LinearSystem();

      /**
       * @brief Copy-assignment operator.
       * @param other Object to copy from.
       * @returns Reference to this object after the operation.
       */
      LinearSystem& operator=(const LinearSystem& other);

      /**
       * @brief Move-assignment operator.
       * @param other Object to move from.
       * @returns Reference to this object after the operation.
       */
      LinearSystem& operator=(LinearSystem&& other) noexcept;

      /**
       * @brief Enforces Dirichlet boundary conditions by eliminating
       *        constrained DOFs from the linear system.
       *
       * For every entry @f$ (i, u_i) @f$ in @p dofs, replaces row @f$i@f$
       * by the unit diagonal row and sets @f$b_i=u_i@f$. Unconstrained rows
       * and constrained columns retain their original entries. This row-only
       * operation gives the prescribed solution for arbitrary square operators,
       * but generally does not preserve matrix symmetry. The solution vector is
       * not changed before the solve.
       *
       * Canonical variational assembly performs row-and-column elimination with
       * right-hand-side lifting, preserving symmetry when the original operator
       * is symmetric. The low-level PETSc operation documented here follows
       * @c MatZeroRows instead.
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

        Vec bc = PETSC_NULLPTR;
        ierr = VecDuplicate(x, &bc);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecCopy(x, bc);
        assert(ierr == PETSC_SUCCESS);

        for (auto const& kv : dofs)
        {
          const PetscInt i = kv.first + offset;
          const auto& ui = kv.second;
          ierr = VecSetValue(bc, i, ui, INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);
        }

        ierr = VecAssemblyBegin(bc);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyEnd(bc);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetOption(a, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroRows(a, rows.size(), rows.data(), 1.0, bc, b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecDestroy(&bc);
        assert(ierr == PETSC_SUCCESS);

        (void)ierr;

        return *this;
      }

      /**
       * @brief Merges periodic DOF constraints into the linear system.
       * @note Currently unimplemented — throws at runtime.
       * @param offset Offset in the indexed data.
       * @returns Reference to this object after the operation.
       * @param dofs DOF-to-functional map to merge.
       */
      template <class DOFScalar>
      LinearSystemBase& merge(
        const IndexMap<std::pair<IndexArray, Math::Vector<DOFScalar>>>& dofs,
        size_t offset = 0)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Periodic DOF constraint merging is not implemented." << Alert::Raise;
        return *this;
      }

      /**
       * @brief Returns the MPI communicator of this linear system.
       * @returns The MPI communicator of this linear system.
       */
      constexpr MPI_Comm getCommunicator() const noexcept
      {
        return m_comm;
      }

      /**
       * @brief Returns a mutable reference to the system matrix @f$ A @f$.
       * @returns A mutable reference to the system matrix @f$ A @f$.
       */
      constexpr MatrixType& getOperator()
      {
        return m_operator;
      }

      /**
       * @brief Returns a read-only reference to the system matrix @f$ A @f$.
       * @returns A read-only reference to the system matrix @f$ A @f$.
       */
      constexpr const MatrixType& getOperator() const
      {
        return m_operator;
      }

      /**
       * @brief Returns a mutable reference to the right-hand side vector @f$ b @f$.
       * @returns A mutable reference to the right-hand side vector @f$ b @f$.
       */
      constexpr VectorType& getVector()
      {
        return m_vector;
      }

      /**
       * @brief Returns a read-only reference to the right-hand side vector @f$ b @f$.
       * @returns A read-only reference to the right-hand side vector @f$ b @f$.
       */
      constexpr const VectorType& getVector() const
      {
        return m_vector;
      }

      /**
       * @brief Returns a mutable reference to the solution vector @f$ x @f$.
       * @returns A mutable reference to the solution vector @f$ x @f$.
       */
      constexpr VectorType& getSolution()
      {
        return m_solution;
      }

      /**
       * @brief Returns a read-only reference to the solution vector @f$ x @f$.
       * @returns A read-only reference to the solution vector @f$ x @f$.
       */
      constexpr const VectorType& getSolution() const
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

      /**
       * @brief Returns the current field-split configuration.
       * @returns The current field-split configuration.
       */
      const FieldSplits& getFieldSplits() const noexcept
      {
        return m_fieldSplits;
      }

      /**
       * @brief Sets the first row of each field on this process.
       *
       * An assembly sets them when it lays out the owned rows of a process
       * field after field, instead of numbering each field from its offset
       * in the system. The owned degrees of freedom of field @f$ k @f$ are
       * then the rows from @p starts[k] on, in the order of its space.
       * @param starts First row of each field, empty for the offsets.
       * @returns Reference to `*this`.
       */
      LinearSystem& setFieldStarts(std::vector<PetscInt> starts)
      {
        m_fieldStarts = std::move(starts);
        return *this;
      }

      /**
       * @brief Gets the first row of each field on this process.
       * @returns First row of each field, empty when each field is numbered
       * from its offset in the system.
       */
      const std::vector<PetscInt>& getFieldStarts() const noexcept
      {
        return m_fieldStarts;
      }

    private:
      MPI_Comm m_comm; ///< The MPI communicator for the linear system.
      MatrixType m_operator; ///< The operator of the linear system.
      VectorType m_solution; ///< The solution vector of the linear system.
      VectorType m_vector; ///< The vector of the linear system.

      FieldSplits m_fieldSplits; ///< The field splits for block preconditioning.

      std::vector<PetscInt> m_fieldStarts; ///< First row of each field on this process.
  };
}

namespace Rodin::PETSc::Math
{
  /**
   * @brief Convenient alias for the PETSc linear system specialization.
   *
   * @see <a href="class_rodin_1_1_math_1_1_linear_system.html">Rodin::Math::LinearSystem</a>
   */
  using LinearSystem = Rodin::Math::LinearSystem<Matrix, Vector>;
}

#endif
