#ifndef RODIN_PETSC_ASSEMBLY_MPI_H
#define RODIN_PETSC_ASSEMBLY_MPI_H

/**
 * @file
 * @brief MPI assembly specializations targeting PETSc objects.
 *
 * In this module, mesh entity pairs `(d, i)` are shard-local indices. PETSc
 * assembly iterates over rank-local MPI mesh polytopes, filters owned entities
 * for rows/contributions, and passes the same local pair to the MPI finite
 * element spaces and integrators. Global bilinear forms keep the test entity
 * owned, but may use ghost trial entities for off-process columns; those trial
 * indices are still local to the rank shard.
 */

#include <petsc.h>
#include <petscmacros.h>
#include <petscmat.h>
#include <cassert>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/MPI/Assembly.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Math/Matrix.h"

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/Assembly/MatrixSetup.h"
#include "Rodin/PETSc/Assembly/VectorSetup.h"

namespace Rodin::Assembly
{
  /**
   * @brief MPI-parallel assembly of a PETSc vector from a linear form.
   *
   * Each MPI rank assembles contributions from its owned polytopes into
   * the distributed PETSc vector.
   *
   * @tparam FES Finite element space type.
   */
  template <class FES>
  class MPI<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      /// @brief Scalar type of the DOF coefficients (`PetscScalar`).
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar"
      );
      /// @brief PETSc vector type (`::Vec`).
      using VectorType     = ::Vec;
      /// @brief Linear form type being assembled.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      /// @brief Parent assembly base class.
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType      = typename Parent::InputType;

      void execute(VectorType& res, const InputType& input) const override
      {
        assert(res);
        const auto& fes = input.getFES();
        const auto& mesh = fes.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx = mesh.getContext();
        const size_t globalSize = fes.getSize();

        size_t begin, end;
        fes.getOwnershipRange(begin, end);
        const size_t ownedSize = end - begin;

        PetscErrorCode ierr;
        ierr = PETSc::Assembly::VectorSetup(res).prepare({
          static_cast<PetscInt>(ownedSize),
          static_cast<PetscInt>(globalSize),
          nullptr,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          MPIIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index idx = it->getIndex();

            if (!shard.isOwned(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            lfi.setPolytope(*it);
            const auto& dofs = fes.getDOFs(d, idx);

            for (PetscInt i = 0; i < static_cast<PetscInt>(dofs.size()); ++i)
            {
              const PetscInt r = static_cast<PetscInt>(dofs[i]);
              const PetscScalar v = lfi.integrate(i);
              ierr = VecSetValue(res, r, v, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }
        }

        ierr = VecAssemblyBegin(res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ierr = VecAssemblyEnd(res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        (void) ierr;
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  /**
   * @brief MPI-parallel assembly of a PETSc matrix from a bilinear form.
   *
   * Each MPI rank assembles contributions from its owned polytopes into
   * the distributed PETSc matrix.
   *
   * @tparam Solution Solution type.
   * @tparam TrialFES Trial finite element space type.
   * @tparam TestFES  Test finite element space type.
   */
  template <class Solution, class TrialFES, class TestFES>
  class MPI<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      /// @brief Scalar type resulting from the dot product of trial and test scalars.
      using DotType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "DotType must be PetscScalar"
      );
      /// @brief PETSc matrix type (`::Mat`).
      using OperatorType     = ::Mat;
      /// @brief Bilinear form type being assembled.
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;
      /// @brief Parent assembly base class.
      using Parent           = AssemblyBase<OperatorType, BilinearFormType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType        = typename Parent::InputType;

      void execute(OperatorType& res, const InputType& input) const override
      {
        assert(res);

        const auto& trialFES = input.getTrialFES();
        const auto& testFES  = input.getTestFES();
        assert(trialFES.getMesh() == testFES.getMesh());
        const auto& mesh = trialFES.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx = mesh.getContext();

        size_t rbegin, rend;
        testFES.getOwnershipRange(rbegin, rend);
        const size_t localRows = rend - rbegin;

        size_t cbegin, cend;
        trialFES.getOwnershipRange(cbegin, cend);
        const size_t localCols = cend - cbegin;

        const size_t globalRows = testFES.getSize();

        const size_t globalCols = trialFES.getSize();

        PetscErrorCode ierr;
        ierr = PETSc::Assembly::MatrixSetup(res).prepare({
          static_cast<PetscInt>(localRows),
          static_cast<PetscInt>(localCols),
          static_cast<PetscInt>(globalRows),
          static_cast<PetscInt>(globalCols),
          nullptr,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MPIIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index idx = it->getIndex();

            if (!shard.isOwned(d, idx))
              continue;

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            bfi.setPolytope(*it);
            const auto& rows = testFES.getDOFs(d, idx);
            const auto& cols = trialFES.getDOFs(d, idx);

            for (Index i = 0; i < rows.size(); ++i)
            {
              const PetscInt r = static_cast<PetscInt>(rows[i]);
              for (Index j = 0; j < cols.size(); ++j)
              {
                const PetscInt c = static_cast<PetscInt>(cols[j]);
                const PetscScalar v = bfi.integrate(j, i);
                ierr = MatSetValue(res, r, c, v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
                (void) ierr;
              }
            }
          }
        }

        ierr = MatAssemblyBegin(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ierr = MatAssemblyEnd(res, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        (void) ierr;
      }

      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };

  /**
   * @brief MPI-parallel assembly of a single-variable PETSc problem.
   *
   * Each MPI rank assembles its local contributions to the system matrix
   * and right-hand side vector of a @ref Rodin::Variational::Problem
   * backed by PETSc distributed objects.
   *
   * @tparam U Trial function type.
   * @tparam V Test function type.
   */
  // MPI assembly for single-variable Problem (PETSc)
  template <class U, class V>
  class MPI<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>>
  {
    public:
      /// @brief PETSc linear system type.
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      /// @brief Problem type being assembled.
      using ProblemType      = Rodin::Variational::Problem<LinearSystemType, U, V>;
      /// @brief Parent assembly base class.
      using Parent           = AssemblyBase<LinearSystemType, ProblemType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType        = typename Parent::InputType;

      /// @brief PETSc matrix type (`::Mat`) for the system operator.
      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      /// @brief PETSc vector type (`::Vec`) for the RHS and solution.
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      /// @brief Finite element space type for the trial function.
      using TrialFESType  = typename Rodin::FormLanguage::Traits<U>::FESType;
      /// @brief Mesh type for the trial finite element space.
      using TrialMeshType = typename Rodin::FormLanguage::Traits<TrialFESType>::MeshType;
      /// @brief Context type (Local or MPI) for the trial mesh.
      using MeshContextType =
        typename Rodin::FormLanguage::Traits<TrialMeshType>::ContextType;

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        execute(axb, input, AssemblyMode::Full);
      }

      void execute(
          LinearSystemType& axb,
          const InputType& input,
          Rodin::Variational::AssemblyTarget target) const
      {
        switch (target)
        {
          case Rodin::Variational::AssemblyTarget::LHS:
            execute(axb, input, AssemblyMode::LHS);
            break;
          case Rodin::Variational::AssemblyTarget::RHS:
            execute(axb, input, AssemblyMode::RHS);
            break;
        }
      }

    private:
      enum class AssemblyMode
      {
        Full,
        LHS,
        RHS
      };

      void execute(
          LinearSystemType& axb,
          const InputType& input,
          AssemblyMode mode) const
      {
        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::MPI>,
          "PETSc MPI assembly should only be used with MPI mesh context.");

        const bool doMatrix = mode != AssemblyMode::RHS;
        const bool doVector = mode != AssemblyMode::LHS;

        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();

        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES  = v.getFiniteElementSpace();
        assert(trialFES.getMesh() == testFES.getMesh());

        const auto& mesh  = trialFES.getMesh();
        const auto& shard = mesh.getShard();
        const auto& ctx   = mesh.getContext();

        const size_t globalCols = trialFES.getSize();
        const size_t globalRows = testFES.getSize();

        // Ownership ranges (global indices)
        size_t rbegin, rend;
        testFES.getOwnershipRange(rbegin, rend);
        const size_t localRows = rend - rbegin;

        size_t cbegin, cend;
        trialFES.getOwnershipRange(cbegin, cend);
        const size_t localCols = cend - cbegin;

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (MPIAIJ); re-use structure across assemblies.
        // ------------------------
        assert(A);
        if (doMatrix)
        {
          ierr = PETSc::Assembly::MatrixSetup(A).prepare({
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(globalRows),
            static_cast<PetscInt>(globalCols),
            nullptr,
            true
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // ------------------------
        // Allocate / reset b (MPI Vec)
        // ------------------------
        assert(b);
        if (doVector)
        {
          ierr = PETSc::Assembly::VectorSetup(b).prepare({
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(globalRows),
            nullptr,
            true,
            true
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        auto& x = axb.getSolution();
        assert(x);
        // The solution vector carries the previous iterate as the solver's
        // initial guess, so it must not be zeroed on reuse.
        ierr = PETSc::Assembly::VectorSetup(x).prepare({
          static_cast<PetscInt>(localCols),
          static_cast<PetscInt>(globalCols),
          nullptr,
          true,
          false
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ConstraintMap<PetscScalar> constraints(std::max(globalRows, globalCols));
        using DBCBaseType = Variational::DirichletBCBase<PetscScalar>;
        using ValueDOFsType = typename DBCBaseType::ValueDOFs;
        using IdentDOFsType = typename DBCBaseType::IdentifiedDOFs;

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;
          dbc.assemble();
          std::visit([&](auto&& dofs)
          {
            using T = std::decay_t<decltype(dofs)>;
            if constexpr (std::is_same_v<T, ValueDOFsType>)
            {
              for (const auto& [local, value] : dofs)
                constraints.setFixed(
                    static_cast<Index>(local),
                    static_cast<PetscScalar>(value));
            }
            else if constexpr (std::is_same_v<T, IdentDOFsType>)
            {
              const auto& affineValues = dbc.getIdentificationValues();
              for (const auto& [slave, pair] : dofs)
              {
                const auto& masters = pair.first;
                const auto& coeffs = pair.second;
                std::vector<typename ConstraintMap<PetscScalar>::Entry> entries;
                entries.reserve(static_cast<size_t>(masters.size()));
                for (Index k = 0; k < masters.size(); k++)
                {
                  entries.push_back({
                      static_cast<Index>(masters[k]),
                      static_cast<PetscScalar>(coeffs[k]) });
                }
                const auto valueIt = affineValues.find(slave);
                const PetscScalar value =
                  valueIt == affineValues.end()
                    ? PetscScalar(0)
                    : static_cast<PetscScalar>(valueIt->second);
                constraints.setIdentification(static_cast<Index>(slave), entries, value);
              }
            }
          }, dbc.getDOFs());
        }

        if (mode != AssemblyMode::Full && !constraints.getIdentifiedRows().empty())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Targeted assembly is not implemented for identification DirichletBCs."
            << Alert::Raise;
        }

        auto matrix_entry = [&](Index row, Index col, PetscScalar val)
        {
          const PetscScalar colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : PetscScalar(0);
          for (const auto& r : constraints.expand(row))
          {
            const PetscInt I = static_cast<PetscInt>(r.index);
            if (colValue != PetscScalar(0))
            {
              const PetscScalar rhsShift = -r.coefficient * val * colValue;
              if (doVector)
              {
                ierr = VecSetValue(b, I, rhsShift, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
                (void) ierr;
              }
            }
            if (doMatrix)
            {
              for (const auto& c : constraints.expand(col))
              {
                const PetscInt J = static_cast<PetscInt>(c.index);
                const PetscScalar v =
                  r.coefficient * val * c.coefficient;
                // Insert unconditionally. PETSc decides whether zero-valued
                // insertions affect the pattern according to its current
                // matrix options.
                ierr = MatSetValue(A, I, J, v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
                (void) ierr;
              }
            }
          }
        };

        auto vector_entry = [&](Index row, PetscScalar val)
        {
          if (val == PetscScalar(0))
            return;
          for (const auto& r : constraints.expand(row))
          {
            const PetscInt I = static_cast<PetscInt>(r.index);
            const PetscScalar v = r.coefficient * val;
            ierr = VecSetValue(b, I, v, ADD_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        };

        // ------------------------
        // Local BFIs (owned elements only)
        // ------------------------
        if (doMatrix)
        {
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto& attrs = bfi.getAttributes();
            MPIIteration seq(mesh, bfi.getRegion());

            for (auto it = seq.getIterator(); it; ++it)
            {
              const size_t d   = it->getDimension();
              const Index  idx = it->getIndex();

              if (!shard.isOwned(d, idx))
                continue;

              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              bfi.setPolytope(*it);

              const auto& rowsDOF = testFES.getDOFs(d, idx);
              const auto& colsDOF = trialFES.getDOFs(d, idx);

              for (Index i = 0; i < rowsDOF.size(); ++i)
              {
                const PetscInt I = rowsDOF[i];
                for (Index j = 0; j < colsDOF.size(); ++j)
                {
                  const PetscInt J = colsDOF[j];
                  const PetscScalar val = bfi.integrate(j, i);
                  matrix_entry(I, J, val);
                }
              }
            }
          }

          // ------------------------
          // Global BFIs (owned test entities only)
          // - This avoids double assembly across ranks.
          // - We allow trial entities to be ghost; PETSc handles off-proc columns.
          // ------------------------
          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();

            MPIIteration trialseq(mesh, bfi.getTrialRegion());
            MPIIteration testseq(mesh,  bfi.getTestRegion());

            for (auto teIt = testseq.getIterator(); teIt; ++teIt)
            {
              const size_t td   = teIt->getDimension();
              const Index  tidx = teIt->getIndex();

              if (!shard.isOwned(td, tidx))
                continue;

              if (!testAttrs.empty())
              {
                const auto a = teIt->getAttribute();
                if (!a || !testAttrs.count(*a))
                  continue;
              }

              const auto& rowsDOF = testFES.getDOFs(td, tidx);

              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                const size_t rd   = trIt->getDimension();
                const Index  ridx = trIt->getIndex();

                if (!trialAttrs.empty())
                {
                  const auto a = trIt->getAttribute();
                  if (!a || !trialAttrs.count(*a))
                    continue;
                }

                // NOTE: do NOT skip ghost trial entities here: they still contribute to
                // owned test rows, producing off-process columns in MatSetValue.
                const auto& colsDOF = trialFES.getDOFs(rd, ridx);

                bfi.setPolytope(*trIt, *teIt);

                for (Index i = 0; i < rowsDOF.size(); ++i)
                {
                  const PetscInt I = rowsDOF[i];
                  for (Index j = 0; j < colsDOF.size(); ++j)
                  {
                    const PetscInt J = colsDOF[j];
                    const PetscScalar val = bfi.integrate(j, i);
                    matrix_entry(I, J, val);
                  }
                }
              }
            }
          }

          // Preassembled bilinear forms
          for (auto& bf : pb.getBFs())
          {
            const auto& op = bf.getOperator();
            PetscInt opStart, opEnd;
            ierr = MatGetOwnershipRange(op, &opStart, &opEnd);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (PetscInt i = opStart; i < opEnd; ++i)
            {
              PetscInt nc;
              const PetscInt* cols;
              const PetscScalar* vals;
              ierr = MatGetRow(op, i, &nc, &cols, &vals);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
              for (PetscInt j = 0; j < nc; ++j)
                matrix_entry(static_cast<Index>(i), static_cast<Index>(cols[j]), vals[j]);
              ierr = MatRestoreRow(op, i, &nc, &cols, &vals);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }
        }

        // Assemble A
        if (doMatrix)
        {
          ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // ------------------------
        // Linear forms (owned elements only)
        // Convention: matches your sequential single-variable version: b += -LF
        // ------------------------
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
          {
            const auto& attrs = lfi.getAttributes();
            MPIIteration seq(mesh, lfi.getRegion());

            for (auto it = seq.getIterator(); it; ++it)
            {
              const size_t d   = it->getDimension();
              const Index  idx = it->getIndex();

              if (!shard.isOwned(d, idx))
                continue;

              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);

              const auto& dofs = testFES.getDOFs(d, idx);
              for (Index l = 0; l < dofs.size(); ++l)
              {
                const PetscInt I = dofs[l];
                const PetscScalar val = lfi.integrate(l);
                if (val != PetscScalar(0))
                  vector_entry(I, -val);
              }
            }
          }

          // Preassembled linear forms: b += LF
          for (auto& lf : pb.getLFs())
          {
            const auto& vec = lf.getVector();
            PetscInt lo, hi;
            ierr = VecGetOwnershipRange(vec, &lo, &hi);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            const PetscScalar* arr;
            ierr = VecGetArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (PetscInt i = lo; i < hi; ++i)
            {
              const PetscScalar val = arr[i - lo];
              if (val != PetscScalar(0))
                vector_entry(static_cast<Index>(i), val);
            }
            ierr = VecRestoreArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        // Assemble b
        if (doVector)
        {
          ierr = VecAssemblyBegin(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyEnd(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        if (doMatrix)
        {
          std::vector<PetscInt> rowsToZero;
          for (const Index gs : constraints.getIdentifiedRows())
            if (rbegin <= gs && gs < rend)
              rowsToZero.push_back(static_cast<PetscInt>(gs));

          ierr = MatZeroRows(
              A,
              static_cast<PetscInt>(rowsToZero.size()),
              rowsToZero.empty() ? nullptr : rowsToZero.data(),
              0.0,
              nullptr,
              nullptr);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (!(rbegin <= gs && gs < rend))
              continue;
            const PetscInt I = static_cast<PetscInt>(gs);
            const PetscScalar one = 1.0;
            ierr = MatSetValue(A, I, I, one, ADD_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (const auto& e : constraints.expand(gs))
            {
              const PetscInt J = static_cast<PetscInt>(e.index);
              const PetscScalar v = -e.coefficient;
              ierr = MatSetValue(A, I, J, v, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
            const PetscScalar rhs = constraints.getIdentificationValue(gs);
            if (doVector)
            {
              ierr = VecSetValue(b, I, rhs, INSERT_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }

          ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          if (doVector)
          {
            ierr = VecAssemblyBegin(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyEnd(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        std::vector<PetscInt> bcIdx;
        std::vector<PetscScalar> bcVals;
        for (Index i = rbegin; i < rend; i++)
        {
          if (constraints.isFixed(i))
          {
            bcIdx.push_back(static_cast<PetscInt>(i));
            bcVals.push_back(constraints.getFixedValue(i));
          }
        }

        {
          if (mode == AssemblyMode::Full)
          {
            Vec bcVec;
            ierr = VecDuplicate(b, &bcVec);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecZeroEntries(bcVec);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecSetValues(
                bcVec,
                static_cast<PetscInt>(bcIdx.size()),
                bcIdx.empty() ? nullptr : bcIdx.data(),
                bcVals.empty() ? nullptr : bcVals.data(),
                INSERT_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyBegin(bcVec);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyEnd(bcVec);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = MatZeroRowsColumns(
                A,
                static_cast<PetscInt>(bcIdx.size()),
                bcIdx.empty() ? nullptr : bcIdx.data(),
                1.0,
                bcVec,
                b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecDestroy(&bcVec);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
          else if (mode == AssemblyMode::LHS)
          {
            ierr = MatZeroRows(
                A,
                static_cast<PetscInt>(bcIdx.size()),
                bcIdx.empty() ? nullptr : bcIdx.data(),
                1.0,
                nullptr,
                nullptr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
          else
          {
            ierr = VecSetValues(
                b,
                static_cast<PetscInt>(bcIdx.size()),
                bcIdx.empty() ? nullptr : bcIdx.data(),
                bcVals.empty() ? nullptr : bcVals.data(),
                INSERT_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyBegin(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyEnd(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        (void) ierr;
      }

    public:
      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };


  /**
   * @brief MPI-parallel assembly of a multi-variable PETSc problem.
   *
   * Assembles a block-structured @ref Rodin::Variational::Problem backed
   * by PETSc distributed objects across MPI ranks.
   *
   * @tparam U1  First trial/test function type.
   * @tparam U2  Second trial/test function type.
   * @tparam U3  Third trial/test function type.
   * @tparam Us  Additional trial/test function types.
   */
  // MPI assembly for multi-variable Problem (PETSc)
  template <class U1, class U2, class U3, class ... Us>
  class MPI<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      /// @brief PETSc linear system type.
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      /// @brief Multi-field problem type being assembled.
      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;
      /// @brief Parent assembly base class.
      using Parent    = AssemblyBase<LinearSystemType, ProblemType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType = typename Parent::InputType;

      /// @brief PETSc matrix type (`::Mat`) for the block system operator.
      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      /// @brief PETSc vector type (`::Vec`) for the block RHS and solution.
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        execute(axb, input, AssemblyMode::Full);
      }

      void execute(
          LinearSystemType& axb,
          const InputType& input,
          Rodin::Variational::AssemblyTarget target) const
      {
        switch (target)
        {
          case Rodin::Variational::AssemblyTarget::LHS:
            execute(axb, input, AssemblyMode::LHS);
            break;
          case Rodin::Variational::AssemblyTarget::RHS:
            execute(axb, input, AssemblyMode::RHS);
            break;
        }
      }

    private:
      enum class AssemblyMode
      {
        Full,
        LHS,
        RHS
      };

      void execute(
          LinearSystemType& axb,
          const InputType& input,
          AssemblyMode mode) const
      {
        auto& A = axb.getOperator();
        auto& b = axb.getVector();
        const bool doMatrix = mode != AssemblyMode::RHS;
        const bool doVector = mode != AssemblyMode::LHS;

        auto& pb            = input.getProblemBody();
        auto& us            = input.getTrialFunctions();
        auto& vs            = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto& trialUUIDMap  = input.getTrialUUIDMap();
        auto& testUUIDMap   = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        using FirstTrialMeshType =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        using MeshContextType =
          typename Rodin::FormLanguage::Traits<FirstTrialMeshType>::ContextType;

        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::MPI>,
          "PETSc MPI assembly should only be used with MPI mesh context.");

        // Mesh + shard
        const auto& mesh = [&]() -> const auto&
        {
          const void* addr = nullptr;
          us.apply([&](const auto& uref)
          {
            if (!addr)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace().getMesh());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>*>(addr);
        }();

        const auto& shard = mesh.getShard();

        // Ownership range for the *global* test vector (assembled rows)
        // Here we assume "ownership range" refers to the assembled global test numbering (with offsets applied),
        // hence we use input.getTotalTestSize() and distribute it by PETSc using Vec/Mat layout already set by caller,
        // OR you keep FES-based ownership and accept that blocks can be owned by different ranks.
        // In your code base you already have FES.getOwnershipRange for each field; we follow that and let PETSc handle.
        //
        // Layout: localRows/localCols computed as sum of owned sizes per block.

        // Helpers (same as your sequential)
        const auto findTrialBlock = [&](const auto& uuid) -> size_t
        {
          auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return it->second;
        };

        const auto findTestBlock = [&](const auto& uuid) -> size_t
        {
          auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return it->second;
        };

        const auto& getTrialFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTrialBlock(uuid);
          const void* addr = nullptr;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        const auto& getTestFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTestBlock(uuid);
          const void* addr = nullptr;
          vs.iapply([&](size_t i, const auto& vref)
          {
            if (i == k)
              addr = static_cast<const void*>(&vref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(vs.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        // Compute localRows/localCols as sum of owned sizes of each block on this rank
        size_t localRows = 0;
        vs.iapply([&](size_t /*i*/, const auto& vref)
        {
          size_t rb, re;
          vref.get().getFiniteElementSpace().getOwnershipRange(rb, re);
          localRows += (re - rb);
        });

        size_t localCols = 0;
        us.iapply([&](size_t /*i*/, const auto& uref)
        {
          size_t cb, ce;
          uref.get().getFiniteElementSpace().getOwnershipRange(cb, ce);
          localCols += (ce - cb);
        });

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (MPIAIJ); re-use structure across assemblies.
        // ------------------------
        assert(A);
        if (doMatrix)
        {
          ierr = PETSc::Assembly::MatrixSetup(A).prepare({
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(localCols),
            static_cast<PetscInt>(nrows),
            static_cast<PetscInt>(ncols),
            nullptr,
            true
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // ------------------------
        // Allocate / reset b (MPI Vec)
        // ------------------------
        assert(b);
        if (doVector)
        {
          ierr = PETSc::Assembly::VectorSetup(b).prepare({
            static_cast<PetscInt>(localRows),
            static_cast<PetscInt>(nrows),
            nullptr,
            true,
            true
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        auto& x = axb.getSolution();
        assert(x);
        // The solution vector carries the previous iterate as the solver's
        // initial guess, so it must not be zeroed on reuse.
        ierr = PETSc::Assembly::VectorSetup(x).prepare({
          static_cast<PetscInt>(localCols),
          static_cast<PetscInt>(ncols),
          nullptr,
          true,
          false
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ConstraintMap<PetscScalar> constraints(std::max(nrows, ncols));
        using DBCBaseType = Variational::DirichletBCBase<PetscScalar>;
        using ValueDOFsType = typename DBCBaseType::ValueDOFs;
        using IdentDOFsType = typename DBCBaseType::IdentifiedDOFs;

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff = trialOffsets[uBlock];

          dbc.assemble();
          std::visit([&](auto&& dofs)
          {
            using T = std::decay_t<decltype(dofs)>;
            if constexpr (std::is_same_v<T, ValueDOFsType>)
            {
              for (const auto& [local, value] : dofs)
                constraints.setFixed(
                    static_cast<Index>(uOff + static_cast<size_t>(local)),
                    static_cast<PetscScalar>(value));
            }
            else if constexpr (std::is_same_v<T, IdentDOFsType>)
            {
              const auto vUUIDOpt = dbc.getValueUUID();
              assert(vUUIDOpt);
              const size_t vBlock = findTrialBlock(*vUUIDOpt);
              const size_t vOff = trialOffsets[vBlock];
              const auto& affineValues = dbc.getIdentificationValues();
              for (const auto& [slave, pair] : dofs)
              {
                const auto& masters = pair.first;
                const auto& coeffs = pair.second;
                std::vector<typename ConstraintMap<PetscScalar>::Entry> entries;
                entries.reserve(static_cast<size_t>(masters.size()));
                for (Index k = 0; k < masters.size(); k++)
                {
                  entries.push_back({
                      static_cast<Index>(vOff + static_cast<size_t>(masters[k])),
                      static_cast<PetscScalar>(coeffs[k]) });
                }
                const auto valueIt = affineValues.find(slave);
                const PetscScalar value =
                  valueIt == affineValues.end()
                    ? PetscScalar(0)
                    : static_cast<PetscScalar>(valueIt->second);
                constraints.setIdentification(
                    static_cast<Index>(uOff + static_cast<size_t>(slave)),
                    entries,
                    value);
              }
            }
          }, dbc.getDOFs());
        }

        if (mode != AssemblyMode::Full && !constraints.getIdentifiedRows().empty())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Targeted assembly is not implemented for identification DirichletBCs."
            << Alert::Raise;
        }

        auto ownsTrialGlobal = [&](Index g) -> bool
        {
          bool owned = false;
          us.iapply([&](size_t i, const auto& uref)
          {
            size_t begin, end;
            uref.get().getFiniteElementSpace().getOwnershipRange(begin, end);
            const size_t off = trialOffsets[i];
            if (off + begin <= g && g < off + end)
              owned = true;
          });
          return owned;
        };

        auto matrix_entry = [&](Index row, Index col, PetscScalar val)
        {
          const PetscScalar colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : PetscScalar(0);
          for (const auto& r : constraints.expand(row))
          {
            const PetscInt I = static_cast<PetscInt>(r.index);
            if (colValue != PetscScalar(0))
            {
              const PetscScalar rhsShift = -r.coefficient * val * colValue;
              ierr = VecSetValue(b, I, rhsShift, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
            for (const auto& c : constraints.expand(col))
            {
              const PetscInt J = static_cast<PetscInt>(c.index);
              const PetscScalar v =
                r.coefficient * val * c.coefficient;
              // Insert unconditionally. PETSc decides whether zero-valued
              // insertions affect the pattern according to its current
              // matrix options.
              ierr = MatSetValue(A, I, J, v, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }
        };

        auto vector_entry = [&](Index row, PetscScalar val)
        {
          if (val == PetscScalar(0))
            return;
          for (const auto& r : constraints.expand(row))
          {
            const PetscInt I = static_cast<PetscInt>(r.index);
            const PetscScalar v = r.coefficient * val;
            ierr = VecSetValue(b, I, v, ADD_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        };

        // ------------------------
        // Assemble bilinear terms into A
        // - For each integrator, only assemble contributions associated to an owned "assembly entity":
        //   - local BFIs: owned cell/facet/etc (same as your MPI Mat/BFI specialization)
        //   - global BFIs: owned test entity (rows owned), trial may be ghost
        // ------------------------

        if (doMatrix)
        {
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto uUUID = bfi.getTrialFunction().getUUID();
            const auto vUUID = bfi.getTestFunction().getUUID();

            const size_t uBlock = findTrialBlock(uUUID);
            const size_t vBlock = findTestBlock(vUUID);

            const size_t uOff = trialOffsets[uBlock];
            const size_t vOff = testOffsets[vBlock];

            const auto& uFES = getTrialFESByUUID(uUUID);
            const auto& vFES = getTestFESByUUID(vUUID);

            const auto& attrs = bfi.getAttributes();
            MPIIteration seq(mesh, bfi.getRegion());

            for (auto it = seq.getIterator(); it; ++it)
            {
              const size_t d   = it->getDimension();
              const Index  idx = it->getIndex();

              if (!shard.isOwned(d, idx))
                continue;

              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              bfi.setPolytope(*it);

              const auto& rows = vFES.getDOFs(d, idx);
              const auto& cols = uFES.getDOFs(d, idx);

              for (Index i = 0; i < rows.size(); ++i)
              {
                const PetscInt I = vOff + rows[i];
                for (Index j = 0; j < cols.size(); ++j)
                {
                  const PetscInt J = uOff + cols[j];
                  const PetscScalar val = bfi.integrate(j, i);
                  matrix_entry(I, J, val);
                }
              }
            }
          }

          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto uUUID = bfi.getTrialFunction().getUUID();
            const auto vUUID = bfi.getTestFunction().getUUID();

            const size_t uBlock = findTrialBlock(uUUID);
            const size_t vBlock = findTestBlock(vUUID);

            const size_t uOff = trialOffsets[uBlock];
            const size_t vOff = testOffsets[vBlock];

            const auto& uFES = getTrialFESByUUID(uUUID);
            const auto& vFES = getTestFESByUUID(vUUID);

            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();

            MPIIteration trialseq(mesh, bfi.getTrialRegion());
            MPIIteration testseq(mesh,  bfi.getTestRegion());

            for (auto teIt = testseq.getIterator(); teIt; ++teIt)
            {
              const size_t td   = teIt->getDimension();
              const Index  tidx = teIt->getIndex();

              if (!shard.isOwned(td, tidx))
                continue;

              if (!testAttrs.empty())
              {
                const auto a = teIt->getAttribute();
                if (!a || !testAttrs.count(*a))
                  continue;
              }

              const auto& rows = vFES.getDOFs(td, tidx);

              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                const size_t rd   = trIt->getDimension();
                const Index  ridx = trIt->getIndex();

                if (!trialAttrs.empty())
                {
                  const auto a = trIt->getAttribute();
                  if (!a || !trialAttrs.count(*a))
                    continue;
                }

                // do not skip ghost trial entity: off-proc columns are fine
                const auto& cols = uFES.getDOFs(rd, ridx);

                bfi.setPolytope(*trIt, *teIt);

                for (Index i = 0; i < rows.size(); ++i)
                {
                  const PetscInt I = vOff + rows[i];
                  for (Index j = 0; j < cols.size(); ++j)
                  {
                    const PetscInt J = uOff + cols[j];
                    const PetscScalar val = bfi.integrate(j, i);
                    matrix_entry(I, J, val);
                  }
                }
              }
            }
          }

          // Preassembled bilinear forms (with block offsets)
          for (auto& bf : pb.getBFs())
          {
            const auto uUUID = bf.getTrialFunction().getUUID();
            const auto vUUID = bf.getTestFunction().getUUID();

            const size_t uBlock = findTrialBlock(uUUID);
            const size_t vBlock = findTestBlock(vUUID);

            const size_t uOff = trialOffsets[uBlock];
            const size_t vOff = testOffsets[vBlock];

            const auto& op = bf.getOperator();
            PetscInt rStart, rEnd;
            ierr = MatGetOwnershipRange(op, &rStart, &rEnd);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;

            for (PetscInt i = rStart; i < rEnd; ++i)
            {
              PetscInt nc;
              const PetscInt* cols;
              const PetscScalar* vals;
              ierr = MatGetRow(op, i, &nc, &cols, &vals);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
              for (PetscInt j = 0; j < nc; ++j)
              {
                matrix_entry(
                  static_cast<Index>(vOff) + static_cast<Index>(i),
                  static_cast<Index>(uOff) + static_cast<Index>(cols[j]),
                  vals[j]);
              }
              ierr = MatRestoreRow(op, i, &nc, &cols, &vals);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }

          // Assemble A
          ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // ------------------------
        // Assemble linear terms into b
        // Convention: negate each LFI value, consistent with all other assembly paths.
        // ------------------------
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
          {
            const auto vUUID = lfi.getTestFunction().getUUID();
            const size_t vBlock = findTestBlock(vUUID);
            const size_t vOff   = testOffsets[vBlock];

            const auto& vFES = getTestFESByUUID(vUUID);

            const auto& attrs = lfi.getAttributes();
            MPIIteration seq(mesh, lfi.getRegion());

            for (auto it = seq.getIterator(); it; ++it)
            {
              const size_t d   = it->getDimension();
              const Index  idx = it->getIndex();

              if (!shard.isOwned(d, idx))
                continue;

              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);

              const auto& dofs = vFES.getDOFs(d, idx);
              for (Index l = 0; l < dofs.size(); ++l)
              {
                const PetscInt I = vOff + dofs[l];
                const PetscScalar val = lfi.integrate(l);
                if (val != PetscScalar(0))
                  vector_entry(I, -val);
              }
            }
          }

          // Preassembled linear forms (with block offsets)
          for (auto& lf : pb.getLFs())
          {
            const auto vUUID = lf.getTestFunction().getUUID();
            const size_t vBlock = findTestBlock(vUUID);
            const size_t vOff   = testOffsets[vBlock];

            const auto& vec = lf.getVector();
            PetscInt lo, hi;
            ierr = VecGetOwnershipRange(vec, &lo, &hi);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;

            const PetscScalar* arr;
            ierr = VecGetArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (PetscInt i = lo; i < hi; ++i)
            {
              const PetscScalar val = arr[i - lo];
              if (val != PetscScalar(0))
                vector_entry(static_cast<Index>(vOff) + static_cast<Index>(i), val);
            }
            ierr = VecRestoreArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }

          // Assemble b
          ierr = VecAssemblyBegin(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyEnd(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        if (doMatrix)
        {
          std::vector<PetscInt> rowsToZero;
          for (const Index gs : constraints.getIdentifiedRows())
            if (ownsTrialGlobal(gs))
              rowsToZero.push_back(static_cast<PetscInt>(gs));

          ierr = MatZeroRows(
              A,
              static_cast<PetscInt>(rowsToZero.size()),
              rowsToZero.empty() ? nullptr : rowsToZero.data(),
              0.0,
              nullptr,
              nullptr);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (!ownsTrialGlobal(gs))
              continue;
            const PetscInt I = static_cast<PetscInt>(gs);
            const PetscScalar one = 1.0;
            ierr = MatSetValue(A, I, I, one, ADD_VALUES);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (const auto& e : constraints.expand(gs))
            {
              const PetscInt J = static_cast<PetscInt>(e.index);
              const PetscScalar v = -e.coefficient;
              ierr = MatSetValue(A, I, J, v, ADD_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
            const PetscScalar rhs = constraints.getIdentificationValue(gs);
            if (doVector)
            {
              ierr = VecSetValue(b, I, rhs, INSERT_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }

          ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          if (doVector)
          {
            ierr = VecAssemblyBegin(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = VecAssemblyEnd(b);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        std::vector<PetscInt> bcIdx;
        std::vector<PetscScalar> bcVals;
        for (Index i = 0; i < constraints.size(); i++)
        {
          if (constraints.isFixed(i) && ownsTrialGlobal(i))
          {
            bcIdx.push_back(static_cast<PetscInt>(i));
            bcVals.push_back(constraints.getFixedValue(i));
          }
        }

        if (mode == AssemblyMode::Full)
        {
          Vec bcVec;
          ierr = VecDuplicate(b, &bcVec);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecZeroEntries(bcVec);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecSetValues(
              bcVec,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.empty() ? nullptr : bcIdx.data(),
              bcVals.empty() ? nullptr : bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyBegin(bcVec);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyEnd(bcVec);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatZeroRowsColumns(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.empty() ? nullptr : bcIdx.data(),
              1.0,
              bcVec,
              b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecDestroy(&bcVec);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }
        else if (mode == AssemblyMode::LHS)
        {
          ierr = MatZeroRows(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.empty() ? nullptr : bcIdx.data(),
              1.0,
              nullptr,
              nullptr);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }
        else
        {
          ierr = VecSetValues(
              b,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.empty() ? nullptr : bcIdx.data(),
              bcVals.empty() ? nullptr : bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyBegin(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyEnd(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        (void) ierr;
      }

    public:
      MPI* copy() const noexcept override
      {
        return new MPI(*this);
      }
  };
}

namespace Rodin::PETSc::Math::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using MPI = Rodin::Assembly::MPI<LinearAlgebraType, Operand>;
}

#endif // RODIN_PETSC_ASSEMBLY_MPI_H
