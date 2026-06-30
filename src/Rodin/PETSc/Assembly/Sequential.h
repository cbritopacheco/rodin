#ifndef RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H
#define RODIN_PETSC_ASSEMBLY_SEQUENTIAL_H

/**
 * @file
 * @brief Sequential assembly specializations targeting PETSc objects.
 */

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <cassert>
#include <optional>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Assembly/Sequential.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/PETSc/Assembly/MatrixSetup.h"
#include "Rodin/PETSc/Assembly/VectorSetup.h"

#include "Rodin/PETSc/Math/Vector.h"
#include "Rodin/PETSc/Math/Matrix.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Assembly
{
  /**
   * @brief Sequential assembly of a PETSc vector from a linear form.
   *
   * Iterates over mesh polytopes on a single thread, evaluates
   * @ref Rodin::Variational::LinearFormIntegrator instances, and inserts
   * entries into the PETSc vector with `VecSetValue`.
   *
   * @tparam FES Finite element space type.
   */
  // Sequential assembly for PETSc Vec (linear form)
  template <class FES>
  class Sequential<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      /// @brief Scalar type of the DOF coefficients (`PetscScalar`).
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief PETSc vector type (`::Vec`).
      using VectorType = ::Vec;
      /// @brief Linear form type being assembled.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      /// @brief Parent assembly base class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType = typename Parent::InputType;

      static_assert(std::is_same_v<ScalarType, PetscScalar>);

      void execute(VectorType& res, const InputType& input) const override
      {
        assert(res);
        const size_t n = input.getFES().getSize();

        PetscErrorCode ierr;

        ierr = PETSc::Assembly::VectorSetup(res).prepare({
          static_cast<PetscInt>(n),
          static_cast<PetscInt>(n),
          nullptr
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            lfi.setPolytope(*it);
            const auto& dofs = input.getFES().getDOFs(it.getDimension(), it->getIndex());
            for (PetscInt l = 0; l < dofs.size(); ++l)
            {
              const PetscScalar v = PetscScalar(lfi.integrate(l));
              ierr = VecSetValue(res, dofs[l], v, ADD_VALUES);
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
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };


  /**
   * @brief Sequential assembly of a PETSc matrix from a bilinear form.
   *
   * Iterates over mesh polytopes on a single thread, evaluates
   * bilinear form integrators, and inserts entries into the PETSc matrix
   * with `MatSetValues`.
   *
   * @tparam Solution Solution type.
   * @tparam TrialFES Trial finite element space type.
   * @tparam TestFES  Test finite element space type.
   */
  // Sequential assembly for PETSc Mat (bilinear form)
  template <class Solution, class TrialFES, class TestFES>
  class Sequential<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      /// @brief Scalar type resulting from the dot product of trial and test scalars.
      using DotType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      /// @brief PETSc matrix type (`::Mat`).
      using OperatorType = ::Mat;
      /// @brief Bilinear form type being assembled.
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;
      /// @brief Parent assembly base class.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType = typename Parent::InputType;

      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly"
      );

      void execute(OperatorType& res, const InputType& input) const override
      {
        assert(res);
        const size_t m = input.getTestFES().getSize();
        const size_t n = input.getTrialFES().getSize();

        PetscErrorCode ierr;
        ierr = PETSc::Assembly::MatrixSetup(res).prepare({
          static_cast<PetscInt>(m),
          static_cast<PetscInt>(n),
          static_cast<PetscInt>(m),
          static_cast<PetscInt>(n),
          nullptr,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        const auto& mesh = input.getTrialFES().getMesh();
        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            const size_t d = it->getDimension();
            const Index i = it->getIndex();

            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            bfi.setPolytope(*it);
            const auto& rows = input.getTestFES().getDOFs(d, i);
            const auto& cols = input.getTrialFES().getDOFs(d, i);
            for (PetscInt i = 0; i < rows.size(); ++i)
            {
              for (PetscInt j = 0; j < cols.size(); ++j)
              {
                const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                ierr = MatSetValue(res, rows[i], cols[j], v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
                (void) ierr;
              }
            }
          }
        }

        // Global contributions
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }
            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }
              bfi.setPolytope(*trIt, *teIt);
              const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
              for (PetscInt i = 0; i < rows.size(); ++i)
                for (PetscInt j = 0; j < cols.size(); ++j)
                {
                  const PetscScalar v = PetscScalar(bfi.integrate(j, i));
                  ierr = MatSetValue(res, rows[i], cols[j], v, ADD_VALUES);
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
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of a single-variable PETSc problem.
   *
   * Assembles a @ref Rodin::Variational::Problem backed by PETSc objects
   * on a single thread, populating both the system matrix and right-hand
   * side vector.
   *
   * @tparam U Trial function type.
   * @tparam V Test function type.
   */
  // Sequential assembly for single-variable Problem (PETSc)
  template <class U, class V>
  class Sequential<
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
      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U, V>;
      /// @brief Parent assembly base class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;
      /// @brief Input data type for the assembly pipeline.
      using InputType = typename Parent::InputType;

      /// @brief PETSc matrix type (`::Mat`) for the system operator.
      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      /// @brief PETSc vector type (`::Vec`) for the RHS and solution.
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec
      /// @brief Scalar type (`PetscScalar`).
      using ScalarType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::ScalarType;   // PetscScalar
      /// @brief Finite element space type for the trial function.
      using TrialFESType = typename Rodin::FormLanguage::Traits<U>::FESType;
      /// @brief Mesh type for the trial finite element space.
      using TrialMeshType = typename Rodin::FormLanguage::Traits<TrialFESType>::MeshType;
      /// @brief Context type (Local or MPI) for the trial mesh.
      using TrialMeshContextType = typename Rodin::FormLanguage::Traits<TrialMeshType>::ContextType;

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
        static_assert(std::is_same_v<TrialMeshContextType, Rodin::Context::Local>,
          "PETSc sequential assembly should only be used with Local mesh context.");

        const bool doMatrix = mode != AssemblyMode::RHS;
        const bool doVector = mode != AssemblyMode::LHS;

        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();

        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES  = v.getFiniteElementSpace();
        const auto& mesh     = trialFES.getMesh();

        const size_t cols = static_cast<size_t>(trialFES.getSize());
        const size_t rows = static_cast<size_t>(testFES.getSize());

        PetscErrorCode ierr;

        assert(A);
        if (doMatrix)
        {
          ierr = PETSc::Assembly::MatrixSetup(A).prepare({
            static_cast<PetscInt>(rows),
            static_cast<PetscInt>(cols),
            static_cast<PetscInt>(rows),
            static_cast<PetscInt>(cols),
            MATSEQAIJ,
            false
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // Vector setup
        assert(b);
        if (doVector)
        {
          ierr = PETSc::Assembly::VectorSetup(b).prepare({
            static_cast<PetscInt>(rows),
            static_cast<PetscInt>(rows),
            VECSEQ,
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
          static_cast<PetscInt>(cols),
          static_cast<PetscInt>(cols),
          VECSEQ,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ConstraintMap<PetscScalar> constraints(std::max(rows, cols));
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

        // Local BFIs
        using MeshType = std::decay_t<decltype(mesh)>;

        if (doMatrix)
        {
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto& attrs = bfi.getAttributes();
            SequentialIteration<MeshType> seq(mesh, bfi.getRegion());
            for (auto it = seq.getIterator(); it; ++it)
            {
              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              const size_t d = it->getDimension();
              const Index  p = it->getIndex();

              bfi.setPolytope(*it);

              const auto& rowsDOF = testFES.getDOFs(d, p);
              const auto& colsDOF = trialFES.getDOFs(d, p);

              for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
              {
                for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
                {
                  const PetscScalar val = PetscScalar(bfi.integrate(j, i));
                  matrix_entry(rowsDOF[i], colsDOF[j], val);
                }
              }
            }
          }

          // Global BFIs
          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();
            SequentialIteration<MeshType> trialseq(mesh, bfi.getTrialRegion());
            SequentialIteration<MeshType> testseq(mesh, bfi.getTestRegion());

            for (auto teIt = testseq.getIterator(); teIt; ++teIt)
            {
              if (!testAttrs.empty())
              {
                const auto a = teIt->getAttribute();
                if (!a || !testAttrs.count(*a))
                  continue;
              }

              const auto& rowsDOF = testFES.getDOFs(teIt->getDimension(), teIt->getIndex());

              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (!trialAttrs.empty())
                {
                  const auto a = trIt->getAttribute();
                  if (!a || !trialAttrs.count(*a))
                    continue;
                }

                const auto& colsDOF = trialFES.getDOFs(trIt->getDimension(), trIt->getIndex());

                bfi.setPolytope(*trIt, *teIt);

                for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
                {
                  for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
                  {
                    const PetscScalar val = PetscScalar(bfi.integrate(j, i));
                    matrix_entry(rowsDOF[i], colsDOF[j], val);
                  }
                }
              }
            }
          }

          // Preassembled bilinear forms
          for (auto& bf : pb.getBFs())
          {
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
                matrix_entry(static_cast<Index>(i), static_cast<Index>(cols[j]), vals[j]);
              ierr = MatRestoreRow(op, i, &nc, &cols, &vals);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }
          }
        }

        // Linear forms
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
          {
            const auto& attrs = lfi.getAttributes();
            SequentialIteration<MeshType> seq(mesh, lfi.getRegion());
            for (auto it = seq.getIterator(); it; ++it)
            {
              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);
              const auto& dofs = testFES.getDOFs(it->getDimension(), it->getIndex());
              for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
              {
                const PetscScalar val = PetscScalar(lfi.integrate(l));
                vector_entry(dofs[l], -val);
              }
            }
          }

          // Preassembled linear forms
          for (auto& lf : pb.getLFs())
          {
            const auto& vec = lf.getVector();
            PetscInt vecSize;
            ierr = VecGetSize(vec, &vecSize);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            const PetscScalar* arr;
            ierr = VecGetArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            for (PetscInt i = 0; i < vecSize; ++i)
              if (arr[i] != PetscScalar(0))
                vector_entry(static_cast<Index>(i), arr[i]);
            ierr = VecRestoreArrayRead(vec, &arr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
          }
        }

        // Assemble A and b
        if (doMatrix)
        {
          ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        if (doVector)
        {
          ierr = VecAssemblyBegin(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          ierr = VecAssemblyEnd(b);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // Identification reconstruction rows
        if (!constraints.getIdentifiedRows().empty())
        {
          std::vector<PetscInt> rowsToZero;
          rowsToZero.reserve(constraints.getIdentifiedRows().size());
          for (const Index gs : constraints.getIdentifiedRows())
            if (static_cast<size_t>(gs) < rows)
              rowsToZero.push_back(static_cast<PetscInt>(gs));

          if (!rowsToZero.empty())
          {
            ierr = MatZeroRows(
                A,
                static_cast<PetscInt>(rowsToZero.size()),
                rowsToZero.data(),
                0.0,
                nullptr,
                nullptr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;

            for (const Index gs : constraints.getIdentifiedRows())
            {
              if (static_cast<size_t>(gs) >= rows)
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
              ierr = VecSetValue(b, I, rhs, INSERT_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }

            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
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

        // Value Dirichlet BCs
        std::vector<PetscInt> bcIdx;
        std::vector<PetscScalar> bcVals;
        for (Index i = 0; i < static_cast<Index>(rows); i++)
        {
          if (constraints.isFixed(i))
          {
            bcIdx.push_back(static_cast<PetscInt>(i));
            bcVals.push_back(constraints.getFixedValue(i));
          }
        }

        if (!bcIdx.empty())
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
                bcIdx.data(),
                bcVals.data(),
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
                bcIdx.data(),
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
                bcIdx.data(),
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
                bcIdx.data(),
                bcVals.data(),
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
      }

    public:
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of a multi-variable PETSc problem.
   *
   * Assembles a block-structured @ref Rodin::Variational::Problem backed
   * by PETSc objects on a single thread.
   *
   * @tparam U1  First trial/test function type.
   * @tparam U2  Second trial/test function type.
   * @tparam U3  Third trial/test function type.
   * @tparam Us  Additional trial/test function types.
   */
  template <class U1, class U2, class U3, class ... Us>
  class Sequential<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      /// @brief PETSc linear system type.
      using LinearSystemType =
        Rodin::PETSc::Math::LinearSystem;

      /// @brief Multi-field problem type being assembled.
      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent assembly base class.
      using Parent =
        AssemblyBase<LinearSystemType, ProblemType>;

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
        const bool doMatrix = mode != AssemblyMode::RHS;
        const bool doVector = mode != AssemblyMode::LHS;

        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        auto&       us           = input.getTrialFunctions();
        auto&       vs           = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto&       trialUUIDMap = input.getTrialUUIDMap();
        auto&       testUUIDMap  = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        using FirstTrialMeshType =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        using MeshContextType =
          typename Rodin::FormLanguage::Traits<FirstTrialMeshType>::ContextType;

        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::Local>,
          "PETSc sequential assembly supports only Local mesh context.");

        PetscErrorCode ierr;

        // ------------------------
        // Allocate / reset A (SeqAIJ); re-use structure across assemblies.
        // ------------------------
        assert(A);
        if (doMatrix)
        {
          ierr = PETSc::Assembly::MatrixSetup(A).prepare({
            static_cast<PetscInt>(nrows),
            static_cast<PetscInt>(ncols),
            static_cast<PetscInt>(nrows),
            static_cast<PetscInt>(ncols),
            MATSEQAIJ,
            false
          });
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        // ------------------------
        // Allocate / reset b (Seq Vec)
        // ------------------------
        assert(b);
        if (doVector)
        {
          ierr = PETSc::Assembly::VectorSetup(b).prepare({
            static_cast<PetscInt>(nrows),
            static_cast<PetscInt>(nrows),
            VECSEQ,
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
          static_cast<PetscInt>(ncols),
          static_cast<PetscInt>(ncols),
          VECSEQ,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        // ------------------------
        // Helpers
        // ------------------------
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

        const auto withTrialFES = [&](const auto& uuid, auto&& fn)
        {
          const size_t k = findTrialBlock(uuid);
          bool found = false;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
            {
              fn(uref.get().getFiniteElementSpace());
              found = true;
            }
          });
          assert(found);
        };

        const auto withTestFES = [&](const auto& uuid, auto&& fn)
        {
          const size_t k = findTestBlock(uuid);
          bool found = false;
          vs.iapply([&](size_t i, const auto& vref)
          {
            if (i == k)
            {
              fn(vref.get().getFiniteElementSpace());
              found = true;
            }
          });
          assert(found);
        };

        using MeshType0 =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        std::optional<std::reference_wrapper<const MeshType0>> meshRef;
        us.apply([&](const auto& uref)
        {
          using MeshT = std::decay_t<decltype(uref.get().getFiniteElementSpace().getMesh())>;
          static_assert(std::is_same_v<MeshT, MeshType0>,
            "Mixed mesh types are not supported in PETSc multi-field sequential assembly.");
          if (!meshRef)
            meshRef = std::cref(uref.get().getFiniteElementSpace().getMesh());
        });
        assert(meshRef.has_value());
        const MeshType0& mesh = meshRef->get();

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
        // ------------------------
        if (doMatrix)
        {
        // Local BFIs
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());

          withTrialFES(uUUID, [&](const auto& uFES)
          {
            withTestFES(vUUID, [&](const auto& vFES)
            {
              for (auto it = seq.getIterator(); it; ++it)
              {
                if (!attrs.empty())
                {
                  const auto a = it->getAttribute();
                  if (!a || !attrs.count(*a))
                    continue;
                }

                const size_t d = it->getDimension();
                const Index  p = it->getIndex();

                bfi.setPolytope(*it);

                const auto& rows = vFES.getDOFs(d, p);
                const auto& cols = uFES.getDOFs(d, p);

                const PetscInt nr = static_cast<PetscInt>(rows.size());
                const PetscInt nc = static_cast<PetscInt>(cols.size());

                for (PetscInt i = 0; i < nr; ++i)
                {
                  const Index I =
                    static_cast<Index>(vOff + static_cast<size_t>(rows[static_cast<size_t>(i)]));
                  for (PetscInt j = 0; j < nc; ++j)
                  {
                    const Index J =
                      static_cast<Index>(uOff + static_cast<size_t>(cols[static_cast<size_t>(j)]));
                    const PetscScalar val =
                      static_cast<PetscScalar>(bfi.integrate(
                            static_cast<size_t>(j), static_cast<size_t>(i)));
                    matrix_entry(I, J, val);
                  }
                }
              }
            });
          });
        }

        // Global BFIs
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          withTrialFES(uUUID, [&](const auto& uFES)
          {
            withTestFES(vUUID, [&](const auto& vFES)
            {
              for (auto teIt = testseq.getIterator(); teIt; ++teIt)
              {
                if (!testAttrs.empty())
                {
                  const auto a = teIt->getAttribute();
                  if (!a || !testAttrs.count(*a))
                    continue;
                }

                const auto& rows = vFES.getDOFs(teIt.getDimension(), teIt->getIndex());

                for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                {
                  if (!trialAttrs.empty())
                  {
                    const auto a = trIt->getAttribute();
                    if (!a || !trialAttrs.count(*a))
                      continue;
                  }

                  const auto& cols = uFES.getDOFs(trIt.getDimension(), trIt->getIndex());

                  bfi.setPolytope(*trIt, *teIt);

                  for (PetscInt i = 0; i < static_cast<PetscInt>(rows.size()); ++i)
                  {
                    const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows[i]));
                    for (PetscInt j = 0; j < static_cast<PetscInt>(cols.size()); ++j)
                    {
                      const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols[j]));
                      const PetscScalar val = static_cast<PetscScalar>(bfi.integrate(j, i));
                      matrix_entry(I, J, val);
                    }
                  }
                }
              }
            });
          });
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
          PetscInt opRows, opCols;
          ierr = MatGetSize(op, &opRows, &opCols);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          for (PetscInt i = 0; i < opRows; ++i)
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

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        } // doMatrix

        // ------------------------
        // Assemble linear terms into b
        // ------------------------
        if (doVector)
        {
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());

          withTestFES(vUUID, [&](const auto& vFES)
          {
            for (auto it = seq.getIterator(); it; ++it)
            {
              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);

              const auto& dofs = vFES.getDOFs(it.getDimension(), it->getIndex());
              for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
              {
                const Index I = static_cast<Index>(vOff + static_cast<size_t>(dofs[l]));
                const PetscScalar val = static_cast<PetscScalar>(lfi.integrate(l));
                vector_entry(I, -val);
              }
            }
          });
        }

        // Preassembled linear forms (with block offsets)
        for (auto& lf : pb.getLFs())
        {
          const auto vUUID = lf.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& vec = lf.getVector();
          PetscInt vecSize;
          ierr = VecGetSize(vec, &vecSize);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;

          const PetscScalar* arr;
          ierr = VecGetArrayRead(vec, &arr);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
          for (PetscInt i = 0; i < vecSize; ++i)
          {
            if (arr[i] != PetscScalar(0))
            {
              vector_entry(static_cast<Index>(vOff) + static_cast<Index>(i), arr[i]);
            }
          }
          ierr = VecRestoreArrayRead(vec, &arr);
          assert(ierr == PETSC_SUCCESS);
          (void) ierr;
        }

        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        } // doVector

        if (!constraints.getIdentifiedRows().empty())
        {
          std::vector<PetscInt> zeroRowsIdx;
          zeroRowsIdx.reserve(constraints.getIdentifiedRows().size());
          for (const Index gs : constraints.getIdentifiedRows())
            if (static_cast<size_t>(gs) < nrows)
              zeroRowsIdx.push_back(static_cast<PetscInt>(gs));

          if (!zeroRowsIdx.empty())
          {
            ierr = MatZeroRows(
                A,
                static_cast<PetscInt>(zeroRowsIdx.size()),
                zeroRowsIdx.data(),
                0.0,
                nullptr, nullptr);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;

            for (const Index gs : constraints.getIdentifiedRows())
            {
              if (static_cast<size_t>(gs) >= nrows)
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
              ierr = VecSetValue(b, I, rhs, INSERT_VALUES);
              assert(ierr == PETSC_SUCCESS);
              (void) ierr;
            }

            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            assert(ierr == PETSC_SUCCESS);
            (void) ierr;
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
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

        std::vector<PetscInt> bcIdx;
        std::vector<PetscScalar> bcVals;
        for (Index i = 0; i < static_cast<Index>(nrows); i++)
        {
          if (constraints.isFixed(i))
          {
            bcIdx.push_back(static_cast<PetscInt>(i));
            bcVals.push_back(constraints.getFixedValue(i));
          }
        }

        if (!bcIdx.empty())
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
                bcIdx.data(),
                bcVals.data(),
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
                bcIdx.data(),
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
                bcIdx.data(),
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
                bcIdx.data(),
                bcVals.data(),
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
      }

    public:
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using Sequential = Rodin::Assembly::Sequential<LinearAlgebraType, Operand>;
}

#endif // RODIN_ASSEMBLY_SEQUENTIAL_PETSC_H
