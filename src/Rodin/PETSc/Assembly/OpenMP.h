/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_OPENMP_H
#define RODIN_PETSC_ASSEMBLY_OPENMP_H

/**
 * @file
 * @brief OpenMP assembly specializations targeting PETSc objects.
 */

#include <omp.h>
#include <petsc.h>
#include <petscerror.h>
#include <algorithm>
#include <cassert>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include "Rodin/Assembly/OpenMP.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/Assembly/MatrixSetup.h"
#include "Rodin/PETSc/Assembly/VectorSetup.h"

namespace Rodin::Assembly
{
  /**
   * @brief OpenMP-parallel assembly of a PETSc vector from a linear form.
   *
   * Uses thread-local buffers and barriers to assemble linear form
   * contributions into a sequential PETSc vector in parallel.
   *
   * @tparam FES Finite element space type.
   */
  template <class FES>
  class OpenMP<
    ::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar for PETSc Vec assembly"
      );

      using VectorType =
        ::Vec;

      using LinearFormType =
        Variational::LinearForm<FES, VectorType>;

      using Parent =
        AssemblyBase<VectorType, LinearFormType>;

      using InputType =
        typename Parent::InputType;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      /// Set number of OpenMP threads
      OpenMP& setThreadCount(size_t tc) noexcept
      {
        m_threadCount = tc;
        return *this;
      }

      /// Get current thread count or max if not set
      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      void execute(VectorType& res, const InputType& input) const override
      {
        assert(res);
        PetscErrorCode ierr;
        const PetscInt n = PetscInt(input.getFES().getSize());

        ierr = PETSc::Assembly::VectorSetup(res).prepare({
          n,
          n,
          nullptr,
          true
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        const auto& mesh = input.getFES().getMesh();
        const int tc     = static_cast<int>(getThreadCount());

        using VectorEntry = std::pair<PetscInt, PetscScalar>;

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());
          const PetscInt cnt = PetscInt(seq.getCount());

          std::vector<std::vector<VectorEntry>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());

            std::vector<VectorEntry> local;
            local.reserve(static_cast<size_t>(cnt));

#pragma omp for
            for (PetscInt i = 0; i < cnt; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty())
              {
                const auto a = mesh.getAttribute(dim, i);
                if (!a || !attrs.count(*a)) continue;
              }

              auto it = seq.getIterator(i);
              integrator->setPolytope(*it);

              const auto& dofs = input.getFES().getDOFs(dim, i);
              for (size_t k = 0; k < dofs.size(); ++k)
              {
                const PetscScalar value = PetscScalar(integrator->integrate(k));
                if (value != PetscScalar(0))
                  local.emplace_back(static_cast<PetscInt>(dofs[k]), value);
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& chunk : chunks)
              {
                if (chunk.empty())
                  continue;
                std::vector<PetscInt> rows;
                std::vector<PetscScalar> vals;
                rows.reserve(chunk.size());
                vals.reserve(chunk.size());
                for (const auto& [row, value] : chunk)
                {
                  rows.push_back(row);
                  vals.push_back(value);
                }
                PetscErrorCode e = VecSetValues(
                    res,
                    static_cast<PetscInt>(rows.size()),
                    rows.data(),
                    vals.data(),
                    ADD_VALUES);
                assert(e == PETSC_SUCCESS);
                (void) e;
              }
            }
          } // end parallel
        }

        PetscErrorCode ierr2 = VecAssemblyBegin(res);
        assert(ierr2 == PETSC_SUCCESS);
        (void) ierr2;
        ierr2 = VecAssemblyEnd(res);
        assert(ierr2 == PETSC_SUCCESS);
        (void) ierr2;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };

  /**
   * @brief OpenMP-parallel assembly of a PETSc matrix from a bilinear form.
   *
   * Uses thread-local buffers and barriers to assemble bilinear form
   * contributions into a sequential PETSc matrix in parallel.
   *
   * @tparam Solution Solution type.
   * @tparam TrialFES Trial finite element space type.
   * @tparam TestFES  Test finite element space type.
   */
  template <class Solution, class TrialFES, class TestFES>
  class OpenMP<
    ::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      using DotType       = typename FormLanguage::Dot<
                             typename FormLanguage::Traits<TrialFES>::ScalarType,
                             typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly"
      );

      using OperatorType =
        ::Mat;

      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<PetscScalar>;

      using Parent =
        AssemblyBase<OperatorType, BilinearFormType>;

      using InputType =
        typename Parent::InputType;

      OpenMP() = default;
      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}
      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      /// Set number of OpenMP threads
      OpenMP& setThreadCount(size_t tc) noexcept
      {
        m_threadCount = tc;
        return *this;
      }

      /// Get current thread count or max if not set
      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      void execute(OperatorType& res, const InputType& input) const override
      {
        assert(res);
        PetscErrorCode ierr;
        const PetscInt m = input.getTestFES().getSize();
        const PetscInt n = input.getTrialFES().getSize();

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

        const auto& mesh = input.getTestFES().getMesh();
        const int tc = static_cast<int>(getThreadCount());
        using MatrixEntry = std::tuple<PetscInt, PetscInt, PetscScalar>;

        auto flush_matrix_entries =
          [&](std::vector<std::vector<MatrixEntry>>& chunks)
        {
          std::vector<MatrixEntry> entries;
          size_t total = 0;
          for (const auto& chunk : chunks)
            total += chunk.size();
          entries.reserve(total);
          for (auto& chunk : chunks)
          {
            entries.insert(
                entries.end(),
                std::make_move_iterator(chunk.begin()),
                std::make_move_iterator(chunk.end()));
            chunk.clear();
          }

          std::sort(entries.begin(), entries.end(),
              [](const MatrixEntry& a, const MatrixEntry& b)
              {
                if (std::get<0>(a) != std::get<0>(b))
                  return std::get<0>(a) < std::get<0>(b);
                return std::get<1>(a) < std::get<1>(b);
              });

          std::vector<PetscInt> cols;
          std::vector<PetscScalar> vals;
          for (size_t pos = 0; pos < entries.size();)
          {
            const PetscInt row = std::get<0>(entries[pos]);
            cols.clear();
            vals.clear();
            while (pos < entries.size() && std::get<0>(entries[pos]) == row)
            {
              cols.push_back(std::get<1>(entries[pos]));
              vals.push_back(std::get<2>(entries[pos]));
              ++pos;
            }
            PetscErrorCode e = MatSetValues(
                res,
                1,
                &row,
                static_cast<PetscInt>(cols.size()),
                cols.data(),
                vals.data(),
                ADD_VALUES);
            assert(e == PETSC_SUCCESS);
            (void) e;
          }
        };

        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          OpenMPIteration seq(mesh, bfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());
          const PetscInt cnt = PetscInt(seq.getCount());

          std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());
            std::vector<MatrixEntry> local;
            local.reserve(static_cast<size_t>(cnt));

#pragma omp for
            for (PetscInt i = 0; i < cnt; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty())
              {
                const auto a = mesh.getAttribute(dim, i);
                if (!a || !attrs.count(*a)) continue;
              }

              auto it = seq.getIterator(i);
              integrator->setPolytope(*it);

              const auto& rows = input.getTestFES().getDOFs(dim, i);
              const auto& cols = input.getTrialFES().getDOFs(dim, i);
              for (size_t r = 0; r < rows.size(); ++r)
                for (size_t c = 0; c < cols.size(); ++c)
                {
                  const PetscScalar v = Math::conj(integrator->integrate(c, r));
                  local.emplace_back(
                      static_cast<PetscInt>(rows[r]),
                      static_cast<PetscInt>(cols[c]),
                      v);
                }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              flush_matrix_entries(chunks);
            }
          } // end parallel
        }

        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          OpenMPIteration trialseq(mesh, bfi.getTrialRegion());
          OpenMPIteration testseq(mesh, bfi.getTestRegion());

          const PetscInt tdim = PetscInt(testseq.getDimension());
          const PetscInt tcnt = PetscInt(testseq.getCount());

          const PetscInt rdim = PetscInt(trialseq.getDimension());
          const PetscInt rcnt = PetscInt(trialseq.getCount());

          std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());
            std::vector<MatrixEntry> local;
            local.reserve(static_cast<size_t>(tcnt));

#pragma omp for
            for (PetscInt te = 0; te < tcnt; ++te)
            {
              if (!testseq.filter(te)) continue;
              if (!testAttrs.empty())
              {
                const auto a = mesh.getAttribute(tdim, te);
                if (!a || !testAttrs.count(*a)) continue;
              }

              auto teIt = testseq.getIterator(te);
              const auto& rows = input.getTestFES().getDOFs(tdim, te);

              for (PetscInt tr = 0; tr < rcnt; ++tr)
              {
                if (!trialseq.filter(tr)) continue;
                if (!trialAttrs.empty())
                {
                  const auto a = mesh.getAttribute(rdim, tr);
                  if (!a || !trialAttrs.count(*a)) continue;
                }

                auto trIt = trialseq.getIterator(tr);
                integrator->setPolytope(*trIt, *teIt);

                const auto& cols = input.getTrialFES().getDOFs(rdim, tr);
                for (size_t r = 0; r < rows.size(); ++r)
                  for (size_t c = 0; c < cols.size(); ++c)
                  {
                    const PetscScalar v = Math::conj(integrator->integrate(c, r));
                    local.emplace_back(
                        static_cast<PetscInt>(rows[r]),
                        static_cast<PetscInt>(cols[c]),
                        v);
                  }
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              flush_matrix_entries(chunks);
            }
          } // end parallel
        }

        PetscErrorCode ierr2 = MatAssemblyBegin(res, MAT_FINAL_ASSEMBLY);
        assert(ierr2 == PETSC_SUCCESS);
        (void) ierr2;
        ierr2 = MatAssemblyEnd(res, MAT_FINAL_ASSEMBLY);
        assert(ierr2 == PETSC_SUCCESS);
        (void) ierr2;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };

  /**
   * @brief OpenMP-parallel assembly of a single-variable PETSc problem.
   *
   * Assembles a @ref Rodin::Variational::Problem backed by PETSc objects
   * using OpenMP threads, populating both the system matrix and right-hand
   * side vector.
   *
   * @tparam U Trial function type.
   * @tparam V Test function type.
   */
  // OpenMP assembly for single-variable Problem (PETSc)
  template <class U, class V>
  class OpenMP<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U, V>>
  {
    public:
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;
      using ProblemType      = Rodin::Variational::Problem<LinearSystemType, U, V>;
      using Parent           = AssemblyBase<LinearSystemType, ProblemType>;
      using InputType        = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec
      using ScalarType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::ScalarType;   // PetscScalar

      using TrialFESType        = typename Rodin::FormLanguage::Traits<U>::FESType;
      using TrialMeshType       = typename Rodin::FormLanguage::Traits<TrialFESType>::MeshType;
      using TrialMeshContextType= typename Rodin::FormLanguage::Traits<TrialMeshType>::ContextType;

      // If you have these base classes in Variational, keep them; otherwise adapt to your hierarchy.
      using LocalBilinearIntegratorBase  = Variational::LocalBilinearFormIntegratorBase<PetscScalar>;
      using GlobalBilinearIntegratorBase = Variational::GlobalBilinearFormIntegratorBase<PetscScalar>;
      using LinearIntegratorBase         = Variational::LinearFormIntegratorBase<PetscScalar>;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other), m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)), m_threadCount(std::move(other.m_threadCount))
      {}

      OpenMP& setThreadCount(size_t tc) noexcept
      {
        m_threadCount = tc;
        return *this;
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

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
          "PETSc OpenMP assembly (sequential objects) supports only Local mesh context.");

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

        const PetscInt ncols = static_cast<PetscInt>(trialFES.getSize());
        const PetscInt nrows = static_cast<PetscInt>(testFES.getSize());

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
          true,
          false
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        ConstraintMap<PetscScalar> constraints(
            static_cast<size_t>(std::max(nrows, ncols)));
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
              {
                constraints.setFixed(
                    static_cast<Index>(local),
                    static_cast<PetscScalar>(value));
              }
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

        using MatrixEntry = std::tuple<PetscInt, PetscInt, PetscScalar>;
        using VectorEntry = std::pair<PetscInt, PetscScalar>;

        auto add_matrix_entries =
          [&](std::vector<MatrixEntry>& local,
              std::vector<VectorEntry>& localRhs,
              Index row, Index col, PetscScalar val)
        {
          const PetscScalar colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : PetscScalar(0);
          for (const auto& r : constraints.expand(row))
          {
            if (colValue != PetscScalar(0))
              localRhs.emplace_back(
                  static_cast<PetscInt>(r.index),
                  -r.coefficient * val * colValue);
            if (doMatrix)
            {
              for (const auto& c : constraints.expand(col))
                local.emplace_back(
                    static_cast<PetscInt>(r.index),
                    static_cast<PetscInt>(c.index),
                    r.coefficient * val * c.coefficient);
            }
          }
        };

        auto add_vector_entries =
          [&](std::vector<VectorEntry>& local, Index row, PetscScalar val)
        {
          if (val == PetscScalar(0))
            return;
          for (const auto& r : constraints.expand(row))
            local.emplace_back(
                static_cast<PetscInt>(r.index),
                r.coefficient * val);
        };

        auto flush_matrix_entries =
          [&](std::vector<std::vector<MatrixEntry>>& chunks)
        {
          std::vector<MatrixEntry> entries;
          size_t total = 0;
          for (const auto& chunk : chunks)
            total += chunk.size();
          entries.reserve(total);
          for (auto& chunk : chunks)
          {
            entries.insert(
                entries.end(),
                std::make_move_iterator(chunk.begin()),
                std::make_move_iterator(chunk.end()));
            chunk.clear();
          }

          std::sort(entries.begin(), entries.end(),
              [](const MatrixEntry& a, const MatrixEntry& b)
              {
                if (std::get<0>(a) != std::get<0>(b))
                  return std::get<0>(a) < std::get<0>(b);
                return std::get<1>(a) < std::get<1>(b);
              });

          std::vector<PetscInt> cols;
          std::vector<PetscScalar> vals;
          for (size_t pos = 0; pos < entries.size();)
          {
            const PetscInt row = std::get<0>(entries[pos]);
            cols.clear();
            vals.clear();
            while (pos < entries.size() && std::get<0>(entries[pos]) == row)
            {
              cols.push_back(std::get<1>(entries[pos]));
              vals.push_back(std::get<2>(entries[pos]));
              ++pos;
            }
            PetscErrorCode e = MatSetValues(
                A,
                1,
                &row,
                static_cast<PetscInt>(cols.size()),
                cols.data(),
                vals.data(),
                ADD_VALUES);
            assert(e == PETSC_SUCCESS);
            (void) e;
          }
        };

        auto flush_vector_entries =
          [&](std::vector<std::vector<VectorEntry>>& chunks)
        {
          for (auto& chunk : chunks)
          {
            if (chunk.empty())
              continue;
            std::vector<PetscInt> rows;
            std::vector<PetscScalar> vals;
            rows.reserve(chunk.size());
            vals.reserve(chunk.size());
            for (const auto& [row, val] : chunk)
            {
              rows.push_back(row);
              vals.push_back(val);
            }
            PetscErrorCode e = VecSetValues(
                b,
                static_cast<PetscInt>(rows.size()),
                rows.data(),
                vals.data(),
                ADD_VALUES);
            assert(e == PETSC_SUCCESS);
            (void) e;
            chunk.clear();
          }
        };

        const int tc = static_cast<int>(getThreadCount());

        // ------------------------
        // Local BFIs (parallel)
        // ------------------------
        if (doMatrix)
        {
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto& attrs = bfi.getAttributes();
            OpenMPIteration seq(mesh, bfi.getRegion());

            const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
            const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

            std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));
            std::vector<std::vector<VectorEntry>> rhsChunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();

              auto integrator =
                std::unique_ptr<LocalBilinearIntegratorBase>(static_cast<LocalBilinearIntegratorBase*>(bfi.copy()));

              std::vector<MatrixEntry> local;
              local.reserve(static_cast<size_t>(cnt));
              std::vector<VectorEntry> localRhs;

#pragma omp for
              for (PetscInt k = 0; k < cnt; ++k)
              {
                if (!seq.filter(k)) continue;
                if (!attrs.empty())
                {
                  const auto a = mesh.getAttribute(dim, k);
                  if (!a || !attrs.count(*a)) continue;
                }

                auto it = seq.getIterator(k);
                integrator->setPolytope(*it);

                const auto& rowsDOF = testFES.getDOFs(dim, k);
                const auto& colsDOF = trialFES.getDOFs(dim, k);

                for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
                {
                  for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
                  {
                    const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(j, i));
                    add_matrix_entries(local, localRhs, rowsDOF[i], colsDOF[j], val);
                  }
                }
              }

              chunks[static_cast<size_t>(tid)] = std::move(local);
              rhsChunks[static_cast<size_t>(tid)] = std::move(localRhs);

#pragma omp barrier
#pragma omp single
              {
                flush_matrix_entries(chunks);
                if (doVector)
                  flush_vector_entries(rhsChunks);
              }
            } // omp parallel
          }

          // ------------------------
          // Global BFIs (parallel over test entities; inner trial loop serial per test entity)
          // ------------------------
          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();

            OpenMPIteration trialseq(mesh, bfi.getTrialRegion());
            OpenMPIteration testseq(mesh,  bfi.getTestRegion());

            const PetscInt tdim = static_cast<PetscInt>(testseq.getDimension());
            const PetscInt tcnt = static_cast<PetscInt>(testseq.getCount());

            const PetscInt rdim = static_cast<PetscInt>(trialseq.getDimension());
            const PetscInt rcnt = static_cast<PetscInt>(trialseq.getCount());

            std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));
            std::vector<std::vector<VectorEntry>> rhsChunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();

              auto integrator =
                std::unique_ptr<GlobalBilinearIntegratorBase>(static_cast<GlobalBilinearIntegratorBase*>(bfi.copy()));

              std::vector<MatrixEntry> local;
              local.reserve(static_cast<size_t>(tcnt));
              std::vector<VectorEntry> localRhs;

#pragma omp for
              for (PetscInt te = 0; te < tcnt; ++te)
              {
                if (!testseq.filter(te)) continue;
                if (!testAttrs.empty())
                {
                  const auto a = mesh.getAttribute(tdim, te);
                  if (!a || !testAttrs.count(*a)) continue;
                }

                auto teIt = testseq.getIterator(te);
                const auto& rowsDOF = testFES.getDOFs(tdim, te);

                for (PetscInt tr = 0; tr < rcnt; ++tr)
                {
                  if (!trialseq.filter(tr)) continue;
                  if (!trialAttrs.empty())
                  {
                    const auto a = mesh.getAttribute(rdim, tr);
                    if (!a || !trialAttrs.count(*a)) continue;
                  }

                  auto trIt = trialseq.getIterator(tr);
                  const auto& colsDOF = trialFES.getDOFs(rdim, tr);

                  integrator->setPolytope(*trIt, *teIt);

                  for (PetscInt i = 0; i < static_cast<PetscInt>(rowsDOF.size()); ++i)
                  {
                    for (PetscInt j = 0; j < static_cast<PetscInt>(colsDOF.size()); ++j)
                    {
                      const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(j, i));
                      add_matrix_entries(local, localRhs, rowsDOF[i], colsDOF[j], val);
                    }
                  }
                }
              }

              chunks[static_cast<size_t>(tid)] = std::move(local);
              rhsChunks[static_cast<size_t>(tid)] = std::move(localRhs);

#pragma omp barrier
#pragma omp single
              {
                flush_matrix_entries(chunks);
                if (doVector)
                  flush_vector_entries(rhsChunks);
              }
            } // omp parallel
          }

          // Preassembled bilinear forms (serial)
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
              {
                std::vector<MatrixEntry> local;
                std::vector<VectorEntry> localRhs;
                add_matrix_entries(
                    local,
                    localRhs,
                    static_cast<Index>(i),
                    static_cast<Index>(cols[j]),
                    vals[j]);
                std::vector<std::vector<MatrixEntry>> matrixChunks(1);
                matrixChunks[0] = std::move(local);
                flush_matrix_entries(matrixChunks);
                if (doVector)
                {
                  std::vector<std::vector<VectorEntry>> vectorChunks(1);
                  vectorChunks[0] = std::move(localRhs);
                  flush_vector_entries(vectorChunks);
                }
              }
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
        // Linear forms (parallel)
        // Note: matches your sequential single-variable sign convention: b += -LF
        // ------------------------
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
          {
            const auto& attrs = lfi.getAttributes();
            OpenMPIteration seq(mesh, lfi.getRegion());

            const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
            const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

            std::vector<std::vector<VectorEntry>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();

              auto integrator =
                std::unique_ptr<LinearIntegratorBase>(static_cast<LinearIntegratorBase*>(lfi.copy()));

              std::vector<VectorEntry> local;

#pragma omp for
              for (PetscInt k = 0; k < cnt; ++k)
              {
                if (!seq.filter(k)) continue;
                if (!attrs.empty())
                {
                  const auto a = mesh.getAttribute(dim, k);
                  if (!a || !attrs.count(*a)) continue;
                }

                auto it = seq.getIterator(k);
                integrator->setPolytope(*it);

                const auto& dofs = testFES.getDOFs(dim, k);
                for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
                {
                  const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(l));
                  add_vector_entries(local, dofs[l], -val);
                }
              }

              chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
              {
                flush_vector_entries(chunks);
              }
            } // omp parallel
          }

          // Preassembled linear forms (serial) : b += LF
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
            std::vector<std::vector<VectorEntry>> chunks(1);
            for (PetscInt i = 0; i < vecSize; ++i)
            {
              add_vector_entries(chunks[0], static_cast<Index>(i), arr[i]);
            }
            flush_vector_entries(chunks);
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

        if (!constraints.getIdentifiedRows().empty())
        {
          std::vector<PetscInt> rowsToZero;
          rowsToZero.reserve(constraints.getIdentifiedRows().size());
          for (const Index gs : constraints.getIdentifiedRows())
            if (static_cast<PetscInt>(gs) < nrows)
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
              if (static_cast<PetscInt>(gs) >= nrows)
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
        }

        std::vector<PetscInt>    bcIdx;
        std::vector<PetscScalar> bcVals;
        for (PetscInt i = 0; i < nrows; i++)
        {
          if (constraints.isFixed(static_cast<Index>(i)))
          {
            bcIdx.push_back(i);
            bcVals.push_back(constraints.getFixedValue(static_cast<Index>(i)));
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
      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };


  /**
   * @brief OpenMP-parallel assembly of a multi-variable PETSc problem.
   *
   * Assembles a block-structured @ref Rodin::Variational::Problem backed
   * by PETSc objects using OpenMP threads.
   *
   * @tparam U1  First trial/test function type.
   * @tparam U2  Second trial/test function type.
   * @tparam U3  Third trial/test function type.
   * @tparam Us  Additional trial/test function types.
   */
  // OpenMP assembly for multi-variable Problem (PETSc)
  template <class U1, class U2, class U3, class ... Us>
  class OpenMP<
      Rodin::PETSc::Math::LinearSystem,
      Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        Rodin::PETSc::Math::LinearSystem,
        Rodin::Variational::Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      using ProblemType =
        Rodin::Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      using Parent    = AssemblyBase<LinearSystemType, ProblemType>;
      using InputType = typename Parent::InputType;

      using OperatorType = typename Rodin::FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename Rodin::FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec

      using LocalBilinearIntegratorBase  = Variational::LocalBilinearFormIntegratorBase<PetscScalar>;
      using GlobalBilinearIntegratorBase = Variational::GlobalBilinearFormIntegratorBase<PetscScalar>;
      using LinearIntegratorBase         = Variational::LinearFormIntegratorBase<PetscScalar>;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other), m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)), m_threadCount(std::move(other.m_threadCount))
      {}

      OpenMP& setThreadCount(size_t tc) noexcept
      {
        m_threadCount = tc;
        return *this;
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

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

        auto& pb            = input.getProblemBody();
        auto& us            = input.getTrialFunctions();
        auto& vs            = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto& trialUUIDMap = input.getTrialUUIDMap();
        auto& testUUIDMap  = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        using FirstTrialMeshType =
          std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>;
        using MeshContextType =
          typename Rodin::FormLanguage::Traits<FirstTrialMeshType>::ContextType;

        static_assert(
          std::is_same_v<MeshContextType, Rodin::Context::Local>,
          "PETSc OpenMP assembly (sequential objects) supports only Local mesh context.");

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
          true,
          false
        });
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;

        // ------------------------
        // Helpers (same as your sequential)
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
            "Mixed mesh types are not supported in PETSc multi-field OpenMP assembly.");
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

        using MatrixEntry = std::tuple<PetscInt, PetscInt, PetscScalar>;
        using VectorEntry = std::pair<PetscInt, PetscScalar>;

        auto add_matrix_entries =
          [&](std::vector<MatrixEntry>& local,
              std::vector<VectorEntry>& localRhs,
              Index row, Index col, PetscScalar val)
        {
          const PetscScalar colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : PetscScalar(0);
          for (const auto& r : constraints.expand(row))
          {
            if (colValue != PetscScalar(0))
              localRhs.emplace_back(
                  static_cast<PetscInt>(r.index),
                  -r.coefficient * val * colValue);
            for (const auto& c : constraints.expand(col))
              local.emplace_back(
                  static_cast<PetscInt>(r.index),
                  static_cast<PetscInt>(c.index),
                  r.coefficient * val * c.coefficient);
          }
        };

        auto add_vector_entries =
          [&](std::vector<VectorEntry>& local, Index row, PetscScalar val)
        {
          if (val == PetscScalar(0))
            return;
          for (const auto& r : constraints.expand(row))
            local.emplace_back(
                static_cast<PetscInt>(r.index),
                r.coefficient * val);
        };

        auto flush_matrix_entries =
          [&](std::vector<std::vector<MatrixEntry>>& chunks)
        {
          std::vector<MatrixEntry> entries;
          size_t total = 0;
          for (const auto& chunk : chunks)
            total += chunk.size();
          entries.reserve(total);
          for (auto& chunk : chunks)
          {
            entries.insert(
                entries.end(),
                std::make_move_iterator(chunk.begin()),
                std::make_move_iterator(chunk.end()));
            chunk.clear();
          }

          std::sort(entries.begin(), entries.end(),
              [](const MatrixEntry& a, const MatrixEntry& b)
              {
                if (std::get<0>(a) != std::get<0>(b))
                  return std::get<0>(a) < std::get<0>(b);
                return std::get<1>(a) < std::get<1>(b);
              });

          std::vector<PetscInt> cols;
          std::vector<PetscScalar> vals;
          for (size_t pos = 0; pos < entries.size();)
          {
            const PetscInt row = std::get<0>(entries[pos]);
            cols.clear();
            vals.clear();
            while (pos < entries.size() && std::get<0>(entries[pos]) == row)
            {
              cols.push_back(std::get<1>(entries[pos]));
              vals.push_back(std::get<2>(entries[pos]));
              ++pos;
            }
            PetscErrorCode e = MatSetValues(
                A,
                1,
                &row,
                static_cast<PetscInt>(cols.size()),
                cols.data(),
                vals.data(),
                ADD_VALUES);
            assert(e == PETSC_SUCCESS);
            (void) e;
          }
        };

        auto flush_vector_entries =
          [&](std::vector<std::vector<VectorEntry>>& chunks)
        {
          for (auto& chunk : chunks)
          {
            if (chunk.empty())
              continue;
            std::vector<PetscInt> rows;
            std::vector<PetscScalar> vals;
            rows.reserve(chunk.size());
            vals.reserve(chunk.size());
            for (const auto& [row, val] : chunk)
            {
              rows.push_back(row);
              vals.push_back(val);
            }
            PetscErrorCode e = VecSetValues(
                b,
                static_cast<PetscInt>(rows.size()),
                rows.data(),
                vals.data(),
                ADD_VALUES);
            assert(e == PETSC_SUCCESS);
            (void) e;
            chunk.clear();
          }
        };

        const int tc = static_cast<int>(getThreadCount());

        // ------------------------
        // Assemble bilinear terms into A (parallel)
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

            const auto& attrs = bfi.getAttributes();
            OpenMPIteration seq(mesh, bfi.getRegion());

            const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
            const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

          std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));
          std::vector<std::vector<VectorEntry>> rhsChunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LocalBilinearIntegratorBase>(static_cast<LocalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<MatrixEntry> local;
            local.reserve(static_cast<size_t>(cnt));
            std::vector<VectorEntry> localRhs;

#pragma omp for
            for (PetscInt k = 0; k < cnt; ++k)
            {
              if (!seq.filter(k)) continue;
              if (!attrs.empty())
              {
                const auto a = mesh.getAttribute(dim, k);
                if (!a || !attrs.count(*a)) continue;
              }

              auto it = seq.getIterator(k);
              withTrialFES(uUUID, [&](const auto& uFES)
              {
                withTestFES(vUUID, [&](const auto& vFES)
                {
                  integrator->setPolytope(*it);
                  const auto& rows = vFES.getDOFs(dim, k);
                  const auto& cols = uFES.getDOFs(dim, k);
                  for (PetscInt i = 0; i < static_cast<PetscInt>(rows.size()); ++i)
                  {
                    const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(rows[i]));
                    for (PetscInt j = 0; j < static_cast<PetscInt>(cols.size()); ++j)
                    {
                      const PetscInt J = static_cast<PetscInt>(uOff + static_cast<size_t>(cols[j]));
                      const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(j, i));
                      add_matrix_entries(local, localRhs, I, J, val);
                    }
                  }
                });
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);
            rhsChunks[static_cast<size_t>(tid)] = std::move(localRhs);

#pragma omp barrier
#pragma omp single
            {
              flush_matrix_entries(chunks);
              flush_vector_entries(rhsChunks);
            }
          } // omp parallel
        }

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

          OpenMPIteration trialseq(mesh, bfi.getTrialRegion());
          OpenMPIteration testseq(mesh,  bfi.getTestRegion());

          const PetscInt tdim = static_cast<PetscInt>(testseq.getDimension());
          const PetscInt tcnt = static_cast<PetscInt>(testseq.getCount());

          const PetscInt rdim = static_cast<PetscInt>(trialseq.getDimension());
          const PetscInt rcnt = static_cast<PetscInt>(trialseq.getCount());

          std::vector<std::vector<MatrixEntry>> chunks(static_cast<size_t>(tc));
          std::vector<std::vector<VectorEntry>> rhsChunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<GlobalBilinearIntegratorBase>(static_cast<GlobalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<MatrixEntry> local;
            local.reserve(static_cast<size_t>(tcnt));
            std::vector<VectorEntry> localRhs;

#pragma omp for
            for (PetscInt te = 0; te < tcnt; ++te)
            {
              if (!testseq.filter(te)) continue;
              if (!testAttrs.empty())
              {
                const auto a = mesh.getAttribute(tdim, te);
                if (!a || !testAttrs.count(*a)) continue;
              }

              auto teIt = testseq.getIterator(te);
              withTrialFES(uUUID, [&](const auto& uFES)
              {
                withTestFES(vUUID, [&](const auto& vFES)
                {
                  const auto& rows = vFES.getDOFs(tdim, te);

                  for (PetscInt tr = 0; tr < rcnt; ++tr)
                  {
                    if (!trialseq.filter(tr)) continue;
                    if (!trialAttrs.empty())
                    {
                      const auto a = mesh.getAttribute(rdim, tr);
                      if (!a || !trialAttrs.count(*a)) continue;
                    }

                    auto trIt = trialseq.getIterator(tr);
                    const auto& cols = uFES.getDOFs(rdim, tr);

                    integrator->setPolytope(*trIt, *teIt);

                    for (PetscInt i = 0; i < static_cast<PetscInt>(rows.size()); ++i)
                    {
                      const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(rows[i]));
                      for (PetscInt j = 0; j < static_cast<PetscInt>(cols.size()); ++j)
                      {
                        const PetscInt J = static_cast<PetscInt>(uOff + static_cast<size_t>(cols[j]));
                        const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(j, i));
                        add_matrix_entries(local, localRhs, I, J, val);
                      }
                    }
                  }
                });
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);
            rhsChunks[static_cast<size_t>(tid)] = std::move(localRhs);

#pragma omp barrier
#pragma omp single
            {
              flush_matrix_entries(chunks);
              flush_vector_entries(rhsChunks);
            }
          } // omp parallel
        }

        // Preassembled bilinear forms (serial, with block offsets)
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
              std::vector<MatrixEntry> local;
              std::vector<VectorEntry> localRhs;
              add_matrix_entries(
                  local,
                  localRhs,
                  static_cast<Index>(vOff) + static_cast<Index>(i),
                  static_cast<Index>(uOff) + static_cast<Index>(cols[j]),
                  vals[j]);
              std::vector<std::vector<MatrixEntry>> matrixChunks(1);
              matrixChunks[0] = std::move(local);
              flush_matrix_entries(matrixChunks);
              std::vector<std::vector<VectorEntry>> vectorChunks(1);
              vectorChunks[0] = std::move(localRhs);
              flush_vector_entries(vectorChunks);
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
        } // doMatrix

        // ------------------------
        // Assemble linear terms into b (parallel)
        // Convention: negate each LFI value, consistent with all other assembly paths.
        // ------------------------
        if (doVector)
        {
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());

          const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
          const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

          std::vector<std::vector<VectorEntry>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LinearIntegratorBase>(static_cast<LinearIntegratorBase*>(lfi.copy()));

            std::vector<VectorEntry> local;

#pragma omp for
            for (PetscInt k = 0; k < cnt; ++k)
            {
              if (!seq.filter(k)) continue;
              if (!attrs.empty())
              {
                const auto a = mesh.getAttribute(dim, k);
                if (!a || !attrs.count(*a)) continue;
              }

              auto it = seq.getIterator(k);
              withTestFES(vUUID, [&](const auto& vFES)
              {
                integrator->setPolytope(*it);
                const auto& dofs = vFES.getDOFs(dim, k);
                for (PetscInt l = 0; l < static_cast<PetscInt>(dofs.size()); ++l)
                {
                  const PetscInt I = static_cast<PetscInt>(vOff + static_cast<size_t>(dofs[l]));
                  const PetscScalar val = static_cast<PetscScalar>(integrator->integrate(l));
                  add_vector_entries(local, I, -val);
                }
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              flush_vector_entries(chunks);
            }
          } // omp parallel
        }

        // Preassembled linear forms (serial, with block offsets)
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
          std::vector<std::vector<VectorEntry>> chunks(1);
          for (PetscInt i = 0; i < vecSize; ++i)
          {
            if (arr[i] != PetscScalar(0))
            {
              add_vector_entries(
                  chunks[0],
                  static_cast<Index>(vOff) + static_cast<Index>(i),
                  arr[i]);
            }
          }
          flush_vector_entries(chunks);
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
        } // doVector

        if (!constraints.getIdentifiedRows().empty())
        {
          std::vector<PetscInt> zeroRowsIdx;
          zeroRowsIdx.reserve(constraints.getIdentifiedRows().size());
          for (const Index gs : constraints.getIdentifiedRows())
            if (gs < nrows)
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
              if (gs >= nrows)
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

        std::vector<PetscInt>    bcIdx;
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
      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };
}

namespace Rodin::PETSc::Assembly
{
  template <class LinearAlgebraType, class Operand>
  using OpenMP = Rodin::Assembly::OpenMP<LinearAlgebraType, Operand>;
}

#endif // RODIN_PETSC_ASSEMBLY_OPENMP_H
