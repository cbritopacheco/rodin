/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_OPENMP_PETSC_H
#define RODIN_ASSEMBLY_OPENMP_PETSC_H

/**
 * @file
 * @brief OpenMP assembly specializations targeting PETSc objects.
 */

#include <omp.h>
#include <petsc.h>
#include <petscerror.h>
#include <optional>

#include "Rodin/Assembly/OpenMP.h"
#include "Rodin/PETSc/Math/LinearSystem.h"

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

        ierr = VecSetSizes(res, n, n);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSet(res, PetscScalar(0));
        assert(ierr == PETSC_SUCCESS);

        const auto& mesh = input.getFES().getMesh();
        const int tc     = static_cast<int>(getThreadCount());

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());
          const PetscInt cnt = PetscInt(seq.getCount());

          // one dense buffer per thread, merged once
          std::vector<std::vector<PetscScalar>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());

            std::vector<PetscScalar> local(static_cast<size_t>(n), PetscScalar(0));

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
                local[dofs[k]] += PetscScalar(integrator->integrate(k));
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              // single serial flush to PETSc, avoids long critical sections
              for (auto& v : chunks)
              {
                for (PetscInt idx = 0; idx < n; ++idx)
                {
                  const PetscScalar val = v[static_cast<size_t>(idx)];
                  if (val != PetscScalar(0))
                  {
                    PetscErrorCode e = VecSetValue(res, idx, val, ADD_VALUES);
                    assert(e == PETSC_SUCCESS);
                  }
                }
              }
            }
          } // end parallel
        }

        PetscErrorCode ierr2 = VecAssemblyBegin(res);
        assert(ierr2 == PETSC_SUCCESS);
        ierr2 = VecAssemblyEnd(res);
        assert(ierr2 == PETSC_SUCCESS);
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

        ierr = MatSetSizes(res, m, n, m, n);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetFromOptions(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(res);
        assert(ierr == PETSC_SUCCESS);

        const auto& mesh = input.getTestFES().getMesh();
        const int tc = static_cast<int>(getThreadCount());

        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          OpenMPIteration seq(mesh, bfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());
          const PetscInt cnt = PetscInt(seq.getCount());

          // per-thread triplet chunks, merged once
          std::vector<std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());
            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
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
                  if (v != PetscScalar(0))
                    local.emplace_back(rows[r], cols[c], v);
                }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              // single serial flush to PETSc
              for (auto& buf : chunks)
              {
                for (auto& [i,j,v] : buf)
                {
                  PetscErrorCode e = MatSetValue(res, i, j, v, ADD_VALUES);
                  assert(e == PETSC_SUCCESS);
                }
              }
            }
          } // end parallel
        }

        PetscErrorCode ierr2 = MatAssemblyBegin(res, MAT_FINAL_ASSEMBLY);
        assert(ierr2 == PETSC_SUCCESS);
        ierr2 = MatAssemblyEnd(res, MAT_FINAL_ASSEMBLY);
        assert(ierr2 == PETSC_SUCCESS);
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
        static_assert(std::is_same_v<TrialMeshContextType, Rodin::Context::Local>,
          "PETSc OpenMP assembly (sequential objects) supports only Local mesh context.");

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
        // Allocate / reset A (SeqAIJ)
        // ------------------------
        assert(A);
        ierr = MatSetSizes(A, nrows, ncols, nrows, ncols);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetType(A, MATSEQAIJ);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Allocate / reset b (Seq Vec)
        // ------------------------
        assert(b);
        ierr = VecSetSizes(b, nrows, nrows);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetType(b, VECSEQ);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(x, ncols, ncols);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetType(x, VECSEQ);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

        const int tc = static_cast<int>(getThreadCount());

        // ------------------------
        // Local BFIs (parallel)
        // ------------------------
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          OpenMPIteration seq(mesh, bfi.getRegion());

          const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
          const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

          std::vector<std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LocalBilinearIntegratorBase>(static_cast<LocalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
            local.reserve(static_cast<size_t>(cnt));

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
                  if (val != PetscScalar(0))
                    local.emplace_back(rowsDOF[i], colsDOF[j], val);
                }
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& buf : chunks)
              {
                for (auto& [I,J,val] : buf)
                {
                  PetscErrorCode e = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(e == PETSC_SUCCESS);
                }
              }
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

          std::vector<std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<GlobalBilinearIntegratorBase>(static_cast<GlobalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
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
                    if (val != PetscScalar(0))
                      local.emplace_back(rowsDOF[i], colsDOF[j], val);
                  }
                }
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& buf : chunks)
              {
                for (auto& [I,J,val] : buf)
                {
                  PetscErrorCode e = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(e == PETSC_SUCCESS);
                }
              }
            }
          } // omp parallel
        }

        // Preassembled bilinear forms (serial)
        for (auto& bf : pb.getBFs())
        {
          const auto& op = bf.getOperator();
          ierr = MatAXPY(A, 1.0, op, DIFFERENT_NONZERO_PATTERN);
          assert(ierr == PETSC_SUCCESS);
        }

        // Assemble A
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Linear forms (parallel)
        // Note: matches your sequential single-variable sign convention: b += -LF
        // ------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());

          const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
          const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

          std::vector<std::vector<PetscScalar>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LinearIntegratorBase>(static_cast<LinearIntegratorBase*>(lfi.copy()));

            std::vector<PetscScalar> local(static_cast<size_t>(nrows), PetscScalar(0));

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
                local[static_cast<size_t>(dofs[l])] += (-val);
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& vecLocal : chunks)
              {
                for (PetscInt i = 0; i < nrows; ++i)
                {
                  const PetscScalar val = vecLocal[static_cast<size_t>(i)];
                  if (val != PetscScalar(0))
                  {
                    PetscErrorCode e = VecSetValue(b, i, val, ADD_VALUES);
                    assert(e == PETSC_SUCCESS);
                  }
                }
              }
            }
          } // omp parallel
        }

        // Preassembled linear forms (serial) : b += -LF
        for (auto& lf : pb.getLFs())
        {
          ierr = VecAXPY(b, -1.0, lf.getVector());
          assert(ierr == PETSC_SUCCESS);
        }

        // Assemble b
        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Dirichlet BCs (serial)
        // ------------------------
        std::vector<PetscInt>    bcIdx;
        std::vector<PetscScalar> bcVals;

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
            bcIdx.push_back(static_cast<PetscInt>(local));
            bcVals.push_back(static_cast<PetscScalar>(value));
          }
        }

        if (!bcIdx.empty())
        {
          Vec bcVec;
          ierr = VecDuplicate(b, &bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetValues(
              bcVec,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecAssemblyBegin(bcVec);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecAssemblyEnd(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = MatZeroRowsColumns(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              1.0,
              bcVec,
              b);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecDestroy(&bcVec);
          assert(ierr == PETSC_SUCCESS);
        }
      }

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
        // Allocate / reset A (SeqAIJ)
        // ------------------------
        assert(A);
        ierr = MatSetSizes(A,
                           nrows,
                           ncols,
                           nrows,
                           ncols);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetType(A, MATSEQAIJ);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Allocate / reset b (Seq Vec)
        // ------------------------
        assert(b);
        ierr = VecSetSizes(b, nrows, nrows);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetType(b, VECSEQ);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecSetFromOptions(b);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecZeroEntries(b);
        assert(ierr == PETSC_SUCCESS);

        auto& x = axb.getSolution();
        assert(x);
        ierr = VecSetSizes(x, ncols, ncols);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetType(x, VECSEQ);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecSetFromOptions(x);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecZeroEntries(x);
        assert(ierr == PETSC_SUCCESS);

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

        const int tc = static_cast<int>(getThreadCount());

        // ------------------------
        // Assemble bilinear terms into A (parallel)
        // ------------------------
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

          std::vector<std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LocalBilinearIntegratorBase>(static_cast<LocalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
            local.reserve(static_cast<size_t>(cnt));

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
                      if (val != PetscScalar(0))
                        local.emplace_back(I, J, val);
                    }
                  }
                });
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& buf : chunks)
              {
                for (auto& [I,J,val] : buf)
                {
                  PetscErrorCode e = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(e == PETSC_SUCCESS);
                }
              }
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

          std::vector<std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<GlobalBilinearIntegratorBase>(static_cast<GlobalBilinearIntegratorBase*>(bfi.copy()));

            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
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
                        if (val != PetscScalar(0))
                          local.emplace_back(I, J, val);
                      }
                    }
                  }
                });
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& buf : chunks)
              {
                for (auto& [I,J,val] : buf)
                {
                  PetscErrorCode e = MatSetValue(A, I, J, val, ADD_VALUES);
                  assert(e == PETSC_SUCCESS);
                }
              }
            }
          } // omp parallel
        }

        // Assemble A
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Assemble linear terms into b (parallel)
        // Note: matches your sequential multi-variable sign convention (no minus)
        // ------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());

          const PetscInt dim = static_cast<PetscInt>(seq.getDimension());
          const PetscInt cnt = static_cast<PetscInt>(seq.getCount());

          std::vector<std::vector<PetscScalar>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();

            auto integrator =
              std::unique_ptr<LinearIntegratorBase>(static_cast<LinearIntegratorBase*>(lfi.copy()));

            std::vector<PetscScalar> local(static_cast<size_t>(nrows), PetscScalar(0));

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
                  local[static_cast<size_t>(I)] += (-val);
                }
              });
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& vecLocal : chunks)
              {
                for (PetscInt i = 0; i < static_cast<PetscInt>(nrows); ++i)
                {
                  const PetscScalar val = vecLocal[static_cast<size_t>(i)];
                  if (val != PetscScalar(0))
                  {
                    PetscErrorCode e = VecSetValue(b, i, val, ADD_VALUES);
                    assert(e == PETSC_SUCCESS);
                  }
                }
              }
            }
          } // omp parallel
        }

        // Assemble b
        ierr = VecAssemblyBegin(b);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecAssemblyEnd(b);
        assert(ierr == PETSC_SUCCESS);

        // ------------------------
        // Impose Dirichlet BCs (serial)
        // ------------------------
        std::vector<PetscInt>    bcIdx;
        std::vector<PetscScalar> bcVals;

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
            bcIdx.push_back(static_cast<PetscInt>(uOff + local));
            bcVals.push_back(static_cast<PetscScalar>(value));
          }
        }

        if (!bcIdx.empty())
        {
          Vec bcVec;
          ierr = VecDuplicate(b, &bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetValues(
              bcVec,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              bcVals.data(),
              INSERT_VALUES);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecAssemblyBegin(bcVec);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecAssemblyEnd(bcVec);
          assert(ierr == PETSC_SUCCESS);

          ierr = MatZeroRowsColumns(
              A,
              static_cast<PetscInt>(bcIdx.size()),
              bcIdx.data(),
              1.0,
              bcVec,
              b);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecDestroy(&bcVec);
          assert(ierr == PETSC_SUCCESS);
        }
      }

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

#endif // RODIN_ASSEMBLY_OPENMP_PETSC_H
