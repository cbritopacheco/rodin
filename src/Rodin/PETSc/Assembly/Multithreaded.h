/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_MULTITHREADED_PETSC_H
#define RODIN_ASSEMBLY_MULTITHREADED_PETSC_H

#include <petsc.h>
#include <petscerror.h>
#include <type_traits>

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Threads/Mutex.h"
#include "Rodin/Threads/ThreadPool.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/LinearFormIntegratorBase.h"
#include "Rodin/Variational/LocalBilinearFormIntegratorBase.h"
#include "Rodin/Variational/GlobalBilinearFormIntegratorBase.h"
#include "Rodin/Utility/Overloaded.h"
#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @brief Multithreaded assembly for PETSc Vec (linear form)
   */
  template <class FES>
  class Multithreaded<
    ::Vec&,
    Variational::LinearForm<FES, ::Vec&>
  > final
    : public AssemblyBase<
        ::Vec&,
        Variational::LinearForm<FES, ::Vec&>
      >
  {
    public:
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      static_assert(
        std::is_same_v<ScalarType, PetscScalar>,
        "FES::ScalarType must be PetscScalar for PETSc Vec assembly"
      );
      using VectorType     = ::Vec&;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      using InputType      = typename Parent::InputType;

#ifdef RODIN_MULTITHREADED
      Multithreaded()
        : Multithreaded(Threads::getGlobalThreadPool())
      {}
#else
      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
      {}
#endif

      Multithreaded(std::reference_wrapper<Threads::ThreadPool> pool)
        : m_pool(pool)
      {}

      Multithreaded(size_t threadCount)
        : m_pool(threadCount)
      {
        assert(threadCount > 0);
      }

      void execute(VectorType& res, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const PetscInt n = PetscInt(input.getFES().getSize());

        // Create and zero the global Vec
        ierr = VecCreate(PETSC_COMM_SELF, &res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecSetSizes(res, n, PETSC_DECIDE);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecSetFromOptions(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecSet(res, PetscScalar(0));
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        const auto& mesh = input.getFES().getMesh();

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          Internal::MultithreadedIteration seq(mesh, lfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());
          const auto  loop = [&](PetscInt start, PetscInt end)
          {
            // Thread‐local accumulation vector
            std::vector<PetscScalar> local(n, PetscScalar(0));
            auto integrator = std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());

            for (PetscInt i = start; i < end; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() && !attrs.count(mesh.getAttribute(dim, i))) continue;

              const auto it = seq.getIterator(i);
              integrator->setPolytope(*it);
              const auto& dofs = input.getFES().getDOFs(dim, i);
              for (size_t k = 0; k < dofs.size(); ++k)
                local[dofs[k]] += PetscScalar(integrator->integrate(k));
            }

            // Merge into global Vec under lock
            m_mutex.lock();
            for (PetscInt idx = 0; idx < n; ++idx)
            {
              if (local[idx] != PetscScalar(0))
              {
                ierr = VecSetValue(res, idx, local[idx], ADD_VALUES);
                PetscCallAbort(PETSC_COMM_SELF, ierr);
              }
            }
            m_mutex.unlock();
          };

          // Dispatch to pool
          auto& pool = getThreadPool();
          pool.pushLoop(0, seq.getCount(), loop);
          pool.waitForTasks();
        }

        // Finalize assembly
        ierr = VecAssemblyBegin(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = VecAssemblyEnd(res);
        PetscCallAbort(PETSC_COMM_SELF, ierr);
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;

      const Threads::ThreadPool& getThreadPool() const
      {
        return std::holds_alternative<Threads::ThreadPool>(m_pool)
          ? std::get<Threads::ThreadPool>(m_pool)
          : std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
      }
  };

  /**
   * @brief Multithreaded assembly for PETSc Mat (bilinear form)
   */
  template <class TrialFES, class TestFES>
  class Multithreaded<
    ::Mat&,
    Variational::BilinearForm<TrialFES, TestFES, ::Mat&>
  > final
    : public AssemblyBase<
        ::Mat&,
        Variational::BilinearForm<TrialFES, TestFES, ::Mat&>
      >
  {
    public:
      using DotType = typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      using OperatorType    = ::Mat&;
      using BilinearFormType= Variational::BilinearForm<TrialFES, TestFES, OperatorType>;
      using Parent          = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType       = typename Parent::InputType;

      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly");

#ifdef RODIN_MULTITHREADED
      Multithreaded()
        : Multithreaded(Threads::getGlobalThreadPool())
      {}
#else
      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
      {}
#endif

      Multithreaded(std::reference_wrapper<Threads::ThreadPool> pool)
        : m_pool(pool)
      {}

      Multithreaded(size_t threadCount)
        : m_pool(threadCount)
      {
        assert(threadCount > 0);
      }

      void execute(OperatorType& A, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const PetscInt m = PetscInt(input.getTestFES().getSize());
        const PetscInt n = PetscInt(input.getTrialFES().getSize());

        // Create and zero the global Mat
        ierr = MatCreate(PETSC_COMM_SELF, &A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatSetSizes(A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatSetFromOptions(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatSetUp(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatZeroEntries(A);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        const auto& mesh = input.getTestFES().getMesh();

        // Local contributions
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::MultithreadedIteration seq(mesh, bfi.getRegion());
          const PetscInt dim = PetscInt(seq.getDimension());

          const auto loop = [&](PetscInt start, PetscInt end)
          {
            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
            auto integrator = std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());

            for (PetscInt i = start; i < end; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() && !attrs.count(mesh.getAttribute(dim, i))) continue;

              const auto it = seq.getIterator(i);
              integrator->setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(dim, i);
              const auto& cols = input.getTrialFES().getDOFs(dim, i);

              for (size_t r = 0; r < rows.size(); ++r)
                for (size_t c = 0; c < cols.size(); ++c)
                {
                  const PetscScalar v = PetscScalar(integrator->integrate(c, r));
                  if (v != PetscScalar(0))
                    local.emplace_back(rows[r], cols[c], v);
                }
            }

            m_mutex.lock();
            for (auto& [i,j,v] : local)
            {
              ierr = MatSetValue(A, i, j, v, ADD_VALUES);
              PetscCallAbort(PETSC_COMM_SELF, ierr);
            }
            m_mutex.unlock();
          };

          auto& pool = getThreadPool();
          pool.pushLoop(0, seq.getCount(), loop);
          pool.waitForTasks();
        }

        // Global contributions
        for (auto& bfi : input.getGlobalBFIs())
        { 
          const auto& tAttrs = bfi.getTestAttributes();
          const auto& rAttrs = bfi.getTrialAttributes();
          Internal::MultithreadedIteration testSeq(mesh, bfi.getTestRegion());
          const PetscInt dimT = PetscInt(testSeq.getDimension());

          const auto loop = [&](PetscInt start, PetscInt end)
          {
            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
            auto integrator = std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());

            for (PetscInt ti = start; ti < end; ++ti)
            {
              if (!testSeq.filter(ti))
                continue;
              if (!tAttrs.empty() && !tAttrs.count(mesh.getAttribute(dimT, ti)))
                continue;

              const auto te = testSeq.getIterator(ti);
              Internal::SequentialIteration trialSeq{ mesh, integrator->getTrialRegion() };
              for (auto tr = trialSeq.getIterator(); tr; ++tr)
              {
                if (!rAttrs.empty() && !rAttrs.count(tr->getAttribute()))
                  continue;
                integrator->setPolytope(*tr, *te);
                const auto& rows = input.getTestFES().getDOFs(dimT, ti);
                const auto& cols = input.getTrialFES().getDOFs(PetscInt(tr->getDimension()), tr->getIndex());

                for (size_t r = 0; r < rows.size(); ++r)
                  for (size_t c = 0; c < cols.size(); ++c)
                  {
                    const PetscScalar v = PetscScalar(integrator->integrate(c, r));
                    if (v != PetscScalar(0))
                      local.emplace_back(rows[r], cols[c], v);
                  }
              }
            }

            m_mutex.lock();
            for (auto& [i,j,v] : local)
            {
              ierr = MatSetValue(A, i, j, v, ADD_VALUES);
              PetscCallAbort(PETSC_COMM_SELF, ierr);
            }
            m_mutex.unlock();
          };

          auto& pool = getThreadPool();
          pool.pushLoop(0, testSeq.getCount(), loop);
          pool.waitForTasks();
        }

        // Finalize assembly
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(PETSC_COMM_SELF, ierr);

        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        PetscCallAbort(PETSC_COMM_SELF, ierr);
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;

      const Threads::ThreadPool& getThreadPool() const
      {
        return std::holds_alternative<Threads::ThreadPool>(m_pool)
          ? std::get<Threads::ThreadPool>(m_pool)
          : std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
      }
  };
}

#endif // RODIN_ASSEMBLY_MULTITHREADED_PETSC_H

