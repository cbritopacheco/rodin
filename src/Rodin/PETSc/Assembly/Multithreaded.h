/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_OPENMP_PETSC_H
#define RODIN_ASSEMBLY_OPENMP_PETSC_H

#include <omp.h>
#include <petsc.h>
#include <petscerror.h>
#include <type_traits>
#include <vector>
#include <memory>
#include <optional>

#include "Rodin/Assembly/OpenMP.h"

namespace Rodin::Assembly
{
  /**
   * @brief OpenMP assembly for PETSc Vec (linear form)
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
      using VectorType     = ::Vec;
      using LinearFormType = Variational::LinearForm<FES, VectorType>;
      using Parent         = AssemblyBase<VectorType, LinearFormType>;
      using InputType      = typename Parent::InputType;

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
        PetscErrorCode ierr;
        const PetscInt n = PetscInt(input.getFES().getSize());

        ierr = VecSetSizes(res, n, PETSC_DECIDE);
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

#pragma omp parallel num_threads(tc)
          {
            // thread-local accumulator
            std::vector<PetscScalar> local(n, PetscScalar(0));
            auto integrator = std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());

#pragma omp for nowait
            for (PetscInt i = 0; i < cnt; ++i)
            {
              if (seq.filter(i) &&
                 (attrs.empty() || attrs.count(mesh.getAttribute(dim, i))))
              {
                auto it = seq.getIterator(i);
                integrator->setPolytope(*it);
                const auto& dofs = input.getFES().getDOFs(dim, i);
                for (size_t k = 0; k < dofs.size(); ++k)
                {
                  local[dofs[k]] += PetscScalar(integrator->integrate(k));
                }
              }
            }

#pragma omp critical
            {
              for (PetscInt idx = 0; idx < n; ++idx)
              {
                if (local[idx] != PetscScalar(0))
                {
                  ierr = VecSetValue(res, idx, local[idx], ADD_VALUES);
                  assert(ierr == PETSC_SUCCESS);
                }
              }
            }
          } // end parallel
        }

        ierr = VecAssemblyBegin(res);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecAssemblyEnd(res);
        assert(ierr == PETSC_SUCCESS);
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };

  /**
   * @brief OpenMP assembly for PETSc Mat (bilinear form)
   */
  template <class TrialFES, class TestFES>
  class OpenMP<
    ::Mat, Variational::BilinearForm<TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat, Variational::BilinearForm<TrialFES, TestFES, ::Mat>>
  {
    public:
      using DotType       = typename FormLanguage::Dot<
                             typename FormLanguage::Traits<TrialFES>::ScalarType,
                             typename FormLanguage::Traits<TestFES>::ScalarType>::Type;
      static_assert(
        std::is_same_v<DotType, PetscScalar>,
        "FES ScalarTypes must yield PetscScalar for PETSc Mat assembly"
      );
      using OperatorType    = ::Mat;
      using BilinearFormType= Variational::BilinearForm<TrialFES, TestFES, OperatorType>;
      using Parent          = AssemblyBase<OperatorType, BilinearFormType>;
      using InputType       = typename Parent::InputType;

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

      void execute(OperatorType& A, const InputType& input) const override
      {
        PetscErrorCode ierr;
        const PetscInt m = PetscInt(input.getTestFES().getSize());
        const PetscInt n = PetscInt(input.getTrialFES().getSize());

        ierr = MatSetSizes(A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetFromOptions(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatSetUp(A);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatZeroEntries(A);
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

#pragma omp parallel num_threads(tc)
          {
            std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
            local.reserve(cnt);
            auto integrator = std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());

#pragma omp for nowait
            for (PetscInt i = 0; i < cnt; ++i)
            {
              if (seq.filter(i) &&
                 (attrs.empty() || attrs.count(mesh.getAttribute(dim, i))))
              {
                auto it = seq.getIterator(i);
                integrator->setPolytope(*it);
                const auto& rows = input.getTestFES().getDOFs(dim, i);
                const auto& cols = input.getTrialFES().getDOFs(dim, i);
                for (size_t r = 0; r < rows.size(); ++r)
                {
                  for (size_t c = 0; c < cols.size(); ++c)
                  {
                    const PetscScalar v = PetscScalar(integrator->integrate(c, r));
                    if (v != PetscScalar(0))
                    {
                      local.emplace_back(rows[r], cols[c], v);
                    }
                  }
                }
              }
            }

#pragma omp critical
            {
              for (auto& [i,j,v] : local)
              {
                ierr = MatSetValue(A, i, j, v, ADD_VALUES);
                assert(ierr == PETSC_SUCCESS);
              }
            }
          }
        }

        // Global contributions
//         for (auto& bfi : input.getGlobalBFIs())
//         {
//           const auto& tAttrs = bfi.getTestAttributes();
//           const auto& rAttrs = bfi.getTrialAttributes();
//           OpenMPIteration testSeq(mesh, bfi.getTestRegion());
//           const PetscInt dimT = PetscInt(testSeq.getDimension());
//           const PetscInt cntT = PetscInt(testSeq.getCount());
// 
// #pragma omp parallel num_threads(tc)
//           {
//             std::vector<std::tuple<PetscInt,PetscInt,PetscScalar>> local;
//             local.reserve(cntT);
//             auto integrator = std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase<PetscScalar>>(bfi.copy());
// 
// #pragma omp for nowait
//             for (PetscInt ti = 0; ti < cntT; ++ti)
//             {
//               if (testSeq.filter(ti) &&
//                  (tAttrs.empty() || tAttrs.count(mesh.getAttribute(dimT, ti))))
//               {
//                 auto teIt = testSeq.getIterator(ti);
//                 SequentialIteration trialSeq{ mesh, integrator->getTrialRegion() };
//                 for (auto trIt = trialSeq.getIterator(); trIt; ++trIt)
//                 {
//                   if (rAttrs.empty() || rAttrs.count(trIt->getAttribute()))
//                   {
//                     integrator->setPolytope(*trIt, *teIt);
//                     const auto& rows = input.getTestFES().getDOFs(dimT, ti);
//                     const auto& cols = input.getTrialFES().getDOFs(PetscInt(trIt->getDimension()), trIt->getIndex());
//                     for (size_t r = 0; r < rows.size(); ++r)
//                     {
//                       for (size_t c = 0; c < cols.size(); ++c)
//                       {
//                         const PetscScalar v = PetscScalar(integrator->integrate(c, r));
//                         if (v != PetscScalar(0))
//                         {
//                           local.emplace_back(rows[r], cols[c], v);
//                         }
//                       }
//                     }
//                   }
//                 }
//               }
//             }
// 
//             #pragma omp critical
//             {
//               for (auto& [i,j,v] : local)
//               {
//                 ierr = MatSetValue(A, i, j, v, ADD_VALUES);
//                 assert(ierr == PETSC_SUCCESS);
//               }
//             }
//           }
//         }

        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);

        ierr = MatAssemblyEnd(A,   MAT_FINAL_ASSEMBLY);
        assert(ierr == PETSC_SUCCESS);
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };
}

#endif // RODIN_ASSEMBLY_OpenMP_PETSC_H
