/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_OPENMP_H
#define RODIN_ASSEMBLY_OPENMP_H

#include <omp.h>

#include "Rodin/Math/Vector.h"

#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Rodin/Assembly/AssemblyBase.h"

#include "Sequential.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @brief OpenMP-based parallel mesh iteration for multi-threaded assembly.
   *
   * OpenMPIteration provides a parallel iteration strategy over mesh elements
   * using OpenMP threading for assembly operations. This class manages the
   * distribution of mesh elements across multiple threads while ensuring
   * thread-safe access to mesh data structures.
   *
   * @tparam MeshType Type of mesh to iterate over (specialized for local context)
   *
   * ## Parallelization Strategy
   * The implementation uses OpenMP to parallelize element-wise assembly:
   * - Each thread processes a subset of mesh elements
   * - Local contributions are computed independently
   * - Synchronization occurs only during global assembly
   */
  template <>
  class OpenMPIteration<Geometry::Mesh<Context::Local>>
  {
    public:
      /// @brief Mesh type for local execution context
      using MeshType = Geometry::Mesh<Context::Local>;

      /**
       * @brief Constructs parallel iteration over a mesh region.
       *
       * @param mesh Mesh to iterate over
       * @param region Geometric region defining the iteration domain
       */
      OpenMPIteration(const MeshType& mesh, const Geometry::Region& region);

      /**
       * @brief Gets an iterator for a specific thread.
       *
       * @param i Thread index
       * @return PolytopeIterator for the specified thread's element subset
       */
      Geometry::PolytopeIterator getIterator(Index i) const;

      /**
       * @brief Gets the topological dimension of elements.
       * @return Dimension of the polytopes being iterated
       */
      size_t getDimension() const;

      /**
       * @brief Gets the total number of elements to process.
       * @return Total element count in the iteration domain
       */
      size_t getCount() const;

      /**
       * @brief Checks if an element should be processed by this iteration.
       *
       * @param i Element index
       * @return true if element should be processed, false otherwise
       */
      bool filter(Index i) const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;  ///< Reference to the mesh
      Geometry::Region m_region;                      ///< Region to iterate over
  };

  /// @brief Template argument deduction guide for OpenMPIteration
  OpenMPIteration(const Geometry::Mesh<Context::Local>& mesh, const Geometry::Region&)
    -> OpenMPIteration<Geometry::Mesh<Context::Local>>;

  /**
   * @brief OpenMP-based parallel assembly for bilinear forms.
   *
   * This class provides a multi-threaded assembly implementation for bilinear
   * forms @f$ a(u,v) @f$ using OpenMP parallelization. The assembly process
   * distributes element computations across multiple threads and uses sparse
   * matrix triplets for thread-safe accumulation.
   *
   * @tparam Solution Solution variable type
   * @tparam TrialFES Trial finite element space type  
   * @tparam TestFES Test finite element space type
   *
   * ## Parallel Assembly Algorithm
   * 1. **Element Distribution**: Distribute elements across OpenMP threads
   * 2. **Local Assembly**: Each thread computes local element matrices independently
   * 3. **Triplet Accumulation**: Store contributions as Eigen triplets for thread safety
   * 4. **Global Assembly**: Combine triplets into final sparse matrix
   *
   * This approach minimizes synchronization overhead while ensuring correctness.
   */
  template <class Solution, class TrialFES, class TestFES>
  class OpenMP<
    std::vector<Eigen::Triplet<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
    Variational::BilinearForm<Solution, TrialFES, TestFES,
      std::vector<Eigen::Triplet<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>> final
    : public AssemblyBase<
        std::vector<Eigen::Triplet<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
        Variational::BilinearForm<Solution, TrialFES, TestFES,
          std::vector<Eigen::Triplet<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        const size_t capacity =
          input.getTestFES().getSize() * std::log(input.getTrialFES().getSize());
        res.clear();
        res.reserve(capacity);

        const auto& mesh = input.getTestFES().getMesh();

        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          OpenMPIteration seq(mesh, bfi.getRegion());
          const Index d = seq.getDimension();
          const Index count = seq.getCount();
          const int tc =
            m_threadCount.has_value() ? static_cast<int>(m_threadCount.value()) : omp_get_max_threads();

          std::vector<OperatorType> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            std::unique_ptr<LocalBilinearFormIntegratorBaseType> integrator(bfi.copy());
            OperatorType local;
            local.reserve(capacity / std::max(1, tc));

#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, i))) continue;

              auto it = seq.getIterator(i);
              integrator->setPolytope(*it);

              const auto& rows = input.getTestFES().getDOFs(d, i);
              const auto& cols = input.getTrialFES().getDOFs(d, i);
              for (size_t r = 0; r < rows.size(); ++r)
                for (size_t c = 0; c < cols.size(); ++c)
                {
                  const ScalarType s = Math::conj(integrator->integrate(c, r));
                  if (s != ScalarType(0))
                    local.emplace_back(rows(r), cols(c), s);
                }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& v : chunks)
              {
                res.insert(res.end(),
                           std::make_move_iterator(v.begin()),
                           std::make_move_iterator(v.end()));
              }
            }
          } // end parallel
        }
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      OpenMP& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };

  /**
   * @brief OpenMP assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class Solution, class TrialFES, class TestFES>
  class OpenMP<
    Math::SparseMatrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
      Solution,
      TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>> final
    : public AssemblyBase<
        Math::SparseMatrix<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Variational::BilinearForm<
          Solution,
          TrialFES, TestFES,
          Math::SparseMatrix<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      OpenMP() = default;

      OpenMP(size_t threadCount)
        : m_assembly(threadCount)
      {
        assert(threadCount > 0);
      }

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_assembly(other.m_assembly)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        std::vector<Eigen::Triplet<ScalarType>> triplets;
        m_assembly.execute(triplets, {
            input.getTrialFES(), input.getTestFES(),
            input.getLocalBFIs(), input.getGlobalBFIs() });
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
      }

      size_t getThreadCount() const noexcept
      {
        return m_assembly.getThreadCount();
      }

      OpenMP& setThreadCount(size_t threadCount) noexcept
      {
        m_assembly.setThreadCount(threadCount);
        return *this;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      OpenMP<
        std::vector<Eigen::Triplet<ScalarType>>,
        Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<ScalarType>>>> m_assembly;
  };

  template <class Solution, class TrialFES, class TestFES>
  class OpenMP<
    Math::Matrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
      Solution,
      TrialFES, TestFES,
      Math::Matrix<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>> final
    : public AssemblyBase<
        Math::Matrix<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Variational::BilinearForm<
          Solution,
          TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = Math::Matrix<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      void execute(OperatorType& res, const InputType& input) const override
      {
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setZero();

        const auto& mesh = input.getTestFES().getMesh();
        const int tc =
          m_threadCount.has_value() ? static_cast<int>(m_threadCount.value()) : omp_get_max_threads();

        // --- Local contributions ---
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          OpenMPIteration seq(mesh, bfi.getRegion());
          const Index d     = seq.getDimension();
          const Index count = seq.getCount();

          std::vector<OperatorType> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto lbfi = std::unique_ptr<LocalBilinearFormIntegratorBaseType>(bfi.copy());
            OperatorType local;
            local.resize(res.rows(), res.cols());
            local.setZero();

#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, i))) continue;

              auto it = seq.getIterator(i);
              lbfi->setPolytope(*it);

              const auto& rows = input.getTestFES().getDOFs(d, i);
              const auto& cols = input.getTrialFES().getDOFs(d, i);
              for (size_t r = 0; r < rows.size(); ++r)
                for (size_t c = 0; c < cols.size(); ++c)
                  local(rows(r), cols(c)) += Math::conj(lbfi->integrate(c, r));
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& M : chunks) res += M;
            }
          } // end parallel
        }

        // --- Global (coupling) contributions ---
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& testAttrs  = bfi.getTestAttributes();
          const auto& trialAttrs = bfi.getTrialAttributes();
          OpenMPIteration testseq(mesh, bfi.getTestRegion());
          const Index d     = testseq.getDimension();
          const Index count = testseq.getCount();

          std::vector<OperatorType> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto gbfi = std::unique_ptr<GlobalBilinearFormIntegratorBaseType>(bfi.copy());
            OperatorType local;
            local.resize(res.rows(), res.cols());
            local.setZero();

#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              if (!testseq.filter(i)) continue;
              if (!testAttrs.empty() && !testAttrs.count(mesh.getAttribute(d, i))) continue;

              auto teIt = testseq.getIterator(i);
              SequentialIteration trialseq{ mesh, gbfi->getTrialRegion() };

              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (!trialAttrs.empty() && !trialAttrs.count(trIt->getAttribute())) continue;

                gbfi->setPolytope(*trIt, *teIt);

                const auto& rows = input.getTestFES().getDOFs(d, teIt->getIndex());
                const auto& cols = input.getTrialFES().getDOFs(d, trIt->getIndex());
                for (size_t r = 0; r < rows.size(); ++r)
                  for (size_t c = 0; c < cols.size(); ++c)
                    local(rows(r), cols(c)) += Math::conj(gbfi->integrate(c, r));
              }
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& M : chunks) res += M;
            }
          } // end parallel
        }
      }

      OpenMP& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };

  /**
   * @brief %OpenMP assembly of the Math::Vector associated to a
   * LinearForm object.
   */
  template <class FES>
  class OpenMP<
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
    Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
    : public AssemblyBase<
        Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
        Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using VectorType = Math::Vector<ScalarType>;

      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      using Parent = AssemblyBase<VectorType, LinearFormType>;

      using InputType = typename Parent::InputType;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      void execute(VectorType& res, const InputType& input) const override
      {
        // initialize global vector
        res.resize(input.getFES().getSize());
        res.setZero();

        const auto& mesh = input.getFES().getMesh();
        const int tc =
          m_threadCount.has_value() ? static_cast<int>(m_threadCount.value()) : omp_get_max_threads();

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());
          const Index d = seq.getDimension();
          const Index count = seq.getCount();

          std::vector<VectorType> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
          {
            const int tid = omp_get_thread_num();
            auto integrator =
              std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());
            VectorType local;
            local.resize(res.size());
            local.setZero();

#pragma omp for
            for (Index i = 0; i < count; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, i))) continue;

              auto it = seq.getIterator(i);
              integrator->setPolytope(*it);

              const auto& dofs = input.getFES().getDOFs(d, i);
              for (size_t k = 0; k < dofs.size(); ++k)
                local(dofs(k)) += integrator->integrate(k);
            }

            chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
            {
              for (auto& v : chunks) res += v;
            }
          } // end parallel
        }
      }

      OpenMP& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };

  template <class Scalar, class Solution, class FES, class ValueDerived>
  class OpenMP<
    IndexMap<Scalar>,
    Variational::DirichletBC<
      Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>> final
    : public AssemblyBase<
        IndexMap<Scalar>,
        Variational::DirichletBC<
          Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>>
  {
    public:
      using FESType = FES;

      using TrialFunctionType = Variational::TrialFunction<Solution, FES>;

      using ValueType = Variational::FunctionBase<ValueDerived>;

      using DirichletBCType = Variational::DirichletBC<TrialFunctionType, ValueType>;

      using Parent = AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      using FESRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using InputType = typename Parent::InputType;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      void execute(IndexMap<Scalar>& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        const auto& value = input.getValue();
        const auto& essBdr = input.getEssentialBoundary();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t faceDim = mesh.getDimension() - 1;
        const size_t faceCount = mesh.getFaceCount();
        const int tc =
          m_threadCount.has_value() ? static_cast<int>(m_threadCount.value()) : omp_get_max_threads();

        res.clear();
        std::vector<IndexMap<Scalar>> chunks(static_cast<size_t>(tc));

#pragma omp parallel num_threads(tc)
        {
          const int tid = omp_get_thread_num();
          IndexMap<Scalar> local;

#pragma omp for
          for (Index i = 0; i < faceCount; ++i)
          {
            if (!mesh.isBoundary(i)) continue;
            if (!essBdr.empty() && !essBdr.count(mesh.getAttribute(faceDim, i))) continue;

            const auto& fe = fes.getFiniteElement(faceDim, i);
            const auto& mapping = fes.getPullback({ faceDim, i }, value);
            for (Index k = 0; k < fe.getCount(); ++k)
            {
              const Index g = fes.getGlobalIndex({ faceDim, i }, k);
              auto it = local.find(g);
              if (it == local.end())
                local.insert(it, std::pair{ g, fe.getLinearForm(k)(mapping) });
            }
          }

          chunks[static_cast<size_t>(tid)] = std::move(local);

#pragma omp barrier
#pragma omp single
          {
            // ordered merge guarded by implicit single barrier
            for (auto& partial : chunks)
              res.insert(boost::container::ordered_unique_range, partial.begin(), partial.end());
          }
        } // end parallel
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      OpenMP& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      OpenMP* copy() const noexcept override
      {
        return new OpenMP(*this);
      }

    private:
      Optional<size_t> m_threadCount;
  };
}

#endif
