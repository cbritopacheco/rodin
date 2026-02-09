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
              assert(rows.size() >= 0);
              for (size_t r = 0; r < static_cast<size_t>(rows.size()); ++r)
              {
                assert(cols.size() >= 0);
                for (size_t c = 0; c < static_cast<size_t>(cols.size()); ++c)
                {
                  const ScalarType s = Math::conj(integrator->integrate(c, r));
                  if (s != ScalarType(0))
                    local.emplace_back(rows(r), cols(c), s);
                }
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
              assert(dofs.size() >= 0);
              for (size_t k = 0; k < static_cast<size_t>(dofs.size()); ++k)
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

  template <class LinearSystem, class TrialFunction, class TestFunction>
  class OpenMP<
    LinearSystem,
    Variational::Problem<LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      using LinearSystemType = LinearSystem;

      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, TrialFunction, TestFunction>>;

      using InputType = typename Parent::InputType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseType =
        Variational::LinearFormIntegratorBase<ScalarType>;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)), m_threadCount(other.m_threadCount)
      {}

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
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

        b.resize(rows);
        b.setZero();

        constexpr bool IsSparse =
          std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>;

        const int tc =
          m_threadCount.has_value()
            ? static_cast<int>(m_threadCount.value())
            : omp_get_max_threads();

        // ---- Dirichlet fixed values (NaN sentinel) ----
        const size_t ndofs = std::max(rows, cols);
        std::vector<ScalarType> fixed(ndofs, Math::nan<ScalarType>());

        auto isFixed = [&](Index i) -> bool
        {
          const size_t k = static_cast<size_t>(i);
          return k < fixed.size() && !Math::isNaN(fixed[k]);
        };

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          for (const auto& [local, value] : dbc.getDOFs())
          {
            const Index I = static_cast<Index>(local);
            fixed[static_cast<size_t>(I)] = static_cast<ScalarType>(value);
          }
        }

        if constexpr (IsSparse)
        {
          // ---- Sparse path: eliminate during assembly (thread-local triplets + RHS) ----
          std::vector<std::vector<Eigen::Triplet<ScalarType>>> tchunks(static_cast<size_t>(tc));
          std::vector<std::vector<ScalarType>> rhsChunks(
            static_cast<size_t>(tc),
            std::vector<ScalarType>(rows, ScalarType(0))
          );

          auto sparse_entry =
            [&](std::vector<Eigen::Triplet<ScalarType>>& localT,
                std::vector<ScalarType>& localRhs,
                Index row, Index col, ScalarType val)
          {
            if (val == ScalarType(0))
              return;

            const bool rowFixed = isFixed(row);
            const bool colFixed = isFixed(col);

            if (rowFixed)
              return;

            if (colFixed && row != col)
            {
              localRhs[static_cast<size_t>(row)] -=
                val * fixed[static_cast<size_t>(col)];
              return;
            }

            localT.emplace_back(row, col, val);
          };

          // ---------------- Local BFIs ----------------
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto& attrs = bfi.getAttributes();
            OpenMPIteration seq(mesh, bfi.getRegion());
            const Index d = seq.getDimension();
            const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& localT = tchunks[static_cast<size_t>(tid)];
              auto& localRhs = rhsChunks[static_cast<size_t>(tid)];

              std::unique_ptr<LocalBilinearFormIntegratorBaseType> integrator(bfi.copy());

#pragma omp for
              for (Index p = 0; p < count; ++p)
              {
                if (!seq.filter(p)) continue;
                if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                auto it = seq.getIterator(p);
                integrator->setPolytope(*it);

                const auto& rowsDOF = testFES.getDOFs(d, p);
                const auto& colsDOF = trialFES.getDOFs(d, p);

                for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
                {
                  const Index I = rowsDOF(i);
                  for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                  {
                    const Index J = colsDOF(j);
                    const ScalarType val = Math::conj(integrator->integrate(j, i));
                    sparse_entry(localT, localRhs, I, J, val);
                  }
                }
              }
            }
          }

          // ---------------- Global BFIs ----------------
          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();

            OpenMPIteration trialseq(mesh, bfi.getTrialRegion());
            OpenMPIteration testseq(mesh, bfi.getTestRegion());

            const Index td = testseq.getDimension();
            const Index tcount = testseq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& localT = tchunks[static_cast<size_t>(tid)];
              auto& localRhs = rhsChunks[static_cast<size_t>(tid)];

              std::unique_ptr<GlobalBilinearFormIntegratorBaseType> integrator(bfi.copy());

#pragma omp for
              for (Index te = 0; te < tcount; ++te)
              {
                if (!testseq.filter(te)) continue;
                if (!testAttrs.empty() && !testAttrs.count(mesh.getAttribute(td, te))) continue;

                const auto& rowsDOF = testFES.getDOFs(td, te);

                const Index rd = trialseq.getDimension();
                const Index rcount = trialseq.getCount();

                for (Index tr = 0; tr < rcount; ++tr)
                {
                  if (!trialseq.filter(tr)) continue;
                  if (!trialAttrs.empty() && !trialAttrs.count(mesh.getAttribute(rd, tr))) continue;

                  const auto& colsDOF = trialFES.getDOFs(rd, tr);

                  auto teIt = testseq.getIterator(te);
                  auto trIt = trialseq.getIterator(tr);

                  integrator->setPolytope(*trIt, *teIt);

                  for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
                  {
                    const Index I = rowsDOF(i);
                    for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                    {
                      const Index J = colsDOF(j);
                      const ScalarType val = Math::conj(integrator->integrate(j, i));
                      sparse_entry(localT, localRhs, I, J, val);
                    }
                  }
                }
              }
            }
          }

          // ---------------- Preassembled bilinear forms ----------------
          for (auto& bf : pb.getBFs())
          {
            const auto& op = bf.getOperator();
            for (int k = 0; k < op.outerSize(); ++k)
              for (typename OperatorType::InnerIterator it(op, k); it; ++it)
                sparse_entry(tchunks[0], rhsChunks[0], it.row(), it.col(), it.value());
          }

          // ---------------- LFIs ----------------
          for (auto& lfi : pb.getLFIs())
          {
            const auto& attrs = lfi.getAttributes();
            OpenMPIteration seq(mesh, lfi.getRegion());
            const Index d = seq.getDimension();
            const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& localRhs = rhsChunks[static_cast<size_t>(tid)];

              std::unique_ptr<LinearFormIntegratorBaseType> integrator(lfi.copy());

#pragma omp for
              for (Index p = 0; p < count; ++p)
              {
                if (!seq.filter(p)) continue;
                if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                auto it = seq.getIterator(p);
                integrator->setPolytope(*it);

                const auto& dofs = testFES.getDOFs(d, p);
                for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
                {
                  const Index I = dofs(l);
                  localRhs[static_cast<size_t>(I)] -= integrator->integrate(l);
                }
              }
            }
          }

          // Preassembled linear forms (serial)
          for (auto& lf : pb.getLFs())
            b -= lf.getVector();

          // ---------------- Reduce RHS chunks into b ----------------
          for (int tid = 0; tid < tc; ++tid)
          {
            const auto& localRhs = rhsChunks[static_cast<size_t>(tid)];
            for (size_t i = 0; i < rows; ++i)
              b.coeffRef(i) += localRhs[i];
          }

          // ---------------- Merge triplets ----------------
          size_t totalTriplets = 0;
          for (auto& v0 : tchunks) totalTriplets += v0.size();

          std::vector<Eigen::Triplet<ScalarType>> all;
          all.reserve(totalTriplets + rows);

          for (auto& v0 : tchunks)
          {
            all.insert(all.end(),
                       std::make_move_iterator(v0.begin()),
                       std::make_move_iterator(v0.end()));
            v0.clear();
          }

          // Inject identity rows for fixed dofs
          for (Index i = 0; i < static_cast<Index>(rows); ++i)
          {
            if (isFixed(i))
            {
              all.emplace_back(i, i, ScalarType(1));
              b.coeffRef(static_cast<size_t>(i)) = fixed[static_cast<size_t>(i)];
            }
          }

          A.resize(rows, cols);
          A.setFromTriplets(all.begin(), all.end());
        }
        else
        {
          // ---- Dense path: assemble then eliminate afterwards ----
          A.resize(rows, cols);
          A.setZero();

          std::vector<OperatorType> Achunks(static_cast<size_t>(tc));
          std::vector<std::vector<ScalarType>> rhsChunks(
            static_cast<size_t>(tc),
            std::vector<ScalarType>(rows, ScalarType(0))
          );

          for (int tid = 0; tid < tc; ++tid)
          {
            Achunks[static_cast<size_t>(tid)].resize(rows, cols);
            Achunks[static_cast<size_t>(tid)].setZero();
          }

          // Local BFIs
          for (auto& bfi : pb.getLocalBFIs())
          {
            const auto& attrs = bfi.getAttributes();
            OpenMPIteration seq(mesh, bfi.getRegion());
            const Index d = seq.getDimension();
            const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& Alocal = Achunks[static_cast<size_t>(tid)];

              std::unique_ptr<LocalBilinearFormIntegratorBaseType> integrator(bfi.copy());

#pragma omp for
              for (Index p = 0; p < count; ++p)
              {
                if (!seq.filter(p)) continue;
                if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                auto it = seq.getIterator(p);
                integrator->setPolytope(*it);

                const auto& rowsDOF = testFES.getDOFs(d, p);
                const auto& colsDOF = trialFES.getDOFs(d, p);

                for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
                {
                  const Index I = rowsDOF(i);
                  for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                  {
                    const Index J = colsDOF(j);
                    const ScalarType val = Math::conj(integrator->integrate(j, i));
                    if (val != ScalarType(0))
                      Alocal(I, J) += val;
                  }
                }
              }
            }
          }

          // Global BFIs
          for (auto& bfi : pb.getGlobalBFIs())
          {
            const auto& trialAttrs = bfi.getTrialAttributes();
            const auto& testAttrs  = bfi.getTestAttributes();

            OpenMPIteration trialseq(mesh, bfi.getTrialRegion());
            OpenMPIteration testseq(mesh, bfi.getTestRegion());

            const Index td = testseq.getDimension();
            const Index tcount = testseq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& Alocal = Achunks[static_cast<size_t>(tid)];

              std::unique_ptr<GlobalBilinearFormIntegratorBaseType> integrator(bfi.copy());

#pragma omp for
              for (Index te = 0; te < tcount; ++te)
              {
                if (!testseq.filter(te)) continue;
                if (!testAttrs.empty() && !testAttrs.count(mesh.getAttribute(td, te))) continue;

                const auto& rowsDOF = testFES.getDOFs(td, te);

                const Index rd = trialseq.getDimension();
                const Index rcount = trialseq.getCount();

                for (Index tr = 0; tr < rcount; ++tr)
                {
                  if (!trialseq.filter(tr)) continue;
                  if (!trialAttrs.empty() && !trialAttrs.count(mesh.getAttribute(rd, tr))) continue;

                  const auto& colsDOF = trialFES.getDOFs(rd, tr);

                  auto teIt = testseq.getIterator(te);
                  auto trIt = trialseq.getIterator(tr);

                  integrator->setPolytope(*trIt, *teIt);

                  for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
                  {
                    const Index I = rowsDOF(i);
                    for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                    {
                      const Index J = colsDOF(j);
                      const ScalarType val = Math::conj(integrator->integrate(j, i));
                      if (val != ScalarType(0))
                        Alocal(I, J) += val;
                    }
                  }
                }
              }
            }
          }

          // Preassembled bilinear forms (serial)
          for (auto& bf : pb.getBFs())
            A += bf.getOperator();

          // Reduce dense matrices
          for (int tid = 0; tid < tc; ++tid)
            A += Achunks[static_cast<size_t>(tid)];

          // LFIs
          for (auto& lfi : pb.getLFIs())
          {
            const auto& attrs = lfi.getAttributes();
            OpenMPIteration seq(mesh, lfi.getRegion());
            const Index d = seq.getDimension();
            const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& localRhs = rhsChunks[static_cast<size_t>(tid)];

              std::unique_ptr<LinearFormIntegratorBaseType> integrator(lfi.copy());

#pragma omp for
              for (Index p = 0; p < count; ++p)
              {
                if (!seq.filter(p)) continue;
                if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                auto it = seq.getIterator(p);
                integrator->setPolytope(*it);

                const auto& dofs = testFES.getDOFs(d, p);
                for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
                  localRhs[static_cast<size_t>(dofs(l))] -= integrator->integrate(l);
              }
            }
          }

          // Reduce RHS chunks into b
          for (int tid = 0; tid < tc; ++tid)
          {
            const auto& localRhs = rhsChunks[static_cast<size_t>(tid)];
            for (size_t i = 0; i < rows; ++i)
              b.coeffRef(i) += localRhs[i];
          }

          // Preassembled linear forms (serial)
          for (auto& lf : pb.getLFs())
            b -= lf.getVector();

          // Dense elimination afterwards
          for (Index idx = 0; idx < static_cast<Index>(rows); ++idx)
          {
            if (!isFixed(idx))
              continue;

            const ScalarType value = fixed[static_cast<size_t>(idx)];

            for (size_t r = 0; r < rows; ++r)
            {
              if (r == static_cast<size_t>(idx)) continue;
              b.coeffRef(r) -= A(r, idx) * value;
              A(r, idx) = ScalarType(0);
            }

            for (size_t c = 0; c < cols; ++c)
            {
              if (c == static_cast<size_t>(idx)) continue;
              A(idx, c) = ScalarType(0);
            }

            A(idx, idx) = ScalarType(1);
            b.coeffRef(static_cast<size_t>(idx)) = value;
          }
        }
      }

      OpenMP* copy() const noexcept override { return new OpenMP(*this); }

    private:
      std::optional<size_t> m_threadCount;
  };

  template <class LinearSystem, class U1, class U2, class U3, class ... Us>
  class OpenMP<
    LinearSystem,
    Variational::Problem<LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      using LinearSystemType = LinearSystem;

      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, U1, U2, U3, Us...>>;

      using InputType = typename Parent::InputType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using LinearFormIntegratorBaseType =
        Variational::LinearFormIntegratorBase<ScalarType>;

      OpenMP() = default;

      OpenMP(const OpenMP& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      OpenMP(OpenMP&& other)
        : Parent(std::move(other)),
          m_threadCount(other.m_threadCount)
      {}

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
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

        b.resize(nrows);
        b.setZero();

        constexpr bool IsSparse =
          std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>;

        const int tc =
          m_threadCount.has_value()
            ? static_cast<int>(m_threadCount.value())
            : omp_get_max_threads();

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

        // Mesh (assume common)
        const auto& mesh = [&]() -> const auto&
        {
          const void* addr = nullptr;
          us.apply([&](const auto& uref)
          {
            if (!addr)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace().getMesh());
          });
          assert(addr);
          return *static_cast<
            const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>*
          >(addr);
        }();

        // ------------------------------------------------------------------
        // Dirichlet: NaN sentinel (global indices!)
        // ------------------------------------------------------------------
        const size_t ndofs = std::max(nrows, ncols);
        std::vector<ScalarType> fixed(ndofs, Math::nan<ScalarType>());

        auto isFixed = [&](Index I) -> bool
        {
          const size_t k = static_cast<size_t>(I);
          return k < fixed.size() && !Math::isNaN(fixed[k]);
        };

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          for (const auto& [local, value] : dbc.getDOFs())
          {
            const Index I = static_cast<Index>(uOff + static_cast<size_t>(local));
            fixed[static_cast<size_t>(I)] = static_cast<ScalarType>(value);
          }
        }

        // ------------------------------------------------------------------
        // Thread-local accumulators
        // ------------------------------------------------------------------
        std::vector<std::vector<Eigen::Triplet<ScalarType>>> tchunks(static_cast<size_t>(tc));
        std::vector<std::vector<ScalarType>> rhsChunks(
          static_cast<size_t>(tc),
          std::vector<ScalarType>(nrows, ScalarType(0))
        );

        std::vector<OperatorType> Achunks;
        if constexpr (!IsSparse)
        {
          Achunks.resize(static_cast<size_t>(tc));
          for (int tid = 0; tid < tc; ++tid)
          {
            Achunks[static_cast<size_t>(tid)].resize(nrows, ncols);
            Achunks[static_cast<size_t>(tid)].setZero();
          }
        }
        else
        {
          if constexpr (IsSparse)
            tchunks.shrink_to_fit(); // no-op typically; keeps symmetry with dense path
        }

        auto sparse_entry = [&](int tid, Index row, Index col, ScalarType val)
        {
          if (val == ScalarType(0))
            return;

          const bool rowFixed = isFixed(row);
          const bool colFixed = isFixed(col);

          if (rowFixed)
            return;

          if (colFixed && row != col)
          {
            rhsChunks[static_cast<size_t>(tid)][static_cast<size_t>(row)] -=
              val * fixed[static_cast<size_t>(col)];
            return;
          }

          tchunks[static_cast<size_t>(tid)].emplace_back(row, col, val);
        };

        // ------------------------------------------------------------------
        // Local BFIs (BUGFIX: type-safe per-block FES access; no casts)
        // ------------------------------------------------------------------
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
          const Index d = seq.getDimension();
          const Index count = seq.getCount();

          us.iapply([&](size_t ui, const auto& uref)
          {
            if (ui != uBlock) return;
            const auto& uFES = uref.get().getFiniteElementSpace();

            vs.iapply([&](size_t vi, const auto& vref)
            {
              if (vi != vBlock) return;
              const auto& vFES = vref.get().getFiniteElementSpace();

#pragma omp parallel num_threads(tc)
              {
                const int tid = omp_get_thread_num();

                std::unique_ptr<LocalBilinearFormIntegratorBaseType> integrator(bfi.copy());
                OperatorType* Alocal = nullptr;
                if constexpr (!IsSparse)
                  Alocal = &Achunks[static_cast<size_t>(tid)];

#pragma omp for
                for (Index p = 0; p < count; ++p)
                {
                  if (!seq.filter(p)) continue;
                  if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                  auto it = seq.getIterator(p);
                  integrator->setPolytope(*it);

                  const auto& rows = vFES.getDOFs(d, p);
                  const auto& cols = uFES.getDOFs(d, p);

                  for (size_t i = 0; i < static_cast<size_t>(rows.size()); ++i)
                  {
                    const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows(i)));
                    for (size_t j = 0; j < static_cast<size_t>(cols.size()); ++j)
                    {
                      const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols(j)));
                      const ScalarType val = Math::conj(integrator->integrate(j, i));
                      if (val == ScalarType(0))
                        continue;

                      if constexpr (IsSparse)
                        sparse_entry(tid, I, J, val);
                      else
                        (*Alocal)(I, J) += val;
                    }
                  }
                }
              } // omp parallel
            });
          });
        }

        // ------------------------------------------------------------------
        // Global BFIs (BUGFIX: type-safe per-block FES access; no casts)
        // ------------------------------------------------------------------
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
          OpenMPIteration testseq(mesh, bfi.getTestRegion());

          const Index td = testseq.getDimension();
          const Index tcount = testseq.getCount();

          us.iapply([&](size_t ui, const auto& uref)
          {
            if (ui != uBlock) return;
            const auto& uFES = uref.get().getFiniteElementSpace();

            vs.iapply([&](size_t vi, const auto& vref)
            {
              if (vi != vBlock) return;
              const auto& vFES = vref.get().getFiniteElementSpace();

#pragma omp parallel num_threads(tc)
              {
                const int tid = omp_get_thread_num();

                std::unique_ptr<GlobalBilinearFormIntegratorBaseType> integrator(bfi.copy());
                OperatorType* Alocal = nullptr;
                if constexpr (!IsSparse)
                  Alocal = &Achunks[static_cast<size_t>(tid)];

#pragma omp for
                for (Index te = 0; te < tcount; ++te)
                {
                  if (!testseq.filter(te)) continue;
                  if (!testAttrs.empty() && !testAttrs.count(mesh.getAttribute(td, te))) continue;

                  const auto& rows = vFES.getDOFs(td, te);

                  const Index rd = trialseq.getDimension();
                  const Index rcount = trialseq.getCount();

                  for (Index tr = 0; tr < rcount; ++tr)
                  {
                    if (!trialseq.filter(tr)) continue;
                    if (!trialAttrs.empty() && !trialAttrs.count(mesh.getAttribute(rd, tr))) continue;

                    const auto& cols = uFES.getDOFs(rd, tr);

                    auto teIt = testseq.getIterator(te);
                    auto trIt = trialseq.getIterator(tr);

                    integrator->setPolytope(*trIt, *teIt);

                    for (size_t i = 0; i < static_cast<size_t>(rows.size()); ++i)
                    {
                      const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows(i)));
                      for (size_t j = 0; j < static_cast<size_t>(cols.size()); ++j)
                      {
                        const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols(j)));
                        const ScalarType val = Math::conj(integrator->integrate(j, i));
                        if (val == ScalarType(0))
                          continue;

                        if constexpr (IsSparse)
                          sparse_entry(tid, I, J, val);
                        else
                          (*Alocal)(I, J) += val;
                      }
                    }
                  }
                }
              } // omp parallel
            });
          });
        }

        // ------------------------------------------------------------------
        // LFIs (type-safe test FES; thread-local rhs)
        // ------------------------------------------------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          OpenMPIteration seq(mesh, lfi.getRegion());
          const Index d = seq.getDimension();
          const Index count = seq.getCount();

          vs.iapply([&](size_t vi, const auto& vref)
          {
            if (vi != vBlock) return;
            const auto& vFES = vref.get().getFiniteElementSpace();

#pragma omp parallel num_threads(tc)
            {
              const int tid = omp_get_thread_num();
              auto& localRhs = rhsChunks[static_cast<size_t>(tid)];

              std::unique_ptr<LinearFormIntegratorBaseType> integrator(lfi.copy());

#pragma omp for
              for (Index p = 0; p < count; ++p)
              {
                if (!seq.filter(p)) continue;
                if (!attrs.empty() && !attrs.count(mesh.getAttribute(d, p))) continue;

                auto it = seq.getIterator(p);
                integrator->setPolytope(*it);

                const auto& dofs = vFES.getDOFs(d, p);
                for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
                {
                  const Index I = static_cast<Index>(vOff + static_cast<size_t>(dofs(l)));
                  localRhs[static_cast<size_t>(I)] -= integrator->integrate(l);
                }
              }
            } // omp parallel
          });
        }

        // ------------------------------------------------------------------
        // Reduce RHS chunks into b
        // ------------------------------------------------------------------
        for (int tid = 0; tid < tc; ++tid)
        {
          const auto& localRhs = rhsChunks[static_cast<size_t>(tid)];
          for (size_t i = 0; i < nrows; ++i)
            b.coeffRef(i) += localRhs[i];
        }

        // Preassembled LFs (serial)
        for (auto& lf : pb.getLFs())
          b -= lf.getVector();

        // ------------------------------------------------------------------
        // Finalize operator
        // ------------------------------------------------------------------
        if constexpr (IsSparse)
        {
          // Preassembled BFs -> apply elimination through sparse_entry (tid=0)
          for (auto& bf : pb.getBFs())
          {
            const auto& op = bf.getOperator();
            for (int k = 0; k < op.outerSize(); ++k)
              for (typename OperatorType::InnerIterator it(op, k); it; ++it)
                sparse_entry(0, it.row(), it.col(), it.value());
          }

          // Merge triplets
          size_t totalTriplets = 0;
          for (auto& v0 : tchunks) totalTriplets += v0.size();

          std::vector<Eigen::Triplet<ScalarType>> all;
          all.reserve(totalTriplets + nrows);

          for (auto& v0 : tchunks)
          {
            all.insert(all.end(),
                       std::make_move_iterator(v0.begin()),
                       std::make_move_iterator(v0.end()));
            v0.clear();
          }

          // Identity rows for fixed dofs
          for (Index I = 0; I < static_cast<Index>(nrows); ++I)
          {
            if (isFixed(I))
            {
              all.emplace_back(I, I, ScalarType(1));
              b.coeffRef(static_cast<size_t>(I)) = fixed[static_cast<size_t>(I)];
            }
          }

          A.resize(nrows, ncols);
          A.setFromTriplets(all.begin(), all.end());
        }
        else
        {
          // Reduce dense thread-local matrices
          A.resize(nrows, ncols);
          A.setZero();
          for (int tid = 0; tid < tc; ++tid)
            A += Achunks[static_cast<size_t>(tid)];

          // Add preassembled dense BFs (serial)
          for (auto& bf : pb.getBFs())
            A += bf.getOperator();

          // Dense elimination after assembly
          for (Index idx = 0; idx < static_cast<Index>(nrows); ++idx)
          {
            if (!isFixed(idx))
              continue;

            const ScalarType value = fixed[static_cast<size_t>(idx)];

            for (size_t r = 0; r < nrows; ++r)
            {
              if (r == static_cast<size_t>(idx)) continue;
              b.coeffRef(r) -= A(r, idx) * value;
              A(r, idx) = ScalarType(0);
            }

            for (size_t c = 0; c < ncols; ++c)
            {
              if (c == static_cast<size_t>(idx)) continue;
              A(idx, c) = ScalarType(0);
            }

            A(idx, idx) = ScalarType(1);
            b.coeffRef(static_cast<size_t>(idx)) = value;
          }
        }
      }


      OpenMP* copy() const noexcept override { return new OpenMP(*this); }

    private:
      std::optional<size_t> m_threadCount;
  };
}

#endif
