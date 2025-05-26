/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_MULTITHREADED_H
#define RODIN_ASSEMBLY_MULTITHREADED_H

#include <omp.h>

#include "Rodin/Math/Vector.h"

#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Rodin/Utility/Overloaded.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"
#include "Sequential.h"

namespace Rodin::Assembly
{
  template <>
  class MultithreadedIteration<Geometry::Mesh<Context::Local>>
  {
    public:
      using MeshType = Geometry::Mesh<Context::Local>;

      MultithreadedIteration(const MeshType& mesh, Variational::Integrator::Region);

      Geometry::PolytopeIterator getIterator(Index i) const;

      size_t getDimension() const;

      size_t getCount() const;

      bool filter(Index i) const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      Variational::Integrator::Region m_region;
  };

  MultithreadedIteration(
      const Geometry::Mesh<Context::Local>& mesh, Variational::Integrator::Region)
    -> MultithreadedIteration<Geometry::Mesh<Context::Local>>;

  template <class TrialFES, class TestFES>
  class Multithreaded<
    std::vector<Eigen::Triplet<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
    Variational::BilinearForm<TrialFES, TestFES,
      std::vector<Eigen::Triplet<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>> final
    : public AssemblyBase<
        std::vector<Eigen::Triplet<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
        Variational::BilinearForm<TrialFES, TestFES,
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
        Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        const size_t capacity = input.getTestFES().getSize() * std::log(input.getTrialFES().getSize());
        res.clear();
        res.reserve(capacity);
        const auto& mesh = input.getTestFES().getMesh();
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MultithreadedIteration seq(mesh, bfi.getRegion());
          const Index d = seq.getDimension();
          const Index count = seq.getCount();
          const int tc = m_threadCount.has_value() ? m_threadCount.value() : omp_get_max_threads();

#pragma omp parallel num_threads(tc)
          {
            // one integrator + triplet buffer per thread
            std::unique_ptr<LocalBilinearFormIntegratorBaseType> integrator(bfi.copy());
            OperatorType local;
            local.reserve(capacity / tc);

#pragma omp for nowait
            for (Index i = 0; i < count; ++i)
            {
              if (seq.filter(i))
              {
                if (attrs.empty() || attrs.count(mesh.getAttribute(d, i)))
                {
                  auto it = seq.getIterator(i);
                  integrator->setPolytope(*it);

                  const auto& rows = input.getTestFES().getDOFs(d, i);
                  const auto& cols = input.getTrialFES().getDOFs(d, i);
                  for (size_t r = 0; r < rows.size(); ++r)
                  {
                    for (size_t c = 0; c < cols.size(); ++c)
                    {
                      const ScalarType s = integrator->integrate(c, r);
                      if (s != ScalarType(0))
                        local.emplace_back(rows(r), cols(c), s);
                    }
                  }
                }
              }
            }

#pragma omp critical
            {
              res.insert(
                res.end(),
                std::make_move_iterator(local.begin()),
                std::make_move_iterator(local.end()));
            }
          } // end parallel
        }
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      Multithreaded& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };

  /**
   * @brief Multithreaded assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::SparseMatrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
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

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Multithreaded() = default;

      Multithreaded(size_t threadCount)
        : m_assembly(threadCount)
      {
        assert(threadCount > 0);
      }

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_assembly(other.m_assembly)
      {}

      Multithreaded(Multithreaded&& other)
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

      Multithreaded& setThreadCount(size_t threadCount) noexcept
      {
        m_assembly.setThreadCount(threadCount);
        return *this;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      Multithreaded<
        std::vector<Eigen::Triplet<ScalarType>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<ScalarType>>>> m_assembly;
  };

  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::Matrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
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
        Variational::BilinearForm<TrialFES, TestFES,
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

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      void execute(OperatorType& res, const InputType& input) const override
      {
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setZero();

        const auto& mesh = input.getTestFES().getMesh();
        const int tc = m_threadCount.has_value() ? m_threadCount.value() : omp_get_max_threads();

        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          MultithreadedIteration seq(mesh, bfi.getRegion());
          const Index d     = seq.getDimension();
          const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
          {
            auto lbfi = std::unique_ptr<LocalBilinearFormIntegratorBaseType>(bfi.copy());
            OperatorType local;
            local.resize(res.rows(), res.cols());
            local.setZero();

#pragma omp for nowait
            for (Index i = 0; i < count; ++i)
            {
              if (!seq.filter(i)) continue;
              if (!attrs.empty() &&
                  !attrs.count(mesh.getAttribute(d, i)))
                continue;

              auto it = seq.getIterator(i);
              lbfi->setPolytope(*it);

              const auto& rows = input.getTestFES().getDOFs(d, i);
              const auto& cols = input.getTrialFES().getDOFs(d, i);
              for (size_t r = 0; r < rows.size(); ++r)
                for (size_t c = 0; c < cols.size(); ++c)
                  local(rows(r), cols(c)) += lbfi->integrate(c, r);
            }

#pragma omp critical
            {
              res += local;
            }
          }
        }

        // --- Global (coupling) contributions ---
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& testAttrs  = bfi.getTestAttributes();
          const auto& trialAttrs = bfi.getTrialAttributes();
          MultithreadedIteration testseq(mesh, bfi.getTestRegion());
          const Index d = testseq.getDimension();
          const Index count = testseq.getCount();

#pragma omp parallel num_threads(tc)
          {
            auto gbfi = std::unique_ptr<GlobalBilinearFormIntegratorBaseType>(bfi.copy());
            OperatorType local;
            local.resize(res.rows(), res.cols());
            local.setZero();

#pragma omp for nowait
            for (Index i = 0; i < count; ++i)
            {
              if (testseq.filter(i))
              {
                if (testAttrs.empty() || testAttrs.count(mesh.getAttribute(d, i)))
                {
                  auto teIt = testseq.getIterator(i);
                  SequentialIteration trialseq{ mesh, gbfi->getTrialRegion() };

                  for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                  {
                    if (trialAttrs.empty() || trialAttrs.count(trIt->getAttribute()))
                    {
                      gbfi->setPolytope(*trIt, *teIt);

                      const auto& rows = input.getTestFES().getDOFs(d, teIt->getIndex());
                      const auto& cols = input.getTrialFES().getDOFs(d, trIt->getIndex());
                      for (size_t r = 0; r < rows.size(); ++r)
                        for (size_t c = 0; c < cols.size(); ++c)
                          local(rows(r), cols(c)) += gbfi->integrate(c, r);
                    }
                  }
                }
              }
            }

#pragma omp critical
            {
              res += local;
            }
          }
        }
      }

      Multithreaded& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };

  /**
   * @brief %Multithreaded assembly of the Math::Vector associated to a
   * LinearForm object.
   */
  template <class FES>
  class Multithreaded<
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

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount))
      {}

      void execute(VectorType& res, const InputType& input) const override
      {
        // initialize global vector
        res.resize(input.getFES().getSize());
        res.setZero();

        const auto& mesh = input.getFES().getMesh();
        const int tc = static_cast<int>(getThreadCount());

        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          MultithreadedIteration seq(mesh, lfi.getRegion());
          const Index d = seq.getDimension();
          const Index count = seq.getCount();

#pragma omp parallel num_threads(tc)
          {
            auto integrator =
              std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>>(lfi.copy());
            VectorType local;
            local.resize(res.size());
            local.setZero();

#pragma omp for nowait
            for (Index i = 0; i < count; ++i)
            {
              if (seq.filter(i))
              {
                if (attrs.empty() || attrs.count(mesh.getAttribute(d, i)))
                {
                  auto it = seq.getIterator(i);
                  integrator->setPolytope(*it);

                  const auto& dofs = input.getFES().getDOFs(d, i);
                  for (size_t k = 0; k < dofs.size(); ++k)
                    local(dofs(k)) += integrator->integrate(k);
                }
              }
            }

#pragma omp critical
            {
              res += local;
            }
          }
        }
      }

      Multithreaded& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };

  template <class Scalar, class FES, class ValueDerived>
  class Multithreaded<
    IndexMap<Scalar>,
    Variational::DirichletBC<
      Variational::TrialFunction<FES>, Variational::FunctionBase<ValueDerived>>> final
    : public AssemblyBase<
        IndexMap<Scalar>,
        Variational::DirichletBC<
          Variational::TrialFunction<FES>, Variational::FunctionBase<ValueDerived>>>
  {
    public:
      using FESType = FES;

      using TrialFunctionType = Variational::TrialFunction<FES>;

      using ValueType = Variational::FunctionBase<ValueDerived>;

      using DirichletBCType = Variational::DirichletBC<TrialFunctionType, ValueType>;

      using Parent = AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      using FESRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using InputType = typename Parent::InputType;

      Multithreaded() = default;

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
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
        const size_t threadCount = m_threadCount.has_value() ? m_threadCount.value() : omp_get_max_threads();

#pragma omp parallel num_threads(threadCount)
        {
          IndexMap<Scalar> partial;

#pragma omp for nowait
          for (Index i = 0; i < faceCount; ++i)
          {
            if (mesh.isBoundary(i))
            {
              if (essBdr.size() == 0 || essBdr.count(mesh.getAttribute(faceDim, i)))
              {
                const auto& fe = fes.getFiniteElement(faceDim, i);
                const auto& mapping =
                  fes.getMapping({ faceDim, i }, value.template cast<FESRangeType>());
                for (Index local = 0; local < fe.getCount(); local++)
                {
                  const Index global = fes.getGlobalIndex({ faceDim, i }, local);
                  auto find = partial.find(global);
                  if (find == partial.end())
                  {
                    const auto& lf = fe.getLinearForm(local);
                    const auto s = lf(mapping);
                    partial.insert(find, std::pair{ global, s });
                  }
                }
              }
            }
          }

#pragma omp critical
          {
            res.insert(boost::container::ordered_unique_range, partial.begin(), partial.end());
          }
        }
      }

      size_t getThreadCount() const noexcept
      {
        return m_threadCount.value_or(omp_get_max_threads());
      }

      Multithreaded& setThreadCount(size_t threadCount) noexcept
      {
        m_threadCount = threadCount;
        return *this;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      std::optional<size_t> m_threadCount;
  };
}

#endif

