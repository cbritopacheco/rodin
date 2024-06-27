/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_MULTITHREADED_H
#define RODIN_ASSEMBLY_MULTITHREADED_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Kernels.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Threads/Mutex.h"
#include "Rodin/Threads/ThreadPool.h"

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
  namespace Internal
  {
    class MultithreadedIteration
    {
      public:
        MultithreadedIteration(const Geometry::MeshBase& mesh, Variational::Integrator::Region);

        Geometry::PolytopeIterator getIterator(Index i) const;

        size_t getDimension() const;

        size_t getCount() const;

        bool filter(Index i) const;

      private:
        std::reference_wrapper<const Geometry::MeshBase> m_mesh;
        Variational::Integrator::Region m_region;
    };
  }

  template <class TrialFES, class TestFES>
  class Multithreaded<
    std::vector<Eigen::Triplet<Scalar>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<
        std::vector<Eigen::Triplet<Scalar>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
  {
    public:
      using Parent =
        AssemblyBase<
          std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>;

      using InputType = typename Parent::InputType;

      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

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

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_pool(
            std::visit(
              [](auto&& arg) -> std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>
              {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::reference_wrapper<Threads::ThreadPool>>)
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(arg);
                else
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(
                      std::in_place_type_t<Threads::ThreadPool>(), arg.getThreadCount());
              }, other.m_pool))
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_pool(std::move(other.getThreadPool()))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const InputType& input) const override
      {
        using TripletVector = std::vector<Eigen::Triplet<Scalar>>;
        const size_t capacity = input.getTestFES().getSize() * std::log(input.getTrialFES().getSize());
        TripletVector res;
        res.clear();
        res.reserve(capacity);
        const size_t threadCount = getThreadPool().getThreadCount();
        const auto& mesh = input.getTestFES().getMesh();
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::MultithreadedIteration seq(mesh, bfi.getRegion());
          const size_t d = seq.getDimension();
          auto loop =
            [&](const Index start, const Index end)
            {
              tl_lbfi.reset(bfi.copy());
              tl_triplets.reserve(capacity / threadCount);
              for (Index i = start; i < end; ++i)
              {
                if (seq.filter(i))
                {
                  if (attrs.size() == 0 || attrs.count(mesh.getAttribute(d, i)))
                  {
                    const auto it = seq.getIterator(i);
                    const auto& trialDOFs = input.getTrialFES().getDOFs(d, i);
                    const auto& testDOFs = input.getTestFES().getDOFs(d, i);
                    tl_lbfi->assemble(*it);
                    Math::Kernels::add(tl_triplets, tl_lbfi->getMatrix(), testDOFs, trialDOFs);
                  }
                }
              }
              m_mutex.lock();
              res.insert(res.end(),
                  std::make_move_iterator(tl_triplets.begin()),
                  std::make_move_iterator(tl_triplets.end()));
              m_mutex.unlock();
              tl_triplets.clear();
            };

          if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          {
            auto& threadPool = std::get<Threads::ThreadPool>(m_pool);
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
          else
          {
            auto& threadPool = std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::MultithreadedIteration testseq(mesh, bfi.getTestRegion());
          const size_t d = testseq.getDimension();
          auto loop =
            [&](const Index start, const Index end)
            {
              tl_gbfi.reset(bfi.copy());
              tl_triplets.reserve(capacity / threadCount);
              for (Index i = start; i < end; ++i)
              {
                if (testseq.filter(i))
                {
                  if (testAttrs.size() == 0 || testAttrs.count(mesh.getAttribute(d, i)))
                  {
                    const auto teIt = testseq.getIterator(i);
                    Internal::SequentialIteration trialseq{ mesh, tl_gbfi->getTrialRegion() };
                    for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                    {
                      if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                      {
                        const auto& trialDOFs = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                        const auto& testDOFs = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                        tl_gbfi->assemble(*trIt, *teIt);
                        Math::Kernels::add(tl_triplets, tl_gbfi->getMatrix(), testDOFs, trialDOFs);
                      }
                    }
                  }
                }
              }
              m_mutex.lock();
              res.insert(res.end(),
                  std::make_move_iterator(tl_triplets.begin()),
                  std::make_move_iterator(tl_triplets.end()));
              m_mutex.unlock();
              tl_triplets.clear();
            };

          if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          {
            auto& threadPool = std::get<Threads::ThreadPool>(m_pool);
            threadPool.pushLoop(0, testseq.getCount(), loop);
            threadPool.waitForTasks();
          }
          else
          {
            auto& threadPool = std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
            threadPool.pushLoop(0, testseq.getCount(), loop);
            threadPool.waitForTasks();
          }
        }
        return res;
      }

      const Threads::ThreadPool& getThreadPool() const
      {
        if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          return std::get<Threads::ThreadPool>(m_pool);
        else
          return std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local std::vector<Eigen::Triplet<Scalar>> tl_triplets;
      static thread_local std::unique_ptr<Variational::LocalBilinearFormIntegratorBase> tl_lbfi;
      static thread_local std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase> tl_gbfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

  template <class TrialFES, class TestFES>
  thread_local
  std::vector<Eigen::Triplet<Scalar>>
  Multithreaded<
    std::vector<Eigen::Triplet<Scalar>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>::tl_triplets;

  template <class TrialFES, class TestFES>
  thread_local
  std::unique_ptr<Variational::LocalBilinearFormIntegratorBase>
  Multithreaded<
    std::vector<Eigen::Triplet<Scalar>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>::tl_lbfi;

  template <class TrialFES, class TestFES>
  thread_local
  std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase>
  Multithreaded<
    std::vector<Eigen::Triplet<Scalar>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>::tl_gbfi;

  /**
   * @brief Multithreaded assembly of the Math::SparseMatrix<Scalar> associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::SparseMatrix<Scalar>,
    Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Scalar>>>
    : public AssemblyBase<
        Math::SparseMatrix<Scalar>,
        Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Scalar>>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::SparseMatrix<Scalar>,
          Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Scalar>>>;

      using InputType = typename Parent::InputType;
      using OperatorType = Math::SparseMatrix<Scalar>;

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
        : m_assembly(pool)
      {}

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
      OperatorType execute(const InputType& input) const override
      {
        const auto triplets = m_assembly.execute({
            input.getTrialFES(), input.getTestFES(),
            input.getLocalBFIs(), input.getGlobalBFIs() });
        OperatorType res(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      Multithreaded<
        std::vector<Eigen::Triplet<Scalar>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> m_assembly;
  };

  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::Matrix<Scalar>,
    Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>
    : public AssemblyBase<
        Math::Matrix<Scalar>,
        Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::Matrix<Scalar>,
          Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>;
      using InputType = typename Parent::InputType;
      using OperatorType = Math::Matrix<Scalar>;

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

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_pool(
            std::visit(
              [](auto&& arg) -> std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>
              {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::reference_wrapper<Threads::ThreadPool>>)
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(arg);
                else
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(
                      std::in_place_type_t<Threads::ThreadPool>(), arg.getThreadCount());
              }, other.m_pool))
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_pool(std::move(other.getThreadPool()))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const InputType& input) const override
      {
        Math::Matrix<Scalar> res(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setZero();
        auto& threadPool = getThreadPool();
        const auto& mesh = input.getTestFES().getMesh();
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::MultithreadedIteration seq(mesh, bfi.getRegion());
          const size_t d = seq.getDimension();
          auto loop =
            [&](const Index start, const Index end)
            {
              tl_lbfi.reset(bfi.copy());
              tl_res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
              tl_res.setZero();
              for (Index i = start; i < end; ++i)
              {
                if (seq.filter(i))
                {
                  if (attrs.size() == 0 || attrs.count(mesh.getAttribute(d, i)))
                  {
                    const auto it = seq.getIterator(i);
                    const auto& trialDOFs = input.getTrialFES().getDOFs(d, i);
                    const auto& testDOFs = input.getTestFES().getDOFs(d, i);
                    tl_lbfi->assemble(*it);
                    Math::Kernels::add(tl_res, tl_lbfi->getMatrix(), testDOFs, trialDOFs);
                  }
                }
              }
              m_mutex.lock();
              res += tl_res;
              m_mutex.unlock();
            };

          if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          {
            auto& threadPool = std::get<Threads::ThreadPool>(m_pool);
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
          else
          {
            auto& threadPool = std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::MultithreadedIteration testseq(mesh, bfi.getTestRegion());
          const size_t d = testseq.getDimension();
          const auto loop =
            [&](const Index start, const Index end)
            {
              tl_gbfi.reset(bfi.copy());
              tl_res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
              tl_res.setZero();
              for (Index i = start; i < end; ++i)
              {
                if (testseq.filter(i))
                {
                  if (testAttrs.size() == 0 || testAttrs.count(mesh.getAttribute(d, i)))
                  {
                    const auto teIt = testseq.getIterator(i);
                    Internal::SequentialIteration trialseq{ mesh, tl_gbfi->getTrialRegion() };
                    for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                    {
                      if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                      {
                        const auto& trialDOFs = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                        const auto& testDOFs = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                        tl_gbfi->assemble(*trIt, *teIt);
                        Math::Kernels::add(tl_res, tl_gbfi->getMatrix(), testDOFs, trialDOFs);
                      }
                    }
                  }
                }
              }
              m_mutex.lock();
              res += tl_res;
              m_mutex.unlock();
            };

          if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          {
            auto& threadPool = std::get<Threads::ThreadPool>(m_pool);
            threadPool.pushLoop(0, testseq.getCount(), loop);
            threadPool.waitForTasks();
          }
          else
          {
            auto& threadPool = std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
            threadPool.pushLoop(0, testseq.getCount(), loop);
            threadPool.waitForTasks();
          }
        }
        return res;
      }

      const Threads::ThreadPool& getThreadPool() const
      {
        if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          return std::get<Threads::ThreadPool>(m_pool);
        else
          return std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local Math::Matrix<Scalar> tl_res;
      static thread_local std::unique_ptr<Variational::LocalBilinearFormIntegratorBase> tl_lbfi;
      static thread_local std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase> tl_gbfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

   template <class TrialFES, class TestFES>
   thread_local
   Math::Matrix<Scalar>
   Multithreaded<Math::Matrix<Scalar>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>::tl_res;

   template <class TrialFES, class TestFES>
   thread_local
   std::unique_ptr<Variational::LocalBilinearFormIntegratorBase>
   Multithreaded<Math::Matrix<Scalar>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>::tl_lbfi;

   template <class TrialFES, class TestFES>
   thread_local
   std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase>
   Multithreaded<Math::Matrix<Scalar>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>>>::tl_gbfi;

  /**
   * @brief %Multithreaded assembly of the Math::Vector<Scalar> associated to a LinearFormBase
   * object.
   */
  template <class FES>
  class Multithreaded<Math::Vector<Scalar>, Variational::LinearForm<FES, Math::Vector<Scalar>>>
    : public AssemblyBase<Math::Vector<Scalar>, Variational::LinearForm<FES, Math::Vector<Scalar>>>
  {
    public:
      using Parent = AssemblyBase<Math::Vector<Scalar>, Variational::LinearForm<FES, Math::Vector<Scalar>>>;
      using InputType = typename Parent::InputType;
      using VectorType = Math::Vector<Scalar>;

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

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_pool(
            std::visit(
              [](auto&& arg) -> std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>
              {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::reference_wrapper<Threads::ThreadPool>>)
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(arg);
                else
                  return std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>>(
                      std::in_place_type_t<Threads::ThreadPool>(), arg.getThreadCount());
              }, other.m_pool))
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_pool(std::move(other.getThreadPool()))
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const InputType& input) const override
      {
        VectorType res(input.getFES().getSize());
        res.setZero();
        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          Internal::MultithreadedIteration seq(mesh, lfi.getRegion());
          const size_t d = seq.getDimension();
          const auto loop =
            [&](const Index start, const Index end)
            {
              tl_lfi.reset(lfi.copy());
              tl_res.resize(input.getFES().getSize());
              tl_res.setZero();
              for (Index i = start; i < end; ++i)
              {
                if (seq.filter(i))
                {
                  if (attrs.size() == 0 || attrs.count(mesh.getAttribute(d, i)))
                  {
                    const auto it = seq.getIterator(i);
                    const auto& dofs = input.getFES().getDOFs(d, i);
                    tl_lfi->assemble(*it);
                    Math::Kernels::add(tl_res, tl_lfi->getVector(), dofs);
                  }
                }
              }
              m_mutex.lock();
              res += tl_res;
              m_mutex.unlock();
            };

          if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          {
            auto& threadPool = std::get<Threads::ThreadPool>(m_pool);
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
          else
          {
            auto& threadPool = std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
            threadPool.pushLoop(0, seq.getCount(), loop);
            threadPool.waitForTasks();
          }
        }
        return res;
      }

      const Threads::ThreadPool& getThreadPool() const
      {
        if (std::holds_alternative<Threads::ThreadPool>(m_pool))
          return std::get<Threads::ThreadPool>(m_pool);
        else
          return std::get<std::reference_wrapper<Threads::ThreadPool>>(m_pool).get();
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local VectorType tl_res;
      static thread_local std::unique_ptr<Variational::LinearFormIntegratorBase> tl_lfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

  template <class FES>
  thread_local
  Math::Vector<Scalar> Multithreaded<Math::Vector<Scalar>, Variational::LinearForm<FES, Math::Vector<Scalar>>>::tl_res;

  template <class FES>
  thread_local
  std::unique_ptr<Variational::LinearFormIntegratorBase>
  Multithreaded<Math::Vector<Scalar>, Variational::LinearForm<FES, Math::Vector<Scalar>>>::tl_lfi;
}

#endif

