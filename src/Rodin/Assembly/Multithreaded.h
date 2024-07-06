/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_MULTITHREADED_H
#define RODIN_ASSEMBLY_MULTITHREADED_H

#include "Rodin/Math/Vector.h"

#include "Rodin/Math/Matrix.h"
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
    std::vector<Eigen::Triplet<Real>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>
    : public AssemblyBase<
        std::vector<Eigen::Triplet<Real>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using BilinearFormType =
        Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

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
        using TripletVector = std::vector<Eigen::Triplet<ScalarType>>;
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
                    tl_lbfi->setPolytope(*it);
                    const auto& rows = input.getTestFES().getDOFs(d, i);
                    const auto& cols = input.getTrialFES().getDOFs(d, i);
                    for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                    {
                      for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                      {
                        const ScalarType s = tl_lbfi->integrate(m, l);
                        if (s != ScalarType(0))
                          tl_triplets.emplace_back(rows(l), cols(m), s);
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
                        tl_gbfi->setPolytope(*trIt, *teIt);
                        const auto& rows = input.getTestFES().getDOFs(d, i);
                        const auto& cols = input.getTrialFES().getDOFs(d, i);
                        for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                        {
                          for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                          {
                            const ScalarType s = tl_gbfi->integrate(m, l);
                            if (s != ScalarType(0))
                              tl_triplets.emplace_back(rows(l), cols(m), s);
                          }
                        }
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
      static thread_local OperatorType tl_triplets;
      static thread_local std::unique_ptr<LocalBilinearFormIntegratorBaseType>  tl_lbfi;
      static thread_local std::unique_ptr<GlobalBilinearFormIntegratorBaseType> tl_gbfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

  template <class TrialFES, class TestFES>
  thread_local
  std::vector<Eigen::Triplet<Real>>
  Multithreaded<
    std::vector<Eigen::Triplet<Real>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>::tl_triplets;

  template <class TrialFES, class TestFES>
  thread_local
  std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<Real>>
  Multithreaded<
    std::vector<Eigen::Triplet<Real>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>::tl_lbfi;

  template <class TrialFES, class TestFES>
  thread_local
  std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase<Real>>
  Multithreaded<
    std::vector<Eigen::Triplet<Real>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>::tl_gbfi;

  /**
   * @brief Multithreaded assembly of the Math::SparseMatrix<Real> associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::SparseMatrix<Real>,
    Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>> final
    : public AssemblyBase<
        Math::SparseMatrix<Real>,
        Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

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
        std::vector<Eigen::Triplet<ScalarType>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<ScalarType>>>> m_assembly;
  };

  template <class TrialFES, class TestFES>
  class Multithreaded<
    Math::Matrix<Real>,
    Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>
    : public AssemblyBase<
        Math::Matrix<Real>,
        Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = Math::Matrix<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

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
        OperatorType res(input.getTestFES().getSize(), input.getTrialFES().getSize());
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
                    tl_lbfi->setPolytope(*it);
                    const auto& rows = input.getTestFES().getDOFs(d, i);
                    const auto& cols = input.getTrialFES().getDOFs(d, i);
                    for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                      for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                        tl_res(rows(l), cols(m)) += tl_lbfi->integrate(m, l);
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
                        tl_gbfi->setPolytope(*trIt, *teIt);
                        const auto& rows = input.getTestFES().getDOFs(d, i);
                        const auto& cols = input.getTrialFES().getDOFs(d, i);
                        for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                          for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                            tl_res(rows(l), cols(m)) += tl_gbfi->integrate(m, l);
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
      static thread_local OperatorType tl_res;
      static thread_local std::unique_ptr<LocalBilinearFormIntegratorBaseType> tl_lbfi;
      static thread_local std::unique_ptr<GlobalBilinearFormIntegratorBaseType> tl_gbfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

   template <class TrialFES, class TestFES>
   thread_local
   Math::Matrix<Real>
   Multithreaded<Math::Matrix<Real>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>::tl_res;

   template <class TrialFES, class TestFES>
   thread_local
   std::unique_ptr<Variational::LocalBilinearFormIntegratorBase<Real>>
   Multithreaded<Math::Matrix<Real>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>::tl_lbfi;

   template <class TrialFES, class TestFES>
   thread_local
   std::unique_ptr<Variational::GlobalBilinearFormIntegratorBase<Real>>
   Multithreaded<Math::Matrix<Real>, Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>::tl_gbfi;

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
                    tl_lfi->setPolytope(*it);
                    const auto& dofs = input.getFES().getDOFs(d, i);
                    for (size_t l = 0; l < static_cast<size_t>(dofs.size()); l++)
                      tl_res(dofs(l)) += tl_lfi->integrate(l);
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
      static thread_local std::unique_ptr<Variational::LinearFormIntegratorBase<ScalarType>> tl_lfi;

      mutable Threads::Mutex m_mutex;
      mutable std::variant<Threads::ThreadPool, std::reference_wrapper<Threads::ThreadPool>> m_pool;
  };

  template <class FES>
  thread_local
  Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>
  Multithreaded<
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
    Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
  ::tl_res;

  template <class FES>
  thread_local
  std::unique_ptr<Variational::LinearFormIntegratorBase<typename FormLanguage::Traits<FES>::ScalarType>>
  Multithreaded<
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
    Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
  ::tl_lfi;
}

#endif

