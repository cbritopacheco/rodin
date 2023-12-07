/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_MULTITHREADED_H
#define RODIN_ASSEMBLY_MULTITHREADED_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Threads/Mutex.h"
#include "Rodin/Threads/ThreadPool.h"

#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES>
  class Multithreaded<Variational::BilinearForm<TrialFES, TestFES,
        std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES,
          std::vector<Eigen::Triplet<Scalar>>>>
  {
    /**
     * @internal
     */
    static void add(
        std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols)
    {
      assert(rows.size() >= 0);
      assert(cols.size() >= 0);
      assert(in.rows() == rows.size());
      assert(in.cols() == cols.size());
      for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
      {
        for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
        {
          const Scalar s = in(i, j);
          if (s != Scalar(0))
            out.emplace_back(rows(i), cols(j), s);
        }
      }
    }

    public:
      using Parent =
        AssemblyBase<Variational::BilinearForm<TrialFES, TestFES,
          std::vector<Eigen::Triplet<Scalar>>>>;

      using Input = typename Parent::Input;

      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
      {}

      Multithreaded(size_t threadCount)
        : m_threadCount(threadCount),
          m_pool(threadCount)
      {
        assert(threadCount > 0);
      }

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount),
          m_pool(m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount)),
          m_pool(m_threadCount)
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const Input& input) const override
      {
        using TripletVector = std::vector<Eigen::Triplet<Scalar>>;
        TripletVector res;
        res.clear();
        res.reserve(input.testFES.getSize() * std::log(input.trialFES.getSize()));
        auto& threadPool = m_pool;
        for (auto& bfi : input.bfis)
        {
          const auto& attrs = bfi.getAttributes();
          switch (bfi.getRegion())
          {
            case Variational::Integrator::Region::Domain:
            {
              const size_t d = input.mesh.getDimension();
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_triplets.reserve(
                      input.testFES.getSize() * std::log(input.trialFES.getSize()) / m_threadCount);
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getCell(i);
                      const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                      const auto& testDOFs = input.testFES.getDOFs(d, i);
                      tl_bfi->assemble(*it);
                      add(tl_triplets, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                    }
                  }
                  m_mutex.lock();
                  res.insert(res.end(),
                      std::make_move_iterator(tl_triplets.begin()),
                      std::make_move_iterator(tl_triplets.end()));
                  m_mutex.unlock();
                  tl_triplets.clear();
                };
              threadPool.pushLoop(0, input.mesh.getCellCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Faces:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_triplets.reserve(
                      input.testFES.getSize() * std::log(input.trialFES.getSize()) / m_threadCount);
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getFace(i);
                      const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                      const auto& testDOFs = input.testFES.getDOFs(d, i);
                      tl_bfi->assemble(*it);
                      add(tl_triplets, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                    }
                  }
                  m_mutex.lock();
                  res.insert(res.end(),
                      std::make_move_iterator(tl_triplets.begin()),
                      std::make_move_iterator(tl_triplets.end()));
                  m_mutex.unlock();
                  tl_triplets.clear();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Boundary:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_triplets.reserve(
                      input.testFES.getSize() * std::log(input.trialFES.getSize()) / m_threadCount);
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isBoundary(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                        const auto& testDOFs = input.testFES.getDOFs(d, i);
                        tl_bfi->assemble(*it);
                        add(tl_triplets, tl_bfi->getMatrix(), testDOFs, trialDOFs);
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
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Interface:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_triplets.reserve(
                      input.testFES.getSize() * std::log(input.trialFES.getSize()) / m_threadCount);
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isInterface(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                        const auto& testDOFs = input.testFES.getDOFs(d, i);
                        tl_bfi->assemble(*it);
                        add(tl_triplets, tl_bfi->getMatrix(), testDOFs, trialDOFs);
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
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
          }
        }
        return res;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local std::vector<Eigen::Triplet<Scalar>> tl_triplets;
      static thread_local std::unique_ptr<Variational::BilinearFormIntegratorBase> tl_bfi;

      const size_t m_threadCount;
      mutable Threads::Mutex m_mutex;
      mutable Threads::ThreadPool m_pool;
  };

  template <class TrialFES, class TestFES>
  thread_local
  std::vector<Eigen::Triplet<Scalar>>
  Multithreaded<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>::tl_triplets;

  template <class TrialFES, class TestFES>
  thread_local
  std::unique_ptr<Variational::BilinearFormIntegratorBase>
  Multithreaded<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>::tl_bfi;

  /**
   * @brief Multithreaded assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Multithreaded<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>;
      using Input = typename Parent::Input;
      using OperatorType = Math::SparseMatrix;

      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
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
      OperatorType execute(const Input& input) const override
      {
        const auto triplets = m_assembly.execute({ input.mesh, input.trialFES, input.testFES, input.bfis });
        OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      Multithreaded<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> m_assembly;
  };

  template <class TrialFES, class TestFES>
  class Multithreaded<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
  {
    /**
     * @internal
     */
    static void add(
        Math::Matrix& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols)
    {
      assert(rows.size() >= 0);
      assert(cols.size() >= 0);
      assert(in.rows() == rows.size());
      assert(in.cols() == cols.size());
      out(rows, cols).noalias() += in;
    }

    public:
      using Parent =
        AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>;
      using Input = typename Parent::Input;
      using OperatorType = Math::Matrix;

      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
      {}

      Multithreaded(size_t threadCount)
        : m_threadCount(threadCount),
          m_pool(threadCount)
      {
        assert(threadCount > 0);
      }

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount),
          m_pool(m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount)),
          m_pool(m_threadCount)
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const Input& input) const override
      {
        Math::Matrix res(input.testFES.getSize(), input.trialFES.getSize());
        res.setZero();
        auto& threadPool = m_pool;
        for (auto& bfi : input.bfis)
        {
          const auto& attrs = bfi.getAttributes();
          switch (bfi.getRegion())
          {
            case Variational::Integrator::Region::Domain:
            {
              const size_t d = input.mesh.getDimension();
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_res.resize(input.testFES.getSize(), input.trialFES.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getCell(i);
                      const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                      const auto& testDOFs = input.testFES.getDOFs(d, i);
                      tl_bfi->assemble(*it);
                      add(tl_res, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getCellCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Faces:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_res.resize(input.testFES.getSize(), input.trialFES.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getFace(i);
                      const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                      const auto& testDOFs = input.testFES.getDOFs(d, i);
                      tl_bfi->assemble(*it);
                      add(tl_res, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Boundary:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_res.resize(input.testFES.getSize(), input.trialFES.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isBoundary(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                        const auto& testDOFs = input.testFES.getDOFs(d, i);
                        tl_bfi->assemble(*it);
                        add(tl_res, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                      }
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Interface:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_bfi.reset(bfi.copy());
                  tl_res.resize(input.testFES.getSize(), input.trialFES.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isInterface(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                        const auto& testDOFs = input.testFES.getDOFs(d, i);
                        tl_bfi->assemble(*it);
                        add(tl_res, tl_bfi->getMatrix(), testDOFs, trialDOFs);
                      }
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
          }
        }
        return res;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local Math::Matrix tl_res;
      static thread_local std::unique_ptr<Variational::BilinearFormIntegratorBase> tl_bfi;

      const size_t m_threadCount;
      mutable Threads::Mutex m_mutex;
      mutable Threads::ThreadPool m_pool;
  };

   template <class TrialFES, class TestFES>
   thread_local
   Math::Matrix
   Multithreaded<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>::tl_res;

   template <class TrialFES, class TestFES>
   thread_local
   std::unique_ptr<Variational::BilinearFormIntegratorBase>
   Multithreaded<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>::tl_bfi;

  /**
   * @brief %Multithreaded assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <class FES>
  class Multithreaded<Variational::LinearForm<FES, Math::Vector>>
    : public AssemblyBase<Variational::LinearForm<FES, Math::Vector>>
  {

    static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s)
    {
      assert(in.size() == s.size());
      out(s).noalias() += in;
    }

    public:
      using Parent = AssemblyBase<Variational::LinearForm<FES, Math::Vector>>;
      using Input = typename Parent::Input;
      using VectorType = Math::Vector;

      Multithreaded()
        : Multithreaded(std::thread::hardware_concurrency())
      {}

      Multithreaded(size_t threadCount)
        : m_threadCount(threadCount),
          m_pool(threadCount)
      {
        assert(threadCount > 0);
      }

      Multithreaded(const Multithreaded& other)
        : Parent(other),
          m_threadCount(other.m_threadCount),
          m_pool(m_threadCount)
      {}

      Multithreaded(Multithreaded&& other)
        : Parent(std::move(other)),
          m_threadCount(std::move(other.m_threadCount)),
          m_pool(m_threadCount)
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const Input& input) const override
      {
        VectorType res(input.fes.getSize());
        res.setZero();
        auto& threadPool = m_pool;
        for (auto& lfi : input.lfis)
        {
          const auto& attrs = lfi.getAttributes();
          switch (lfi.getRegion())
          {
            case Variational::Integrator::Region::Domain:
            {
              const size_t d = input.mesh.getDimension();
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_lfi.reset(lfi.copy());
                  tl_res.resize(input.fes.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getCell(i);
                      const auto& dofs = input.fes.getDOFs(d, i);
                      tl_lfi->assemble(*it);
                      add(tl_res, tl_lfi->getVector(), dofs);
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getCellCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Faces:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_lfi.reset(lfi.copy());
                  tl_res.resize(input.fes.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                    {
                      const auto it = input.mesh.getFace(i);
                      const auto& dofs = input.fes.getDOFs(d, i);
                      tl_lfi->assemble(*it);
                      add(tl_res, tl_lfi->getVector(), dofs);
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Boundary:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_lfi.reset(lfi.copy());
                  tl_res.resize(input.fes.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isBoundary(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& dofs = input.fes.getDOFs(d, i);
                        tl_lfi->assemble(*it);
                        add(tl_res, tl_lfi->getVector(), dofs);
                      }
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
            case Variational::Integrator::Region::Interface:
            {
              const size_t d = input.mesh.getDimension() - 1;
              auto loop =
                [&](const Index start, const Index end)
                {
                  tl_lfi.reset(lfi.copy());
                  tl_res.resize(input.fes.getSize());
                  tl_res.setZero();
                  for (Index i = start; i < end; ++i)
                  {
                    if (input.mesh.isInterface(i))
                    {
                      if (attrs.size() == 0 || attrs.count(input.mesh.getAttribute(d, i)))
                      {
                        const auto it = input.mesh.getFace(i);
                        const auto& dofs = input.fes.getDOFs(d, i);
                        tl_lfi->assemble(*it);
                        add(tl_res, tl_lfi->getVector(), dofs);
                      }
                    }
                  }
                  m_mutex.lock();
                  res += tl_res;
                  m_mutex.unlock();
                };
              threadPool.pushLoop(0, input.mesh.getFaceCount(), loop);
              threadPool.waitForTasks();
              break;
            }
          }
        }
        return res;
      }

      Multithreaded* copy() const noexcept override
      {
        return new Multithreaded(*this);
      }

    private:
      static thread_local VectorType tl_res;
      static thread_local std::unique_ptr<Variational::LinearFormIntegratorBase> tl_lfi;

      const size_t m_threadCount;
      mutable Threads::Mutex m_mutex;
      mutable Threads::ThreadPool m_pool;
  };

  template <class FES>
  thread_local
  Math::Vector Multithreaded<Variational::LinearForm<FES, Math::Vector>>::tl_res;

  template <class FES>
  thread_local
  std::unique_ptr<Variational::LinearFormIntegratorBase>
  Multithreaded<Variational::LinearForm<FES, Math::Vector>>::tl_lfi;
}

#endif

