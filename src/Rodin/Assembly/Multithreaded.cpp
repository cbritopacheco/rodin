/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>

#include "BS_thread_pool.hpp"

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Multithreaded.h"

namespace Rodin::Assembly
{
  Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::Multithreaded()
    : Multithreaded(std::thread::hardware_concurrency())
  {}

  Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::Multithreaded(size_t threadCount)
    : m_threadCount(threadCount)
  {
    assert(threadCount > 0);
  }

  Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::Multithreaded(const Multithreaded& other)
    : Parent(other),
      m_threadCount(other.m_threadCount)
  {}

  Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::Multithreaded(Multithreaded&& other)
    : Parent(std::move(other)),
      m_threadCount(std::move(other.m_threadCount))
  {}

  Math::SparseMatrix
  Multithreaded<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::execute(const BilinearAssemblyInput& input) const
  {
    Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>> assembly(m_threadCount);
    const auto triplets = assembly.execute(input);
    OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
    res.setFromTriplets(triplets.begin(), triplets.end());
    return res;
  }

  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::Multithreaded()
    : Multithreaded(std::thread::hardware_concurrency())
  {}

  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::Multithreaded(size_t threadCount)
    : m_threadCount(threadCount)
  {
    assert(threadCount > 0);
  }

  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::Multithreaded(const Multithreaded& other)
    : Parent(other),
      m_threadCount(other.m_threadCount)
  {}

  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::Multithreaded(Multithreaded&& other)
    : Parent(std::move(other)),
      m_threadCount(std::move(other.m_threadCount))
  {}

  void
  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::add(std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
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

  thread_local
  std::vector<Eigen::Triplet<Scalar>>
  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>::tl_triplets;

  thread_local
  std::unique_ptr<Variational::BilinearFormIntegratorBase>
  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>::tl_bfi;

  std::vector<Eigen::Triplet<Scalar>>
  Multithreaded<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::execute(const BilinearAssemblyInput& input) const
  {
    using TripletVector = std::vector<Eigen::Triplet<Scalar>>;
    TripletVector res;
    res.clear();
    res.reserve(input.testFES.getSize() * std::log(input.trialFES.getSize()));
    BS::thread_pool threadPool(m_threadCount);
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
          threadPool.push_loop(0, input.mesh.getCellCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
          break;
        }
      }
    }
    return res;
  }

  thread_local
  Multithreaded<Variational::LinearFormBase<Math::Vector>>::VectorType
  Multithreaded<Variational::LinearFormBase<Math::Vector>>::tl_res;

  thread_local
  std::unique_ptr<Variational::LinearFormIntegratorBase>
  Multithreaded<Variational::LinearFormBase<Math::Vector>>::tl_lfi;

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded()
    : Multithreaded(std::thread::hardware_concurrency())
  {}

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded(size_t threadCount)
    : m_threadCount(threadCount)
  {
    assert(threadCount > 0);
  }

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded(const Multithreaded& other)
    : Parent(other),
      m_threadCount(other.m_threadCount)
  {}

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded(Multithreaded&& other)
    : Parent(std::move(other)),
      m_threadCount(std::move(other.m_threadCount))
  {}

  void
  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::add(Math::Vector& out, const Math::Vector& in, const IndexArray& s)
  {
    assert(in.size() == s.size());
    size_t i = 0;
    for (const auto& global : s)
      out.coeffRef(global) += in.coeff(i++);
  }

  Math::Vector
  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::execute(const Input& input) const
  {
    BS::thread_pool threadPool(m_threadCount);
    VectorType res(input.fes.getSize());
    res.setZero();
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
          threadPool.push_loop(0, input.mesh.getCellCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
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
          threadPool.push_loop(0, input.mesh.getFaceCount(), loop);
          threadPool.wait_for_tasks();
          break;
        }
      }
    }
    return res;
  }
}



