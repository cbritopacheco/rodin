/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"

#include "Multithreaded.h"

namespace Rodin::Assembly
{
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
    : m_threadCount(threadCount),
      m_pool(threadCount)
  {
    assert(threadCount > 0);
  }

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded(const Multithreaded& other)
    : Parent(other),
      m_threadCount(other.m_threadCount),
      m_pool(m_threadCount)
  {}

  Multithreaded<Variational::LinearFormBase<Math::Vector>>
  ::Multithreaded(Multithreaded&& other)
    : Parent(std::move(other)),
      m_threadCount(std::move(other.m_threadCount)),
      m_pool(m_threadCount)
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
}



