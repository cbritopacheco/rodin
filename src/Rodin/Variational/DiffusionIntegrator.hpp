/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_HPP
#define RODIN_VARIATIONAL_DIFFUSIONINTEGRATOR_HPP

#include "Problem.h"

#include "DiffusionIntegrator.h"

namespace Rodin::Variational
{
   template <class T>
   DiffusionIntegrator<T>::DiffusionIntegrator(const T& lambda)
      : m_lambda(ScalarCoefficient<T>(lambda))
   {}

   template <class T>
   DiffusionIntegrator<T>::DiffusionIntegrator(const DiffusionIntegrator& other)
      : m_lambda(other.m_lambda), m_bf(other.m_bf)
   {}

   template <class T>
   DiffusionIntegrator<T>&
   DiffusionIntegrator<T>::setBilinearForm(BilinearFormBase& bf)
   {
      m_bf.emplace(bf);
      return *this;
   }

   template <class T>
   void DiffusionIntegrator<T>::eval()
   {
      assert(m_bf);
      m_lambda.eval();
      m_bf->get()
           .getHandle()
           .AddDomainIntegrator(
                 new mfem::DiffusionIntegrator(m_lambda.coeff()));
   }

   template <class T>
   DiffusionIntegrator<T>& DiffusionIntegrator<T>::toggleSign()
   {
      m_lambda.toggleSign();
      return *this;
   }

   template <class T>
   template <class ... Args>
   DiffusionIntegrator<T>*
   DiffusionIntegrator<T>::create(Args&&... args) noexcept
   {
      return new DiffusionIntegrator(std::forward<Args>(args)...);
   }

   template <class T>
   DiffusionIntegrator<T>*
   DiffusionIntegrator<T>::copy() const noexcept
   {
      return new DiffusionIntegrator(*this);
   }
}

#endif
