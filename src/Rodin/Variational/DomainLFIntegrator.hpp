/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_HPP
#define RODIN_VARIATIONAL_DOMAINLFINTEGRATOR_HPP

#include "Problem.h"

#include "DomainLFIntegrator.h"

namespace Rodin::Variational
{
   template <class T>
   DomainLFIntegrator<T>::DomainLFIntegrator(const T& f)
      : m_f(ScalarCoefficient<T>(f))
   {
      // Make the coefficient negative so that we have -(f, v) on the right
      // left hand side of the variational problem
      m_f.toggleSign();
   }

   template <class T>
   DomainLFIntegrator<T>::DomainLFIntegrator(const DomainLFIntegrator& other)
      : m_f(other.m_f), m_lf(other.m_lf)
   {}

   template <class T>
   DomainLFIntegrator<T>& DomainLFIntegrator<T>::setLinearForm(LinearFormBase& lf)
   {
      m_lf.emplace(lf);
      return *this;
   }

   template <class T>
   void DomainLFIntegrator<T>::eval()
   {
      m_f.buildMFEMCoefficient();
      m_lf->get()
         .getHandle()
         .AddDomainIntegrator(
               new mfem::DomainLFIntegrator(
                  m_f.getMFEMCoefficient()));
   }

   template <class T>
   DomainLFIntegrator<T>& DomainLFIntegrator<T>::toggleSign()
   {
      m_f.toggleSign();
      return *this;
   }

   template <class T>
   template <class ... Args>
   DomainLFIntegrator<T>*
   DomainLFIntegrator<T>::create(Args&&... args)
   noexcept
   {
      return new DomainLFIntegrator(std::forward<Args>(args)...);
   }

   template <class T>
   DomainLFIntegrator<T>* DomainLFIntegrator<T>::copy()
   const noexcept
   {
      return new DomainLFIntegrator(*this);
   }
}

#endif
