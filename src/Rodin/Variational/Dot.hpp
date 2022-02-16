/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_HPP
#define RODIN_VARIATIONAL_DOT_HPP

#include "MatrixCoefficient.h"

#include "Dot.h"

namespace Rodin::Variational
{
   // ---- Dot<MatrixCoefficientBase> ----------------------------------------
   template <class A, class B>
   constexpr
   Dot<A, B>::Dot(const A& a, const B& b)
      : m_a(a.copy()), m_b(b.copy())
   {}

   template <class A, class B>
   constexpr
   Dot<A, B>::Dot(const Dot& other)
      :  m_a(other.m_a->copy()), m_b(other.m_b->copy())
   {}

   template <class A, class B>
   void Dot<A, B>::build()
   {
      m_a->build();
      m_b->build();

      m_mfemCoefficient.emplace(
            m_a->get(), m_b->get());
   }

   template <class A, class B>
   mfem::Coefficient& Dot<A, B>::get()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
#endif
