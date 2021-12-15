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
   // ------------------------------------------------------------------------
   template <class A>
   Dot<A,
      std::enable_if_t<std::is_base_of_v<MatrixCoefficientBase, A>>>
   ::Dot(const A& a, const A& b)
      : m_a(a.copy()), m_b(b.copy())
   {}

   template <class A>
   Dot<A,
      std::enable_if_t<std::is_base_of_v<MatrixCoefficientBase, A>>>
   ::Dot(const Dot& other)
      : m_a(other.m_a->copy()), m_b(other.m_b->copy())
   {}

   template <class A>
   void
   Dot<A,
      std::enable_if_t<std::is_base_of_v<MatrixCoefficientBase, A>>>
   ::buildMFEMCoefficient()
   {
      m_a->buildMFEMMatrixCoefficient();
      m_b->buildMFEMMatrixCoefficient();

      m_mfemCoefficient.emplace(
            m_a->getMFEMMatrixCoefficient(), m_b->getMFEMMatrixCoefficient());
   }

   template <class A>
   mfem::Coefficient&
   Dot<A,
      std::enable_if_t<std::is_base_of_v<MatrixCoefficientBase, A>>>
   ::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}
#endif
