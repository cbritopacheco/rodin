/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARCOEFFICIENT_HPP
#define RODIN_VARIATIONAL_SCALARCOEFFICIENT_HPP

#include "Rodin/Alert.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   // ---- T (Arithmetic type) -----------------------------------------------
   // ------------------------------------------------------------------------
   template <class T>
   void
   ScalarCoefficient<T>::buildMFEMCoefficient()
   {
      m_mfemCoefficient.emplace(m_x);
   }

   template <class T>
   mfem::Coefficient&
   ScalarCoefficient<T>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   // ---- GridFunction<FEC> -------------------------------------------------
   // ------------------------------------------------------------------------
   template <class FEC>
   constexpr
   ScalarCoefficient<GridFunction<FEC>>
   ::ScalarCoefficient(GridFunction<FEC>& u)
      : m_u(u)
   {
      assert(u.getFiniteElementSpace().getDimension() == 1);
   }

   template <class FEC>
   constexpr
   ScalarCoefficient<GridFunction<FEC>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_u(other.m_u)
   {}


   template <class FEC>
   void
   ScalarCoefficient<GridFunction<FEC>>::buildMFEMCoefficient()
   {
      m_mfemCoefficient.emplace(&m_u.getHandle());
   }

   template <class FEC>
   mfem::Coefficient& ScalarCoefficient<GridFunction<FEC>>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}

#endif
