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

#include "FormLanguage/ScalarCoefficientSum.h"
#include "FormLanguage/ScalarCoefficientUnaryMinus.h"

namespace Rodin::Variational
{
   // ---- T (Arithmetic type) -----------------------------------------------
   // ------------------------------------------------------------------------
   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
   ::ScalarCoefficient(const T& x)
      : m_x(x)
   {}

   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_x(other.m_x),
         m_mfemCoefficient(other.m_mfemCoefficient)
   {}

   template <class T>
   void
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::buildMFEMCoefficient()
   {
      m_mfemCoefficient.emplace(m_x);
   }

   template <class T>
   mfem::Coefficient&
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }

   // ---- GridFunction<FEC> -------------------------------------------------
   // ------------------------------------------------------------------------
   template <class FEC>
   ScalarCoefficient<GridFunction<FEC>>
   ::ScalarCoefficient(GridFunction<FEC>& u)
      : m_u(u)
   {
      if (u.getFiniteElementSpace().getDimension() != 1)
      {
         (Alert::Exception() << "ScalarCoefficient can only be initialized "
                             << "with a scalar valued GridFunction").raise();

      }
   }

   template <class FEC>
   ScalarCoefficient<GridFunction<FEC>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_u(other.m_u),
         m_mfemCoefficient(other.m_mfemCoefficient)
   {}


   template <class FEC>
   void
   ScalarCoefficient<GridFunction<FEC>>::buildMFEMCoefficient()
   {
      m_mfemCoefficient.emplace(m_u);
   }

   template <class FEC>
   mfem::Coefficient& ScalarCoefficient<GridFunction<FEC>>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}

#endif
