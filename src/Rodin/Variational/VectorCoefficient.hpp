/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORCOEFFICIENT_HPP
#define RODIN_VARIATIONAL_VECTORCOEFFICIENT_HPP

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   // ---- std::initializer_list<T> ------------------------------------------
   template <class T>
   constexpr
   VectorCoefficient<std::initializer_list<T>>
   ::VectorCoefficient(std::initializer_list<T> values)
      : m_dimension(values.size())
   {
      m_values.reserve(m_dimension);
      for (auto& v : values)
         m_values.emplace_back(new ScalarCoefficient<T>(v));
   }

   template <class T>
   constexpr
   VectorCoefficient<std::initializer_list<T>>
   ::VectorCoefficient(const VectorCoefficient& other)
      :  m_dimension(other.m_dimension)
   {
      m_values.reserve(m_dimension);
      for (auto& v : other.m_values)
         m_values.emplace_back(v->copy());
   }

   template <class T>
   size_t
   VectorCoefficient<std::initializer_list<T>>
   ::getDimension() const
   {
      return m_dimension;
   }

   template <class T>
   void
   VectorCoefficient<std::initializer_list<T>>
   ::buildMFEMVectorCoefficient()
   {
      m_mfemVectorCoefficient.emplace(m_dimension);
      for (size_t i = 0; i < m_values.size(); i++)
      {
         m_values[i]->buildMFEMCoefficient();
         m_mfemVectorCoefficient->Set(i, &m_values[i]->getMFEMCoefficient(), false);
      }
   }

   template <class T>
   mfem::VectorCoefficient&
   VectorCoefficient<std::initializer_list<T>>
   ::getMFEMVectorCoefficient()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }

   // ---- GridFunction<FEC> -------------------------------------------------
   template <class FEC>
   constexpr
   VectorCoefficient<GridFunction<FEC>>
   ::VectorCoefficient(GridFunction<FEC>& u)
      :  m_dimension(u.getFiniteElementSpace().getDimension()),
         m_u(u)
   {}

   template <class FEC>
   constexpr
   VectorCoefficient<GridFunction<FEC>>
   ::VectorCoefficient(const VectorCoefficient& other)
      :  m_dimension(other.m_dimension),
         m_u(other.m_u)
   {}

   template <class FEC>
   size_t
   VectorCoefficient<GridFunction<FEC>>::getDimension() const
   {
      return m_dimension;
   }

   template <class FEC>
   void
   VectorCoefficient<GridFunction<FEC>>::buildMFEMVectorCoefficient()
   {
      m_mfemVectorCoefficient.emplace(&m_u.getHandle());
   }

   template <class FEC>
   mfem::VectorCoefficient&
   VectorCoefficient<GridFunction<FEC>>::getMFEMVectorCoefficient()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }
}

#endif
