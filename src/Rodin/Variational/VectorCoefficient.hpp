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
   template <class ... Values>
   size_t
   VectorCoefficient<Values...>
   ::getDimension() const
   {
      return m_dimension;
   }

   template <class ... Values>
   void
   VectorCoefficient<Values...>
   ::build()
   {
      m_mfemCoefficients.reserve(m_dimension);
      m_mfemVectorCoefficient.emplace(m_dimension);

      makeCoefficientsFromTuple(m_values);

      for (size_t i = 0; i < m_dimension; i++)
      {
         m_mfemCoefficients[i]->build();
         m_mfemVectorCoefficient->Set(
               i, &m_mfemCoefficients[i]->get(), false);
      }
   }

   template <class ... Values>
   mfem::VectorCoefficient&
   VectorCoefficient<Values...>
   ::get()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }

   template <class ... Values>
   template<std::size_t I, class ... Tp>
   typename std::enable_if_t<I == sizeof...(Tp)>
   VectorCoefficient<Values...>
   ::makeCoefficientsFromTuple(const std::tuple<Tp...>&)
   {}

   template <class ... Values>
   template<std::size_t I, class ... Tp>
   typename std::enable_if_t<I < sizeof...(Tp)>
   VectorCoefficient<Values...>
   ::makeCoefficientsFromTuple(const std::tuple<Tp...>& t)
   {
      assert(m_mfemVectorCoefficient);
      auto& s = m_mfemCoefficients.emplace_back(new ScalarCoefficient(std::get<I>(t)));
      s->build();
      m_mfemVectorCoefficient->Set(I, &s->get(), false);
      makeCoefficientsFromTuple<I + 1, Tp...>(t);
   }

   // ---- GridFunction<FEC> -------------------------------------------------
   template <class FEC>
   constexpr
   VectorCoefficient<GridFunction<FEC>>
   ::VectorCoefficient(GridFunction<FEC>& u)
      :  m_dimension(u.getFiniteElementSpace().getRangeDimension()),
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
   VectorCoefficient<GridFunction<FEC>>::build()
   {
      m_mfemVectorCoefficient.emplace(&m_u.getHandle());
   }

   template <class FEC>
   mfem::VectorCoefficient&
   VectorCoefficient<GridFunction<FEC>>::get()
   {
      assert(m_mfemVectorCoefficient);
      return *m_mfemVectorCoefficient;
   }
}

#endif
