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
      m_mfemCoefficients.emplace_back(new ScalarCoefficient(std::get<I>(t)));
      makeCoefficientsFromTuple<I + 1, Tp...>(t);
   }

   // ---- GridFunction<FEC> -------------------------------------------------
   template <class FEC>
   constexpr
   VectorCoefficient<GridFunction<FEC>>
   ::VectorCoefficient(GridFunction<FEC>& u)
      :  m_dimension(u.getFiniteElementSpace().getRangeDimension()),
         m_u(u),
         m_mfemVectorCoefficient(&u.getHandle())
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
}

#endif
