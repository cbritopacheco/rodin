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
   // ---- ... Values (Variadic vector) --------------------------------------
   template <class ... Values>
   VectorCoefficient::VectorCoefficient(const Values&... values)
      :  m_dimension(sizeof...(Values)),
         m_mfemVectorArrayCoefficient(m_dimension)
   {
      m_values.reserve(sizeof...(Values));
      makeCoefficientsFromTuple(std::forward_as_tuple(values...));
   }

   template<std::size_t I, class ... Tp>
   typename std::enable_if_t<I == sizeof...(Tp)>
   VectorCoefficient::makeCoefficientsFromTuple(const std::tuple<Tp...>&)
   {}

   template<std::size_t I, class ... Tp>
   typename std::enable_if_t<I < sizeof...(Tp)>
   VectorCoefficient::makeCoefficientsFromTuple(const std::tuple<Tp...>& t)
   {
      auto& s = m_values.emplace_back(new ScalarCoefficient(std::get<I>(t)));
      s->buildMFEMCoefficient();
      m_mfemVectorArrayCoefficient.Set(I, &s->getMFEMCoefficient(), false);
      makeCoefficientsFromTuple<I + 1, Tp...>(t);
   }
}

#endif
