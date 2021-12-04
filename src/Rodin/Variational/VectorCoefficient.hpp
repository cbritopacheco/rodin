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
   VectorCoefficient<Values...>::VectorCoefficient(const Values&... values)
      :  m_values(values...),
         m_sign(false),
         m_coeff(sizeof...(Values))
   {}

   template <class ... Values>
   VectorCoefficient<Values...>::VectorCoefficient(const VectorCoefficient& other)
      : m_values(other.m_values), m_sign(other.m_sign), m_coeff(other.m_coeff)
   {}

   template <class ... Values>
   void VectorCoefficient<Values...>::eval()
   {
      setCoefficients(m_values);
   }

   template <class ... Values>
   mfem::VectorCoefficient& VectorCoefficient<Values...>::coeff()
   {
      return m_coeff;
   }

   template <class ... Values>
   template<std::size_t I, typename... Tp>
   typename std::enable_if_t<I == sizeof...(Tp)>
   VectorCoefficient<Values...>::setCoefficients(std::tuple<Tp...>&)
   {}

   template <class ... Values>
   template<std::size_t I, typename... Tp>
   typename std::enable_if_t<I < sizeof...(Tp)>
   VectorCoefficient<Values...>::setCoefficients(std::tuple<Tp...>& t)
   {
      auto& coeff = m_coeffs.emplace_back(new ScalarCoefficient(std::get<I>(t)));
      if (m_sign)
         coeff->toggleSign().eval();
      else
         coeff->eval();
      m_coeff.Set(I, &coeff->coeff(), false);
      setCoefficients<I + 1, Tp...>(t);
   }

   template <class ... Values>
   VectorCoefficient<Values...>& VectorCoefficient<Values...>::toggleSign()
   {
      m_sign = !m_sign;
      return *this;
   }

   template <class ... Values>
   template <class ... Args>
   VectorCoefficient<Values...>*
   VectorCoefficient<Values...>::create(Args&&... args) noexcept
   {
      return new VectorCoefficient(std::forward<Args>(args)...);
   }

   template <class ... Values>
   VectorCoefficient<Values...>*
   VectorCoefficient<Values...>::copy() const noexcept
   {
      return new VectorCoefficient(*this);
   }

}

#endif
