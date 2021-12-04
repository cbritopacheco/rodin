/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORCOEFFICIENT_H
#define RODIN_VARIATIONAL_VECTORCOEFFICIENT_H

#include <memory>
#include <tuple>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/RodinBase.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   template <class ... Values>
   class VectorCoefficient : public FormLanguage::RodinBase
   {
      public:
         VectorCoefficient(const Values&... values);

         /**
          * Copy constructor
          *
          * @param other Other coefficient to copy
          */
         VectorCoefficient(const VectorCoefficient& other);

         /**
          * Evaluates the coefficient. After this you may obtain the value with
          * the coeff() function.
          *
          * @see coeff()
          */
         void eval() override;

         size_t dimension() const
         {
            return sizeof...(Values);
         }

         mfem::VectorCoefficient& coeff();

         /**
          * Toggles the sign of the coefficient
          *
          * @returns Reference to self (for method chaining)
          */
         VectorCoefficient& toggleSign();

         /**
          * @param args Arguments to pass to the Coeff constructor
          * @returns Non-owning pointer to the new Coeff object
          */
         template <class ... Args>
         static VectorCoefficient* create(Args&&... args) noexcept;

         /**
          * @returns Non-owning pointer to the copy of the Coeff object
          */
         virtual VectorCoefficient* copy() const noexcept override;

      private:
         template<std::size_t I = 0, typename... Tp>
         typename std::enable_if_t<I == sizeof...(Tp)>
         setCoefficients(std::tuple<Tp...>& t);

         template<std::size_t I = 0, typename... Tp>
         typename std::enable_if_t<I < sizeof...(Tp)>
         setCoefficients(std::tuple<Tp...>& t);

         bool m_sign;
         std::tuple<Values...> m_values;
         std::vector<std::unique_ptr<ScalarCoefficientBase>> m_coeffs;
         mfem::VectorArrayCoefficient m_coeff;
   };
}

#include "VectorCoefficient.hpp"

#endif
