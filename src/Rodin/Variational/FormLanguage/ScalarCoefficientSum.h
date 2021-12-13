/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SCALARCOEFFSUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SCALARCOEFFSUM_H

#include <memory>
#include <utility>

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   class ScalarCoefficientSum : public Base
   {
      public:
         ScalarCoefficientSum(
               const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

         ScalarCoefficientSum(const ScalarCoefficientSum& other);

         ScalarCoefficientBase& getLHS();

         ScalarCoefficientBase& getRHS();

         virtual ScalarCoefficientSum* copy() const noexcept override
         {
            return new ScalarCoefficientSum(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };

   ScalarCoefficientSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif
