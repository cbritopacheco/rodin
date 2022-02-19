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

#include "Rodin/Variational/ScalarCoefficient.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class ScalarSum : public ScalarCoefficientBase
   {
      public:
         ScalarSum(
               const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

         ScalarSum(const ScalarSum& other);

         ScalarCoefficientBase& getLHS();

         ScalarCoefficientBase& getRHS();

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         virtual ScalarSum* copy() const noexcept override
         {
            return new ScalarSum(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };

   ScalarSum
      operator+(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif
