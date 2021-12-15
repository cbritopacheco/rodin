/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_COEFFUNARYMINUS_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_COEFFUNARYMINUS_H

#include <memory>
#include <optional>

#include "Rodin/Variational/ScalarCoefficient.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class ScalarCoefficientUnaryMinus : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficientUnaryMinus(const ScalarCoefficientBase& s);

         ScalarCoefficientUnaryMinus(const ScalarCoefficientUnaryMinus& other);

         ScalarCoefficientUnaryMinus(ScalarCoefficientUnaryMinus&&) = default;

         ScalarCoefficientBase& getScalarCoefficient();

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         virtual ScalarCoefficientUnaryMinus* copy() const noexcept override
         {
            return new ScalarCoefficientUnaryMinus(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_s;
         std::optional<mfem::SumCoefficient> m_mfemCoefficient;
   };

   ScalarCoefficientUnaryMinus operator-(const ScalarCoefficientBase& s);

   ScalarSum
      operator-(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif
