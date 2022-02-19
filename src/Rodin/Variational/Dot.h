/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

#include <utility>
#include <optional>

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the double-dot product between two matrices.
    *
    * For two given @f$ n \times m @f$ matrices @f$ A @f$ and @f$ B @f$, their
    * double-dot product @f$ A \colon B @f$ is defined as
    * @f[
    *    A \colon B = \sum_{i = 1}^n \sum_{j = 1}^m A_{ij} B_{ij} = \mathrm{tr}(B^T A)
    * @f]
    *
    * @tparam A Derived type from MatrixCoefficientBase
    * @tparam B Derived type from MatrixCoefficientBase
    */
   template <>
   class Dot<MatrixCoefficientBase> : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Dot product between two given matrices.
          * @param[in] a Derived instance of MatrixCoefficientBase
          * @param[in] b Derived instance of MatrixCoefficientBase
          */
         Dot(const MatrixCoefficientBase& a, const MatrixCoefficientBase& b);

         Dot(const Dot& other);

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<MatrixCoefficientBase> m_a, m_b;
   };

   Dot(const MatrixCoefficientBase&, const MatrixCoefficientBase&) -> Dot<MatrixCoefficientBase>;
}

#endif
