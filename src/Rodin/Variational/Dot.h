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
   namespace Internal
   {
      /**
       * @internal
       */
      class MatrixDotProductCoefficient : public mfem::Coefficient
      {
         public:
            MatrixDotProductCoefficient(
                  mfem::MatrixCoefficient& a, mfem::MatrixCoefficient& b)
               : m_a(a), m_b(b)
            {}

            virtual double Eval(
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               mfem::DenseMatrix ma, mb;
               m_a.Eval(ma, T, ip);
               m_b.Eval(mb, T, ip);
               return ma * mb;
            }

         private:
            mfem::MatrixCoefficient& m_a;
            mfem::MatrixCoefficient& m_b;
      };
   }

   template <class A, class B>
   Dot(const A&, const B&) -> Dot<A, B>;

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
   template <class A, class B>
   class Dot<A, B> : public ScalarCoefficientBase
   {
      static_assert(
            std::is_base_of_v<MatrixCoefficientBase, A> &&
            std::is_base_of_v<MatrixCoefficientBase, B>
            );

      public:
         /**
          * @brief Constructs the Dot product between two given matrices.
          * @param[in] a Derived instance of MatrixCoefficientBase
          * @param[in] b Derived instance of MatrixCoefficientBase
          */
         constexpr
         Dot(const A& a, const B& b);

         constexpr
         Dot(const Dot& other);

         void build() override;

         mfem::Coefficient& get() override;

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<MatrixCoefficientBase> m_a, m_b;
         std::optional<Internal::MatrixDotProductCoefficient> m_mfemCoefficient;
   };
}

#include "Dot.hpp"

#endif
