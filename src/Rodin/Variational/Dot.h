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

   template <class A, class B, class Enable>
   class Dot
   {
      public:
         Dot(const A&, const B&);
   };

   template <class A, class B>
   class Dot<A, B,
         std::enable_if_t<
            std::is_base_of_v<MatrixCoefficientBase, A> &&
            std::is_base_of_v<MatrixCoefficientBase, B>>>
      : public ScalarCoefficientBase
   {
      public:
         Dot(const A& a, const B& b);

         Dot(const Dot& other);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

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
