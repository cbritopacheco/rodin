/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class TraceCoefficient : public mfem::Coefficient
      {
         public:
            TraceCoefficient(mfem::MatrixCoefficient& matrix)
               : m_matrix(matrix)
            {}

            virtual double Eval(
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               mfem::DenseMatrix mat;
               m_matrix.Eval(mat, T, ip);
               return mat.Trace();
            }

         private:
            mfem::MatrixCoefficient& m_matrix;
      };
   }

   class Trace : public ScalarCoefficientBase
   {
      public:
         Trace(const MatrixCoefficientBase& m);

         Trace(const Trace& other);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         Trace* copy() const noexcept override
         {
            return new Trace(*this);
         }
      private:
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
         std::optional<Internal::TraceCoefficient> m_mfemCoefficient;
   };
}

#endif
