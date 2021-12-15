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
      /**
       * @internal
       */
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

   /**
    * @brief Represents the trace @f$ \mathrm{tr}(A) @f$ of some square matrix
    * @f$ A @f$.
    *
    * The trace of an @f$ n \times n @f$ matrix @f$ A @f$ is defined as
    * @f[
    *    \mathrm{tr}(A) = \sum_{i = 1}^n A_{ii}
    * @f]
    */
   class Trace : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Trace of the given matrix
          * @param[in] m Square matrix
          */
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
