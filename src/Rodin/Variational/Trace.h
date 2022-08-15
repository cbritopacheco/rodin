/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "ScalarFunction.h"
#include "BasisOperator.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the trace @f$ \mathrm{tr}(A) @f$ of some square matrix
    * @f$ A @f$.
    *
    * The trace of an @f$ n \times n @f$ matrix @f$ A @f$ is defined as
    * @f[
    *    \mathrm{tr}(A) = \sum_{i = 1}^n A_{ii}
    * @f]
    */
   template <>
   class Trace<FunctionBase> : public ScalarFunctionBase
   {
      public:
         /**
          * @brief Constructs the Trace of the given matrix
          * @param[in] m Square matrix
          */
         Trace(const FunctionBase& m);

         Trace(const Trace& other);

         Trace(Trace&& other);

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         Trace& traceOf(const std::set<int>& attrs) override
         {
            ScalarFunctionBase::traceOf(attrs);
            m_matrix->traceOf(attrs);
            return *this;
         }

         Trace* copy() const noexcept override
         {
            return new Trace(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_matrix;
   };
   Trace(const FunctionBase&) -> Trace<FunctionBase>;


   /**
    * @f[
    *    \text{tr} \ A(u)
    * @f]
    * with @f$ A(u) \in \mathbb{R}^{p \times p} @f$.
    */
   template <>
   class Trace<BasisOperator> : public DenseBasisOperator
   {
      public:
         Trace(std::unique_ptr<BasisOperator> op)
            : DenseBasisOperator(1, 1, op->getDOFs()),
              m_op(std::move(op))
         {
            if (getOperator().isDense())
            {
               const auto& dense = static_cast<const DenseBasisOperator&>(getOperator());
               for (int i = 0; i < getDOFs(); i++)
                  operator()(0, 0, i) = dense(i).Trace();
            }
            else if (getOperator().isSparse())
            {
               const auto& sparse = static_cast<const SparseBasisOperator&>(getOperator());
               for (int i = 0; i < getDOFs(); i++)
               {
                  mfem::Vector diag;
                  sparse(i).GetDiag(diag);
                  operator()(0, 0, i) = diag.Sum();
               }
            }
            else
            {
               assert(false);
            }
         }

         BasisOperator& getOperator()
         {
            return *m_op;
         }

      private:
         std::unique_ptr<BasisOperator> m_op;
   };
   Trace(std::unique_ptr<BasisOperator>) -> Trace<BasisOperator>;
}

#endif
