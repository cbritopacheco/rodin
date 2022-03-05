/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_MULTIPLICATION_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_MULTIPLICATION_H

#include <memory>
#include <type_traits>

#include "Rodin/Alert.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ScalarCoefficient.h"
#include "Rodin/Variational/VectorCoefficient.h"
#include "Rodin/Variational/MatrixCoefficient.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/TrialFunction.h"
#include "FormLanguage/Base.h"
#include "ForwardDecls.h"

namespace Rodin::Variational
{
   /**
    * @brief Multiplication of two ScalarCoefficientBase instances.
    */
   template <>
   class Mult<ScalarCoefficientBase, ScalarCoefficientBase>
      : public ScalarCoefficientBase
   {
      public:
         Mult(double lhs, const ScalarCoefficientBase& rhs)
            : Mult(ScalarCoefficient(lhs), rhs)
         {}

         Mult(const ScalarCoefficientBase& rhs, double lhs)
            : Mult(ScalarCoefficient(lhs), rhs)
         {}

         Mult(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  ScalarCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  ScalarCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ScalarCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            return getLHS().getValue(trans, ip) * getRHS().getValue(trans, ip);
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };
   Mult(const ScalarCoefficientBase&, const ScalarCoefficientBase&)
      -> Mult<ScalarCoefficientBase, ScalarCoefficientBase>;
   Mult(double, const ScalarCoefficientBase&)
      -> Mult<ScalarCoefficientBase, ScalarCoefficientBase>;
   Mult(const ScalarCoefficientBase&, double)
      -> Mult<ScalarCoefficientBase, ScalarCoefficientBase>;
   // Mult<ScalarCoefficientBase, ScalarCoefficientBase>
   // operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   // template <class T>
   // std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ScalarCoefficientBase>>
   // operator*(T lhs, const ScalarCoefficientBase& rhs)
   // {
   //    return Mult(ScalarCoefficient(lhs), rhs);
   // }

   // template <class T>
   // std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ScalarCoefficientBase>>
   // operator*(const ScalarCoefficientBase& rhs, T lhs)
   // {
   //    return Mult(rhs, ScalarCoefficient(lhs));
   // }

   /**
    * @brief Multiplication of ScalarCoefficientBase and MatrixCoefficientBase.
    */
   template <>
   class Mult<ScalarCoefficientBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Mult(const MatrixCoefficientBase& rhs, double lhs)
            : Mult(ScalarCoefficient(lhs), rhs)
         {}

         Mult(double lhs, const MatrixCoefficientBase& rhs)
            : Mult(ScalarCoefficient(lhs), rhs)
         {}

         Mult(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
            : Mult(rhs, lhs)
         {}

         Mult(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  MatrixCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  MatrixCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const MatrixCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         int getRows() const override
         {
            return getRHS().getRows();
         }

         int getColumns() const override
         {
            return getRHS().getColumns();
         }

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            getRHS().getValue(value, trans, ip);
            value *= getLHS().getValue(trans, ip);
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };
   Mult(const ScalarCoefficientBase&, const MatrixCoefficientBase&)
      -> Mult<ScalarCoefficientBase, MatrixCoefficientBase>;
   Mult(const MatrixCoefficientBase&, const ScalarCoefficientBase&)
      -> Mult<ScalarCoefficientBase, MatrixCoefficientBase>;
   Mult(double, const MatrixCoefficientBase&)
      -> Mult<ScalarCoefficientBase, MatrixCoefficientBase>;
   Mult(const MatrixCoefficientBase&, double)
      -> Mult<ScalarCoefficientBase, MatrixCoefficientBase>;

   // template <class FEC>
   // Mult<ScalarCoefficientBase, GridFunction<FEC>>
   // operator*(const ScalarCoefficientBase& lhs, const GridFunction<FEC>& rhs)
   // {
   //    return Mult(lhs, ScalarCoefficient(rhs));
   // }
}

#endif
