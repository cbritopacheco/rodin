/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARCOEFFICIENT_H
#define RODIN_VARIATIONAL_SCALARCOEFFICIENT_H

#include <memory>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/RodinBase.h"

namespace Rodin::Variational
{
   class ScalarCoefficientBase
   {
      public:
         /**
          * @internal
          * @brief Toggles the sign of the coefficient.
          *
          * If @f$ S @f$ is the scalar coefficient then this method applies the
          * following rule:
          *
          * @f[
          *    S \leftarrow - S
          * @f]
          *
          * @returns Reference to self (for method chaining)
          */
         virtual ScalarCoefficientBase& toggleSign() = 0;

         /**
          * @internal
          * @brief Builds the underlying mfem::Coefficient object.
          */
         virtual void buildMFEMCoefficient() = 0;

         /**
          * @internal
          * @brief Returns the underlying mfem::Coefficient object.
          * @note Typically one should only call this after one has called
          * buildMFEMCoefficient().
          */
         virtual mfem::Coefficient& getMFEMCoefficient() = 0;

         /**
          * @internal
          * @brief Builds a copy of the object and returns a non-owning
          * pointer to the new object.
          */
         virtual ScalarCoefficientBase* copy() const noexcept = 0;
   };

   /**
    * @brief A ScalarCoefficient represents the continuous functions that
    * represent the scalar coefficients in a PDE.
    */
   template <class T, class Enable>
   class ScalarCoefficient
   {
      public:
         ScalarCoefficient(T&)
         {}

         ScalarCoefficient(const T&)
         {}
   };

   /**
    * @brief Represents a scalar coefficient of arithmetic type `T`.
    *
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class T>
   class ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const T& x);
         ScalarCoefficient(const ScalarCoefficient& other);

         void buildMFEMCoefficient() override;

         bool isEvaluated() const;
         mfem::Coefficient& getMFEMCoefficient() override;
         ScalarCoefficient& toggleSign() override;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;
         virtual ScalarCoefficient* copy() const noexcept override;

      protected:
         void track(mfem::Coefficient* ptr);

      private:
         T m_x;
         std::unique_ptr<mfem::Coefficient> m_coeff;
   };

   /**
    * Represents the construction from another `Coeff`. The evaluation
    * behaviour is to the copy the nested coefficient and return its
    * evaluation.
    *
    */
   template <class T>
   class ScalarCoefficient<ScalarCoefficient<T>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const ScalarCoefficient<T>& nested);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void buildMFEMCoefficient() override;
         bool isEvaluated() const;
         mfem::Coefficient& getMFEMCoefficient() override;
         ScalarCoefficient& toggleSign() override;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      private:
         std::unique_ptr<ScalarCoefficient<T>> m_nested;
   };

   template <class Lhs, class Rhs>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientSum<Lhs, Rhs>& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void buildMFEMCoefficient() override;
         mfem::Coefficient& getMFEMCoefficient() override;

         ScalarCoefficient& toggleSign() override;

         bool isEvaluated() const;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      protected:
         void track(mfem::Coefficient* ptr);

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>> m_expr;

         std::unique_ptr<ScalarCoefficient<Lhs>> m_lhs;
         std::unique_ptr<ScalarCoefficient<Rhs>> m_rhs;

         std::unique_ptr<mfem::Coefficient> m_coeff;
   };

   template <class T>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientUnaryMinus<T>& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void buildMFEMCoefficient() override;
         mfem::Coefficient& getMFEMCoefficient() override;

         ScalarCoefficient& toggleSign() override;

         bool isEvaluated() const;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientUnaryMinus<T>> m_expr;
   };
}

#include "ScalarCoefficient.hpp"

#endif
