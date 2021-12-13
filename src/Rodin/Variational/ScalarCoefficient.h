/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARCOEFFICIENT_H
#define RODIN_VARIATIONAL_SCALARCOEFFICIENT_H

#include <memory>
#include <optional>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"

#include "FormLanguage/Base.h"
#include "FormLanguage/ForwardDecls.h"

namespace Rodin::Variational
{
   class ScalarCoefficientBase : public FormLanguage::Base
   {
      public:
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
         virtual ScalarCoefficientBase* copy() const noexcept override = 0;
   };

   /**
    * @internal
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

         mfem::Coefficient& getMFEMCoefficient() override;

         virtual ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         T m_x;
         std::optional<mfem::ConstantCoefficient> m_mfemCoefficient;
   };

   template <class FEC>
   class ScalarCoefficient<GridFunction<FEC>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(GridFunction<FEC>& u);

         ScalarCoefficient(const ScalarCoefficient& other);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         virtual ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         GridFunction<FEC>& m_u;
         std::optional<mfem::GridFunctionCoefficient> m_mfemCoefficient;
   };

   /**
    * @internal
    * @brief Represents the sum of two ScalarCoefficient objects, which itself
    * is a ScalarCoefficient
    */
   template <>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientSum>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientSum& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void buildMFEMCoefficient() override;
         mfem::Coefficient& getMFEMCoefficient() override;

         virtual ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientSum> m_expr;
         std::optional<mfem::SumCoefficient> m_mfemCoefficient;
   };

   /**
    * @internal
    * @brief Represents the negation of the scalar coefficient.
    */
   template <>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientUnaryMinus& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void buildMFEMCoefficient() override;
         mfem::Coefficient& getMFEMCoefficient() override;

         virtual ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientUnaryMinus> m_expr;
         std::optional<mfem::SumCoefficient> m_mfemCoefficient;
   };
}

#include "ScalarCoefficient.hpp"

#endif
