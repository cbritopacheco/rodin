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
          * @brief Builds the underlying mfem::Coefficient object.
          */
         virtual void buildMFEMCoefficient() = 0;

         /**
          * @brief Returns the underlying mfem::Coefficient object.
          * @note Typically one should only call this after one has called
          * buildMFEMCoefficient().
          */
         virtual mfem::Coefficient& getMFEMCoefficient() = 0;

         /**
          * @brief Builds a copy of the object and returns a non-owning
          * pointer to the new object.
          */
         virtual ScalarCoefficientBase* copy() const noexcept override = 0;
   };

   /**
    * @brief A ScalarCoefficient represents a continuous function that
    * represent some scalar coefficient in a PDE.
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
    * @brief Represents a ScalarCoefficient of arithmetic type `T`.
    *
    * @tparam T Arithmetic type
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class T>
   class ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
      : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs a ScalarCoefficient from an arithmetic value.
          * @param[in] x Arithmetic value
          */
         ScalarCoefficient(const T& x);

         /**
          * @internal
          */
         ScalarCoefficient(const ScalarCoefficient& other);

         void buildMFEMCoefficient() override;
         mfem::Coefficient& getMFEMCoefficient() override;
         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         T m_x;
         std::optional<mfem::ConstantCoefficient> m_mfemCoefficient;
   };

   /**
    * @brief Represents a scalar coefficient which is built from a
    * GridFunction.
    *
    * @tparam FEC Finite element collection
    */
   template <class FEC>
   class ScalarCoefficient<GridFunction<FEC>>
      : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs a ScalarCoefficient from a GridFunction u
          * @param[in] u GridFunction which belongs to the finite element
          * collection FEC
          */
         ScalarCoefficient(GridFunction<FEC>& u);

         /**
          * @internal
          */
         ScalarCoefficient(const ScalarCoefficient& other);

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         GridFunction<FEC>& m_u;
         std::optional<mfem::GridFunctionCoefficient> m_mfemCoefficient;
   };
}

#include "ScalarCoefficient.hpp"

#endif
