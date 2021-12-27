/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORCOEFFICIENT_H
#define RODIN_VARIATIONAL_VECTORCOEFFICIENT_H

#include <memory>
#include <optional>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/Base.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @internal
    * @brief Abstract base class for VectorCoefficient objects.
    */
   class VectorCoefficientBase : public FormLanguage::Base
   {
      public:
         /**
          * @brief Gets the dimension of the vector object.
          * @returns Dimension of vector.
          */
         virtual size_t getDimension() const = 0;

         /**
          * @internal
          * @brief Builds the underlying mfem::VectorCoefficient object.
          */
         virtual void buildMFEMVectorCoefficient() = 0;

         /**
          * @internal
          * @brief Returns the underlying mfem::VectorCoefficient object.
          * @note Typically one should only call this after one has called
          * buildMFEMVectorCoefficient().
          */
         virtual mfem::VectorCoefficient& getMFEMVectorCoefficient() = 0;

         /**
          * @internal
          * @brief Builds a copy of the object and returns a non-owning
          * pointer to the new object.
          */
         virtual VectorCoefficientBase* copy() const noexcept override = 0;
   };

   template <class T>
   VectorCoefficient(std::initializer_list<T>)
      -> VectorCoefficient<std::initializer_list<T>>;

   template <class T>
   class VectorCoefficient<std::initializer_list<T>>
      : public VectorCoefficientBase
   {
      static_assert(std::is_convertible_v<T, ScalarCoefficient<T>>);

      public:
         /**
          * @brief Constructs a VectorCoefficient from an initializer list
          */
         constexpr
         VectorCoefficient(std::initializer_list<T> values);

         /**
          * @brief Copies the data.
          *
          * @param other Other coefficient to copy
          */
         constexpr
         VectorCoefficient(const VectorCoefficient& other);

         size_t getDimension() const override;

         void buildMFEMVectorCoefficient() override;

         mfem::VectorCoefficient& getMFEMVectorCoefficient() override;

         VectorCoefficient* copy() const noexcept override
         {
            return new VectorCoefficient(*this);
         }

      private:
         size_t m_dimension;
         std::vector<std::unique_ptr<ScalarCoefficientBase>> m_values;
         std::optional<mfem::VectorArrayCoefficient> m_mfemVectorCoefficient;
   };

   template <class FEC>
   VectorCoefficient(GridFunction<FEC>&)
      -> VectorCoefficient<GridFunction<FEC>>;

   template <class FEC>
   class VectorCoefficient<GridFunction<FEC>>
      : public VectorCoefficientBase
   {
      public:
         /**
          * @brief Constructs a VectorCoefficient from an initializer list
          */
         constexpr
         VectorCoefficient(GridFunction<FEC>& u);

         /**
          * @brief Copies the data.
          *
          * @param other Other coefficient to copy
          */
         constexpr
         VectorCoefficient(const VectorCoefficient& other);

         size_t getDimension() const override;

         void buildMFEMVectorCoefficient() override;

         mfem::VectorCoefficient& getMFEMVectorCoefficient() override;

         VectorCoefficient* copy() const noexcept override
         {
            return new VectorCoefficient(*this);
         }

      private:
         size_t m_dimension;
         GridFunction<FEC>& m_u;
         std::optional<mfem::VectorGridFunctionCoefficient> m_mfemVectorCoefficient;
   };
}

#include "VectorCoefficient.hpp"

#endif
