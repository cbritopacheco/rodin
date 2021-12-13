/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_VECTORCOEFFICIENT_H
#define RODIN_VARIATIONAL_VECTORCOEFFICIENT_H

#include <memory>
#include <tuple>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/RodinBase.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @internal
    * @brief Abstract base class for VectorCoefficient objects.
    */
   class VectorCoefficientBase
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
         virtual VectorCoefficientBase* copy() const noexcept = 0;
   };

   /**
    * @brief Represents a vector valued coefficient.
    */
   class VectorCoefficient : public VectorCoefficientBase
   {
      public:
         /**
          * @brief Constructs a VectorCoefficient from a set of values.
          */
         template <class ... Values>
         VectorCoefficient(const Values&... values);

         /**
          * @brief Copies the data.
          *
          * @param other Other coefficient to copy
          */
         VectorCoefficient(const VectorCoefficient& other);

         /**
          * @internal
          * @brief Creates a new object of type VectorCoefficient and returns a
          * non-owning pointer to the new object.
          * @param args Parameters which will be forwarded to the
          * VectorCoefficient constructor.
          * @returns Non-owning pointer to the new object.
          */
         template <class ... Args>
         static VectorCoefficient* create(Args&&... args) noexcept;

         size_t getDimension() const override
         {
            return m_dimension;
         }

         void buildMFEMVectorCoefficient() override;

         virtual VectorCoefficient* copy() const noexcept override;

         mfem::VectorCoefficient& getMFEMVectorCoefficient() override;

      private:
         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I == sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>&);

         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I < sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>& t);

         size_t m_dimension;
         std::vector<std::unique_ptr<ScalarCoefficientBase>> m_values;
         mfem::VectorArrayCoefficient m_mfemVectorArrayCoefficient;
   };
}

#include "VectorCoefficient.hpp"

#endif
