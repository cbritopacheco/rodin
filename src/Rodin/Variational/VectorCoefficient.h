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
          * @brief Builds the underlying mfem::VectorCoefficient object.
          */
         virtual void buildMFEMVectorCoefficient() = 0;

         /**
          * @brief Returns the underlying mfem::VectorCoefficient object.
          * @note Typically one should only call this after one has called
          * buildMFEMVectorCoefficient().
          */
         virtual mfem::VectorCoefficient& getMFEMVectorCoefficient() = 0;

         /**
          * @brief Builds a copy of the object and returns a non-owning
          * pointer to the new object.
          */
         virtual VectorCoefficientBase* copy() const noexcept override = 0;
   };

   template <class ... Values>
   VectorCoefficient(Values&&...) -> VectorCoefficient<Values...>;

   template <class ... Values>
   class VectorCoefficient : public VectorCoefficientBase
   {
      public:
         constexpr
         VectorCoefficient(Values... values)
            :  m_dimension(sizeof...(Values)),
               m_values(std::forward_as_tuple(values...))
         {}

         constexpr
         VectorCoefficient(const VectorCoefficient& other)
            :  m_dimension(other.m_dimension),
               m_values(other.m_values)
         {}

         size_t getDimension() const override;

         void buildMFEMVectorCoefficient() override;

         mfem::VectorCoefficient& getMFEMVectorCoefficient() override;

         VectorCoefficient* copy() const noexcept override
         {
            return new VectorCoefficient(*this);
         }

      private:
         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I == sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>&);

         template<std::size_t I = 0, class ... Tp>
         typename std::enable_if_t<I < sizeof...(Tp)>
         makeCoefficientsFromTuple(const std::tuple<Tp...>& t);

         size_t m_dimension;
         std::tuple<Values...> m_values;
         std::vector<std::unique_ptr<ScalarCoefficientBase>> m_mfemCoefficients;
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
          * @brief Constructs a VectorCoefficient from a vector valued grid
          * function.
          */
         constexpr
         VectorCoefficient(GridFunction<FEC>& u);

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
