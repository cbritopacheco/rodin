/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARCOEFFICIENT_H
#define RODIN_VARIATIONAL_SCALARCOEFFICIENT_H

#include <map>
#include <set>
#include <memory>
#include <optional>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"

#include "FormLanguage/Base.h"
#include "FormLanguage/ForwardDecls.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects representing scalar coefficients.
    */
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

         virtual Restriction<ScalarCoefficientBase> restrictTo(int attr);

         virtual Restriction<ScalarCoefficientBase> restrictTo(
               const std::set<int>& attrs);
   };

   template <class T>
   ScalarCoefficient(const T&)
      -> ScalarCoefficient<std::enable_if_t<std::is_arithmetic_v<T>, T>>;

   /**
    * @brief Represents a ScalarCoefficient of arithmetic type `T`.
    *
    * @tparam T Arithmetic type
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class T>
   class ScalarCoefficient : public ScalarCoefficientBase
   {
      public:
         /**
          * @brief Constructs a ScalarCoefficient from an arithmetic value.
          * @param[in] x Arithmetic value
          */
         constexpr
         ScalarCoefficient(const T& x)
            : m_x(x)
         {}

         constexpr
         ScalarCoefficient(const ScalarCoefficient& other)
            : m_x(other.m_x)
         {}

         constexpr
         ScalarCoefficient(ScalarCoefficient&& other)
            : m_x(std::move(other.m_x))
         {}

         T getValue() const
         {
            return m_x;
         }

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

   template <class FEC>
   ScalarCoefficient(const GridFunction<FEC>&)
      -> ScalarCoefficient<GridFunction<FEC>>;

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
         constexpr
         ScalarCoefficient(GridFunction<FEC>& u);

         constexpr
         ScalarCoefficient(const ScalarCoefficient& other);

         const GridFunction<FEC>& getValue() const
         {
            return m_u;
         }

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

   ScalarCoefficient(std::function<double(const double*)>)
      -> ScalarCoefficient<std::function<double(const double*)>>;

   template <>
   class ScalarCoefficient<std::function<double(const double*)>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(std::function<double(const double*)> f)
            : m_f(f)
         {}

         ScalarCoefficient(const ScalarCoefficient& other)
            : m_f(other.m_f)
         {}

         std::function<double(const double*)> getValue() const
         {
            return m_f;
         }

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override
         {
            assert(m_mfemCoefficient);
            return *m_mfemCoefficient;
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::function<double(const double*)> m_f;
         std::optional<mfem::FunctionCoefficient> m_mfemCoefficient;
   };


   ScalarCoefficient(const std::map<int, double>&)
      -> ScalarCoefficient<std::map<int, double>>;

   ScalarCoefficient(std::initializer_list<std::pair<int, double>>&)
      -> ScalarCoefficient<std::map<int, double>>;

   template <>
   class ScalarCoefficient<std::map<int, double>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const std::map<int, double>& pieces)
            : m_pieces(pieces)
         {}

         ScalarCoefficient(const ScalarCoefficient& other)
            : m_pieces(other.m_pieces)
         {}

         const std::map<int, double>& getValue() const
         {
            return m_pieces;
         }

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override
         {
            assert(m_mfemCoefficient);
            return *m_mfemCoefficient;
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::map<int, double> m_pieces;
         std::optional<mfem::PWConstCoefficient> m_mfemCoefficient;
   };
}

#include "ScalarCoefficient.hpp"

#endif
