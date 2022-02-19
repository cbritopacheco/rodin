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
   namespace Internal
   {
      class ScalarCoefficient : public mfem::Coefficient
      {
         public:
            ScalarCoefficient(const ScalarCoefficientBase& s);
            double Eval(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         private:
            std::unique_ptr<ScalarCoefficientBase> m_s;
      };
   }

   /**
    * @brief Abstract base class for objects representing scalar coefficients.
    */
   class ScalarCoefficientBase
      : public FormLanguage::Buildable<Internal::ScalarCoefficient>
   {
      public:
         constexpr
         ScalarCoefficientBase() = default;

         constexpr
         ScalarCoefficientBase(const ScalarCoefficientBase&) = default;

         constexpr
         ScalarCoefficientBase& setTraceDomain(int domain)
         {
            m_traceDomain = domain;
            return *this;
         }

         constexpr
         std::optional<int> getTraceDomain() const
         {
            return m_traceDomain;
         }

         virtual ~ScalarCoefficientBase() = default;

         virtual Restriction<ScalarCoefficientBase> restrictTo(int attr);

         virtual Restriction<ScalarCoefficientBase> restrictTo(
               const std::set<int>& attrs);

         virtual double getValueOnInteriorBoundary(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip);

         virtual double getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) = 0;

         std::unique_ptr<Internal::ScalarCoefficient> build() const override
         {
            return std::make_unique<Internal::ScalarCoefficient>(*this);
         }

         virtual ScalarCoefficientBase* copy() const noexcept override = 0;

      private:
         std::optional<int> m_traceDomain;
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
         ScalarCoefficient(const ScalarCoefficient& other) = default;

         constexpr
         ScalarCoefficient(ScalarCoefficient&&) = default;

         constexpr
         T getValue() const
         {
            return m_x;
         }

         double getValue(mfem::ElementTransformation&, const mfem::IntegrationPoint&) override
         {
            return m_x;
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         const T m_x;
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
         ScalarCoefficient(const GridFunction<FEC>& u)
            :  m_u(u),
               m_mfemCoefficient(&u.getHandle())
         {
            assert(u.getFiniteElementSpace().getRangeDimension() == 1);
         }

         constexpr
         ScalarCoefficient(const ScalarCoefficient& other) = default;

         const GridFunction<FEC>& getValue() const
         {
            return m_u;
         }

         double getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            return m_mfemCoefficient.Eval(trans, ip);
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         const GridFunction<FEC>& m_u;
         mfem::GridFunctionCoefficient m_mfemCoefficient;
   };

   ScalarCoefficient(std::function<double(const double*)>)
      -> ScalarCoefficient<std::function<double(const double*)>>;

   template <>
   class ScalarCoefficient<std::function<double(const double*)>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(std::function<double(const double*)> f)
            : m_f(f),
              m_mfemCoefficient(
                [this](const mfem::Vector& v)
                {
                  return m_f(v.GetData());
                })
         {}

         ScalarCoefficient(const ScalarCoefficient& other) = default;

         std::function<double(const double*)> getValue() const
         {
            return m_f;
         }

         double getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            return m_mfemCoefficient.Eval(trans, ip);
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::function<double(const double*)> m_f;
         mfem::FunctionCoefficient m_mfemCoefficient;
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
            : m_pieces(pieces),
              m_mfemCoefficient(pieces.rbegin()->first) // Maximum attribute
         {
            int maxAttr = m_pieces.rbegin()->first;
            for (int i = 1; i <= maxAttr; i++)
            {
               auto v = m_pieces.find(i);
               if (v != m_pieces.end())
                  m_mfemCoefficient(i) = v->second;
               else
                  m_mfemCoefficient(i) = 0.0;
            }
         }

         ScalarCoefficient(const ScalarCoefficient& other)
            : m_pieces(other.m_pieces)
         {}

         const std::map<int, double>& getValue() const
         {
            return m_pieces;
         }

         double getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            return m_mfemCoefficient.Eval(trans, ip);
         }

         ScalarCoefficient* copy() const noexcept override
         {
            return new ScalarCoefficient(*this);
         }

      private:
         std::map<int, double> m_pieces;
         mfem::PWConstCoefficient m_mfemCoefficient;
   };
}

#include "ScalarCoefficient.hpp"

#endif
