/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MAGNITUDE_H
#define RODIN_VARIATIONAL_MAGNITUDE_H

#include <utility>
#include <optional>

#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      /**
       * @internal
       */
      class NormL2Coefficient : public mfem::Coefficient
      {
         public:
            NormL2Coefficient(mfem::VectorCoefficient& v)
               : m_v(v)
            {}

            virtual double Eval(
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               mfem::Vector tmp;
               m_v.Eval(tmp, T, ip);
               return tmp.Norml2();
            }

         private:
            mfem::VectorCoefficient& m_v;
      };
   }

   /**
    * @brief Computes the magnitude of the vector.
    *
    * For a vector valued function @f$ u : \Omega \rightarrow \mathbb{R}^d @f$
    * the magnitude at each point @f$ x @f$ is defined by the Euclidean norm:
    * @f[
    *    | u(x) | = \sqrt{ \sum_{i=1}^d u_i(x)^2 }
    * @f]
    */
   class Magnitude : public ScalarCoefficientBase
   {
      public:
         template <class FEC>
         constexpr
         Magnitude(GridFunction<FEC>& u)
            : m_v(new VectorCoefficient(u))
         {}

         Magnitude(const VectorCoefficientBase& v)
            : m_v(v.copy())
         {}

         Magnitude(const Magnitude& other)
            : m_v(other.m_v->copy())
         {}

         void build() override
         {
            m_v->build();
            m_mfemCoefficient.emplace(m_v->get());
         }

         mfem::Coefficient& get() override
         {
            assert(m_mfemCoefficient);
            return *m_mfemCoefficient;
         }

         Magnitude* copy() const noexcept override
         {
            return new Magnitude(*this);
         }
      private:
         std::unique_ptr<VectorCoefficientBase> m_v;
         std::optional<Internal::NormL2Coefficient> m_mfemCoefficient;
   };
}

#include "Dot.hpp"

#endif
