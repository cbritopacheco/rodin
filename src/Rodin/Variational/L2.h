/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_L2_H
#define RODIN_VARIATIONAL_L2_H

#include <functional>

#include "Rodin/Mesh.h"

#include "ForwardDecls.h"

#include "FiniteElementCollection.h"

namespace Rodin::Variational
{
   class L2 : public FiniteElementCollectionBase
   {
      public:
         /**
          * @brief Possible types of bases for the H1 finite element space.
          */
         enum class Basis
         {
            /**
             * @brief Gauss Legendre basis (endpoints are not included).
             */
            GaussLegendre        = mfem::BasisType::GaussLegendre,

            /**
             * @brief Gauss Lobatto basis (endpoints are included).
             */
            GaussLobato          = mfem::BasisType::GaussLobatto,

            /**
             * @brief Bernstein polynomial basis.
             */
            Bernstein            = mfem::BasisType::Positive,

            /**
             * @brief Open uniform basis.
             *
             * The nodes @f$ x_i @f$ are defined by:
             * @f[
             *    x_i := \dfrac{i + 1}{n + 1}
             * @f]
             * for @f$ i = 0, \ldots, n - 1 @f$.
             */
            OpenUniform          = mfem::BasisType::OpenUniform,

            /**
             * @brief Closed uniform basis.
             *
             * The nodes @f$ x_i @f$ are defined by:
             * @f[
             *    x_i := \dfrac{i}{n - 1}
             * @f]
             * for @f$ i = 0, \ldots, n - 1 @f$.
             */
            ClosedUniform        = mfem::BasisType::ClosedUniform,

            /**
             * @brief Open-half uniform basis.
             *
             * The nodes @f$ x_i @f$ are defined by:
             * @f[
             *    x_i := \dfrac{i + \frac{1}{2}}{n}
             * @f]
             * for @f$ i = 0, \ldots, n - 1 @f$.
             */
            OpenHalfUniform      = mfem::BasisType::OpenHalfUniform,

            /**
             * @brief Serendipity basis (squares / cubes).
             * @todo Find out exactly which basis this is.
             */
            Serendipity          = mfem::BasisType::Serendipity,

            /**
             * @brief Closed Gauss legendre basis.
             * @todo Find out exactly which basis this is.
             */
            ClosedGaussLegendre  = mfem::BasisType::ClosedGL,

            /**
             * @brief Integrated GLL indicator functions.
             * @todo Find out exactly which basis this is.
             */
            IntegratedGLL        = mfem::BasisType::IntegratedGLL
         };

         static constexpr Basis DefaultBasis = Basis::GaussLegendre;

         L2(const int order, const int elemDim, Basis basis)
            : m_fec(order, elemDim, static_cast<int>(basis))
         {
            assert(order >= 0);
         }

         Basis getBasisType() const
         {
            return m_basis;
         }

         mfem::FiniteElementCollection& getHandle()
         {
            return m_fec;
         }

         const mfem::FiniteElementCollection& getHandle() const
         {
            return m_fec;
         }

      private:
         mfem::L2_FECollection m_fec;
         Basis m_basis;
   };
}

#endif

