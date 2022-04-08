/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H
#define RODIN_VARIATIONAL_H1_H

#include <functional>

#include "Rodin/Mesh.h"

#include "ForwardDecls.h"

#include "FiniteElementCollection.h"

namespace Rodin::Variational
{
   /**
    * @brief Arbitrary order @f$ H^1(\Omega)^d @f$ conforming (continuous) finite
    * element space.
    *
    * Given some discretization @f$ \mathcal{T}_h @f$ (e.g. a triangulation)
    * of @f$ \Omega @f$, instances of this class will represent the finite
    * element space
    * @f[
    *    V_h := \left\{ v : \overline{\Omega} \rightarrow \mathbb{R}^d \mid
    *       v_{|\tau} \in \mathcal{P}_\tau,
    *    \ \forall \tau \in \mathcal{T}_h \right\}
    * @f]
    * where @f$ \mathcal{P}_\tau \subset H^1(\tau) @f$ and @f$ V_h \subset
    * C^0(\Omega) @f$ so that @f$ V_h \subset H^1(\Omega)^d @f$, i.e. the
    * elements are @f$ H^1 @f$ conforming. The space @f$ P_\tau @f$ depends on
    * the kind of basis chosen.
    *
    */
   class H1 : public FiniteElementCollectionBase
   {
      public:
         /**
          * @brief Possible types of bases for the H1 finite element space.
          */
         enum Basis
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

         static constexpr Basis DefaultBasis = Basis::GaussLobato;

         H1(const int order, const int elemDim, Basis basis)
            : m_fec(order, elemDim, basis)
         {
            assert(order >= 1);
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
         mfem::H1_FECollection m_fec;
         Basis m_basis;
   };
}

#endif
