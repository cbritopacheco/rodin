/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H
#define RODIN_VARIATIONAL_H1_H

#include <functional>

#include <mfem.hpp>

#include "Rodin/Mesh.h"

#include "ForwardDecls.h"

#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   /**
    * @brief Arbitrary order @f$ H^1(\Omega)^d @f$ conforming (continuous) finite
    * element space supported on the domain @f$ \Omega @f$.
    *
    * Given some discretization @f$ \mathcal{T}_h @f$ (e.g. a triangulation)
    * of @f$ \Omega @f$, instances of this class will represent the space
    * @f[
    *    V_h := \left\{ v \in H^1(\Omega)^d \mid v_{|\tau} \in \mathcal{P},
    *    \quad \forall \tau \in \mathcal{T}_h \right\}
    * @f]
    * where @f$ \mathcal{P} @f$ denotes a vector space of functions from
    * @f$ \tau @f$ to @f$ \mathbb{R}^d @f$.
    *
    */
   class H1 : public FiniteElementSpace<H1>
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

         /**
          * @brief Constructs a finite element space of H1 functions supported
          * on the given mesh.
          *
          * @param[in] mesh Reference to mesh.
          * @param[in] dim Dimension @f$ d @f$ of the range.
          * @param[in] order Order of approximation. Must be greater than or
          * equal to 1.
          * @param[in] basis Type of basis to use.
          */
         H1(Mesh& mesh, const int dim=1, const int order=1, Basis basis = GaussLobato)
            :  m_mesh(mesh),
               m_fec(order, m_mesh.get().getDimension()),
               m_fes(&m_mesh.get().getHandle(), &m_fec, dim),
               m_basis(basis)
         {
            assert(order >= 1);
         }

         /**
          * @brief Gets the type of basis
          * @returns Type of basis
          */
         Basis getBasis() const
         {
            return m_basis;
         }

         /**
          * @brief Returns the reference to the mesh.
          */
         Mesh& getMesh()
         {
            return m_mesh;
         }

         /**
          * @brief Gets dimension of the range space.
          */
         int getDimension() const
         {
            return m_fes.GetVDim();
         }

         /**
          * @internal
          * @brief Returns the underlying mfem::H1_FECollection object.
          * @returns Reference to the underlying mfem::H1_FECollection.
          */
         mfem::H1_FECollection& getFEC()
         {
            return m_fec;
         }

         /**
          * @internal
          * @brief Returns the underlying mfem::FiniteElementSpace object.
          * @returns Reference to the underlying mfem::FiniteElementSpace.
          */
         mfem::FiniteElementSpace& getFES()
         {
            return m_fes;
         }

      private:
         std::reference_wrapper<Mesh> m_mesh;
         mfem::H1_FECollection m_fec;
         mfem::FiniteElementSpace m_fes;
         Basis m_basis;
   };
}

#endif
