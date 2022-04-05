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
   class H1 : public FiniteElementSpaceBase
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
         H1(Mesh& mesh,
               const int dim = 1, const int order = 1, Basis basis = GaussLobato)
            :  m_mesh(mesh),
               m_fec(order, m_mesh.get().getDimension(), basis),
               m_fes(&m_mesh.get().getHandle(), &m_fec, dim),
               m_basis(basis)
         {
            assert(order >= 1);
         }

         /**
          * @brief Gets the type of basis
          * @returns Type of basis
          */
         Basis getBasisType() const
         {
            return m_basis;
         }

         /**
          * @brief Returns the reference to the mesh on which the space is
          * supported.
          */
         Mesh& getMesh() override
         {
            return m_mesh;
         }

         /**
          * @brief Returns a constant reference to the mesh on which the space is
          * supported.
          */
         const Mesh& getMesh() const override
         {
            return m_mesh;
         }

         /**
          * @brief Gets the dimension @f$ d @f$ of the range space.
          */
         int getVectorDimension() const override
         {
            return m_fes.GetVDim();
         }

         int getNumberOfDofs() const override
         {
            return m_fes.GetNDofs();
         }

         void update() override
         {
            m_fes.Update();
         }

         /**
          * @internal
          * @brief Returns the underlying mfem::H1_FECollection object.
          * @returns Reference to the underlying mfem::H1_FECollection.
          */
         mfem::H1_FECollection& getFEC() override
         {
            return m_fec;
         }

         const mfem::H1_FECollection& getFEC() const override
         {
            return m_fec;
         }

         /**
          * @internal
          * @brief Returns the underlying mfem::FiniteElementSpace object.
          * @returns Reference to the underlying mfem::FiniteElementSpace.
          */
         mfem::FiniteElementSpace& getFES() override
         {
            return m_fes;
         }

         const mfem::FiniteElementSpace& getFES() const override
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
