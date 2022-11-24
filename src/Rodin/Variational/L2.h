/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_L2_H
#define RODIN_VARIATIONAL_L2_H

#include <functional>

#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"

#include "FiniteElementSpace.h"
#include "FiniteElementCollection.h"

namespace Rodin::Variational
{
   template <class TraitTag>
   class L2 : public FiniteElementSpaceBase
   {
      public:
         using Context = TraitTag;

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

         class FEC : public FiniteElementCollectionBase
         {
            public:
               constexpr
               FEC(const int order, const int elemDim, Basis basis)
                  : m_fec(new mfem::L2_FECollection(order, elemDim, static_cast<int>(basis)))
               {
                  assert(order >= 0);
               }

               constexpr
               FEC(FEC&& other)
                  :  FiniteElementCollectionBase(std::move(other)),
                     m_fec(std::move(other.m_fec)),
                     m_basis(std::move(other.m_basis))
               {}

               constexpr
               FEC& operator=(FEC&& other)
               {
                  FiniteElementCollectionBase::operator=(std::move(other));
                  m_fec = std::move(other.m_fec);
                  m_basis = std::move(other.m_basis);
                  return *this;
               }

               constexpr
               Basis getBasisType() const
               {
                  return m_basis;
               }

               mfem::FiniteElementCollection& getHandle() override
               {
                  assert(m_fec);
                  return *m_fec;
               }

               const mfem::FiniteElementCollection& getHandle() const override
               {
                  assert(m_fec);
                  return *m_fec;
               }

            private:
               std::unique_ptr<mfem::L2_FECollection> m_fec;
               Basis m_basis;
         };

         constexpr
         L2(Geometry::Mesh<Context>& mesh,
               int vdim = 1, int order = 0, Basis basis = DefaultBasis)
            :  m_fec(order, mesh.getDimension(), basis),
               m_mesh(mesh),
               m_fes(new mfem::FiniteElementSpace(
                        &mesh.getHandle(), &m_fec.getHandle(), vdim))
         {}

         constexpr
         L2(const L2& other)
            :  FiniteElementSpaceBase(other),
               m_fec(other.m_fec),
               m_mesh(other.m_mesh),
               m_fes(new mfem::FiniteElementSpace(*other.m_fes))
         {}

         constexpr
         L2(L2&& other)
            :  FiniteElementSpaceBase(std::move(other)),
               m_fec(std::move(other.m_fec)),
               m_mesh(std::move(other.m_mesh)),
               m_fes(std::move(other.m_fes))
         {}

         constexpr
         L2& operator=(L2&& other)
         {
            FiniteElementSpaceBase::operator=(std::move(other));
            m_fec = std::move(other.m_fec);
            m_mesh = std::move(other.m_mesh);
            m_fes = std::move(other.m_fes);
            return *this;
         }

         int getSize() const override
         {
            return getHandle().GetVSize();
         }

         bool isParallel() const override
         {
            return false;
         }

         Geometry::Mesh<Context>& getMesh() override
         {
            return m_mesh;
         }

         const Geometry::Mesh<Context>& getMesh() const override
         {
            return m_mesh;
         }

         const FEC& getFiniteElementCollection() const override
         {
            return m_fec;
         }

         mfem::Array<int> getDOFs(
               const Geometry::SimplexBase& element) const override
         {
            mfem::Array<int> res;
            switch (element.getRegion())
            {
               case Geometry::Region::Domain:
               {
                  m_fes->GetElementVDofs(element.getIndex(), res);
                  break;
               }
               case Geometry::Region::Boundary:
               {
                  m_fes->GetBdrElementVDofs(element.getIndex(), res);
                  break;
               }
               case Geometry::Region::Interface:
               {
                  m_fes->GetFaceVDofs(element.getIndex(), res);
                  break;
               }
            }
            return res;
         }

         const mfem::FiniteElement& getFiniteElement(
               const Geometry::SimplexBase& element) const override
         {
            switch (element.getRegion())
            {
               case Geometry::Region::Domain:
                  return *m_fes->GetFE(element.getIndex());
               case Geometry::Region::Interface:
                  return *m_fes->GetFaceElement(element.getIndex());
               case Geometry::Region::Boundary:
                  return *m_fes->GetBE(element.getIndex());
            }
         }

         mfem::FiniteElementSpace& getHandle() override
         {
            assert(m_fes);
            return *m_fes;
         }

         const mfem::FiniteElementSpace& getHandle() const override
         {
            assert(m_fes);
            return *m_fes;
         }

      private:
         FEC m_fec;
         std::reference_wrapper<Geometry::Mesh<Context>> m_mesh;
         std::unique_ptr<mfem::FiniteElementSpace> m_fes;
   };

   template <class Context>
   L2(Geometry::Mesh<Context>& mesh,
      int vdim = 1,
      int order = 1, typename L2<Context>::Basis basis = L2<Context>::DefaultBasis) -> L2<Context>;
}

#endif

