/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H
#define RODIN_VARIATIONAL_H1_H

#include <memory>
#include <functional>

#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"

#include "FiniteElementSpace.h"
#include "FiniteElementCollection.h"

namespace Rodin::Variational
{
  template <class TraitTag>
  class H1 : public FiniteElementSpaceBase
  {
    public:
      using Context = TraitTag;

      /**
       * @brief Possible types of bases for the H1 finite element space.
       */
      enum class Basis
      {
        /**
         * @brief Gauss Lobatto basis (endpoints are included).
         */
        GaussLobato       = mfem::BasisType::GaussLobatto,

        /**
         * @brief Bernstein polynomial basis.
         */
        Bernstein        = mfem::BasisType::Positive,

        /**
         * @brief Closed uniform basis.
         *
         * The nodes @f$ x_i @f$ are defined by:
         * @f[
         *   x_i := \dfrac{i}{n - 1}
         * @f]
         * for @f$ i = 0, \ldots, n - 1 @f$.
         */
        ClosedUniform      = mfem::BasisType::ClosedUniform,

        /**
         * @brief Serendipity basis (squares / cubes).
         * @todo Find out exactly which basis this is.
         */
        Serendipity       = mfem::BasisType::Serendipity,
      };

      static constexpr Basis DefaultBasis = Basis::GaussLobato;

      class FEC : public FiniteElementCollectionBase
      {
        public:
          constexpr
          FEC(const int order, const int elemDim, Basis basis)
            : m_fec(new mfem::H1_FECollection(order, elemDim, static_cast<int>(basis)))
          {
            assert(order >= 1);
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
          std::unique_ptr<mfem::H1_FECollection> m_fec;
          Basis m_basis;
      };

      constexpr
      H1(Geometry::Mesh<Context>& mesh,
          int vdim = 1, int order = 1, Basis basis = DefaultBasis)
        :  m_fec(order, mesh.getDimension(), basis),
          m_mesh(mesh),
          m_fes(new mfem::FiniteElementSpace(
                &mesh.getHandle(), &m_fec.getHandle(), vdim))
      {}

      constexpr
      H1(const H1& other)
        :  FiniteElementSpaceBase(other),
          m_fec(other.m_fec),
          m_mesh(other.m_mesh),
          m_fes(new mfem::FiniteElementSpace(*other.m_fes))
      {}

      constexpr
      H1(H1&& other)
        :  FiniteElementSpaceBase(std::move(other)),
          m_fec(std::move(other.m_fec)),
          m_mesh(std::move(other.m_mesh)),
          m_fes(std::move(other.m_fes))
      {}

      constexpr
      H1& operator=(H1&& other)
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

      mfem::Array<int> getDOFs(const Geometry::Simplex& element) const override
      {
        mfem::Array<int> res;
        if (element.getDimension() == getMesh().getDimension())
        {
          m_fes->GetElementVDofs(element.getIndex(), res);
        }
        else if (element.getDimension() == getMesh().getDimension() - 1)
        {
          m_fes->GetFaceVDofs(element.getIndex(), res);
        }
        else
        {
          assert(false);
        }
        return res;
      }

      const mfem::FiniteElement& getFiniteElement(const Geometry::Simplex& element) const override
      {
        if (element.getDimension() == getMesh().getDimension())
        {
          return *m_fes->GetFE(element.getIndex());
        }
        else if (element.getDimension() == getMesh().getDimension() - 1)
        {
          return *m_fes->GetFaceElement(element.getIndex());
        }
        else
        {
          assert(false);
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

  template <class Trait>
  H1(Geometry::Mesh<Trait>& mesh,
    int vdim = 1,
    int order = 1, typename H1<Trait>::Basis basis = H1<Trait>::DefaultBasis) -> H1<Trait>;
}

#endif
