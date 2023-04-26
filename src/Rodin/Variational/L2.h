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
  template <class Derived, class ContextType>
  class L2Base
    : public FiniteElementSpaceBase
  {
    public:
      using Context = ContextType;
      using Parent = FiniteElementSpaceBase;

      /**
       * @brief Possible types of bases for the H1 finite element space.
       */
      enum class Basis
      {
        /**
         * @brief Gauss Legendre basis (endpoints are not included).
         */
        GaussLegendre      = mfem::BasisType::GaussLegendre,

        /**
         * @brief Gauss Lobatto basis (endpoints are included).
         */
        GaussLobato       = mfem::BasisType::GaussLobatto,

        /**
         * @brief Bernstein polynomial basis.
         */
        Bernstein        = mfem::BasisType::Positive,

        /**
         * @brief Open uniform basis.
         *
         * The nodes @f$ x_i @f$ are defined by:
         * @f[
         *   x_i := \dfrac{i + 1}{n + 1}
         * @f]
         * for @f$ i = 0, \ldots, n - 1 @f$.
         */
        OpenUniform       = mfem::BasisType::OpenUniform,

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
         * @brief Open-half uniform basis.
         *
         * The nodes @f$ x_i @f$ are defined by:
         * @f[
         *   x_i := \dfrac{i + \frac{1}{2}}{n}
         * @f]
         * for @f$ i = 0, \ldots, n - 1 @f$.
         */
        OpenHalfUniform    = mfem::BasisType::OpenHalfUniform,

        /**
         * @brief Serendipity basis (squares / cubes).
         * @todo Find out exactly which basis this is.
         */
        Serendipity       = mfem::BasisType::Serendipity,

        /**
         * @brief Closed Gauss legendre basis.
         * @todo Find out exactly which basis this is.
         */
        ClosedGaussLegendre  = mfem::BasisType::ClosedGL,

        /**
         * @brief Integrated GLL indicator functions.
         * @todo Find out exactly which basis this is.
         */
        IntegratedGLL      = mfem::BasisType::IntegratedGLL
      };

      static constexpr const Basis DefaultBasis = Basis::GaussLegendre;

      class FEC : public FiniteElementCollectionBase
      {
        public:
          constexpr
          FEC(const size_t order, const size_t elemDim, Basis basis)
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

      L2Base(Geometry::Mesh<Context>& mesh,
          const size_t vdim, const size_t order, Basis basis = DefaultBasis)
        : m_fec(order, mesh.getDimension(), basis),
          m_mesh(mesh),
          m_fes(new mfem::FiniteElementSpace(&mesh.getHandle(), &m_fec.getHandle(), vdim))
      {}

      L2Base(const L2Base& other)
        : FiniteElementSpaceBase(other),
          m_fec(other.m_fec),
          m_mesh(other.m_mesh),
          m_fes(new mfem::FiniteElementSpace(*other.m_fes))
      {}

      L2Base(L2Base&& other)
        :  FiniteElementSpaceBase(std::move(other)),
          m_fec(std::move(other.m_fec)),
          m_mesh(std::move(other.m_mesh)),
          m_fes(std::move(other.m_fes))
      {}

      L2Base& operator=(L2Base&& other)
      {
        FiniteElementSpaceBase::operator=(std::move(other));
        m_fec = std::move(other.m_fec);
        m_mesh = std::move(other.m_mesh);
        m_fes = std::move(other.m_fes);
        return *this;
      }

      inline
      size_t getSize() const final override
      {
        return getHandle().GetVSize();
      }

      inline
      Geometry::Mesh<Context>& getMesh() final override
      {
        return m_mesh;
      }

      inline
      const Geometry::Mesh<Context>& getMesh() const final override
      {
        return m_mesh;
      }

      inline
      const FEC& getFiniteElementCollection() const final override
      {
        return m_fec;
      }

      inline
      FiniteElement<Derived>
      getFiniteElement(const Geometry::Polytope& simplex) const
      {
        if (simplex.getDimension() == getMesh().getDimension())
        {
          return FiniteElement<Derived>(simplex, m_fes->GetFE(simplex.getIndex()));
        }
        else if (simplex.getDimension() == getMesh().getDimension() - 1)
        {
          return FiniteElement<Derived>(simplex, m_fes->GetFaceElement(simplex.getIndex()));
        }
        else
        {
          assert(false);
          return FiniteElement<Derived>(simplex, nullptr);
        }
      }

      inline
      mfem::Array<int> getDOFs(const Geometry::Polytope& element) const final override
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

      mfem::FiniteElementSpace& getHandle() const final override
      {
        assert(m_fes);
        return *m_fes;
      }

    private:
      FEC m_fec;
      std::reference_wrapper<Geometry::Mesh<Context>> m_mesh;
      std::unique_ptr<mfem::FiniteElementSpace> m_fes;
  };

  template <class ContextType>
  class L2<ContextType, Scalar>
    : public L2Base<L2<ContextType, Scalar>, ContextType>
  {
    public:
      using Context = ContextType;
      using Parent = L2Base<L2<ContextType, Scalar>, ContextType>;
      using Range = Scalar;
      using Basis = typename Parent::Basis;
      using Parent::operator=;

      static constexpr const Basis DefaultBasis = Basis::GaussLobato;

      L2(Geometry::Mesh<Context>& mesh,
          FiniteElementOrder order = FiniteElementOrder(0),
          Basis basis = DefaultBasis)
        : Parent(mesh, 1, order, basis)
      {}

      L2(const L2& other)
        : Parent(other)
      {}

      L2(L2&& other)
        : Parent(std::move(other))
      {}
  };

  template <class Context>
  L2(Geometry::Mesh<Context>& mesh,
      FiniteElementOrder order = FiniteElementOrder(0),
      typename L2Base<L2<Context, Scalar>, Context>::Basis basis =
        L2Base<L2<Context, Scalar>, Context>::DefaultBasis)
  -> L2<Context, Scalar>;

  template <class ContextType>
  class L2<ContextType, Math::Vector>
    : public L2Base<L2<ContextType, Math::Vector>, ContextType>
  {
    public:
      using Context = ContextType;
      using Parent = L2Base<L2<ContextType, Math::Vector>, ContextType>;
      using Range = Math::Vector;
      using Basis = typename Parent::Basis;
      using Parent::operator=;

      static constexpr const Basis DefaultBasis = Basis::GaussLobato;

      L2(Geometry::Mesh<Context>& mesh,
          size_t vdim,
          FiniteElementOrder order = FiniteElementOrder(0),
          Basis basis = DefaultBasis)
        : Parent(mesh, vdim, order, basis)
      {}

      L2(const L2& other)
        : Parent(other)
      {}

      L2(L2&& other)
        : Parent(std::move(other))
      {}
  };

  template <class Context>
  L2(Geometry::Mesh<Context>& mesh,
      size_t vdim,
      FiniteElementOrder order = FiniteElementOrder(0),
      typename L2Base<L2<Context, Math::Vector>, Context>::Basis basis =
        L2Base<L2<Context, Math::Vector>, Context>::DefaultBasis)
  -> L2<Context, Math::Vector>;
}

#endif

