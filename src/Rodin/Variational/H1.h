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
  template <class Derived, class ContextType>
  class H1Base
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

      static constexpr const Basis DefaultBasis = Basis::GaussLobato;

      class FEC : public FiniteElementCollectionBase
      {
        public:
          constexpr
          FEC(const size_t order, const size_t elemDim, Basis basis)
            : m_fec(new mfem::H1_FECollection(order, elemDim, static_cast<int>(basis)))
          {
            assert(order >= 1);
          }

          constexpr
          FEC(FEC&& other)
            : FiniteElementCollectionBase(std::move(other)),
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
      H1Base(const Geometry::Mesh<Context>& mesh,
          const size_t vdim, const size_t order, Basis basis = DefaultBasis)
        : m_fec(order, mesh.getDimension(), basis),
          m_mesh(mesh),
          m_fes(new mfem::FiniteElementSpace(
                &mesh.getHandle(), &m_fec.getHandle(), vdim))
      {
        assert(order >= 1);
      }

      constexpr
      H1Base(const H1Base& other)
        : FiniteElementSpaceBase(other),
          m_fec(other.m_fec),
          m_mesh(other.m_mesh),
          m_fes(new mfem::FiniteElementSpace(*other.m_fes))
      {}

      constexpr
      H1Base(H1Base&& other)
        : FiniteElementSpaceBase(std::move(other)),
          m_fec(std::move(other.m_fec)),
          m_mesh(std::move(other.m_mesh)),
          m_fes(std::move(other.m_fes))
      {}

      constexpr
      H1Base& operator=(H1Base&& other)
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
      const Geometry::Mesh<Context>& getMesh() const final override
      {
        return m_mesh.get();
      }

      inline
      const FEC& getFiniteElementCollection() const final override
      {
        return m_fec;
      }

      inline
      auto getFiniteElement(const Geometry::Simplex& element) const
      {
        return static_cast<const Derived&>(*this).getFiniteElement(element);
      }

      inline
      mfem::Array<int> getDOFs(const Geometry::Simplex& element) const final override
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
      std::reference_wrapper<const Geometry::Mesh<Context>> m_mesh;
      std::unique_ptr<mfem::FiniteElementSpace> m_fes;
  };

  template <class ContextType>
  class H1<Scalar, ContextType> final
    : public H1Base<H1<Scalar, ContextType>, ContextType>
  {
    public:
      using Context = ContextType;
      using Parent = H1Base<H1<Scalar, ContextType>, ContextType>;
      using RangeType = Scalar;
      using Basis = typename Parent::Basis;
      using Parent::operator=;
      static constexpr const Basis DefaultBasis = Basis::GaussLobato;

      constexpr
      H1(const Geometry::Mesh<Context>& mesh,
          FiniteElementOrder order = FiniteElementOrder(1),
          Basis basis = DefaultBasis)
        : Parent(mesh, 1, order, basis)
      {}

      constexpr
      H1(const H1& other)
        : Parent(other)
      {}

      constexpr
      H1(H1&& other)
        : Parent(std::move(other))
      {}

      constexpr
      H1& operator=(H1&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      inline
      FiniteElement<H1<Scalar, Context>>
      getFiniteElement(const Geometry::Simplex& element) const
      {
        if (element.getDimension() == this->getMesh().getDimension())
        {
          return FiniteElement<H1<Scalar, Context>>(
              element, this->getHandle().GetFE(element.getIndex()));
        }
        else if (element.getDimension() == this->getMesh().getDimension() - 1)
        {
          return FiniteElement<H1<Scalar, Context>>(
              element, this->getHandle().GetFaceElement(element.getIndex()));
        }
        else
        {
          assert(false);
          return FiniteElement<H1<Scalar, Context>>(element, nullptr);
        }
      }
  };

  template <class Context>
  H1(Geometry::Mesh<Context>& mesh,
      FiniteElementOrder order = FiniteElementOrder(1),
      typename H1Base<H1<Scalar, Context>, Context>::Basis basis =
      H1Base<H1<Scalar, Context>, Context>::DefaultBasis) -> H1<Scalar, Context>;

  template <class ContextType>
  class H1<Math::Vector, ContextType> final
    : public H1Base<H1<Math::Vector, ContextType>, ContextType>
  {
    public:
      using Context = ContextType;
      using Parent = H1Base<H1<Math::Vector, ContextType>, ContextType>;
      using RangeType = Math::Vector;
      using Basis = typename Parent::Basis;
      using Parent::operator=;
      static constexpr const Basis DefaultBasis = Basis::GaussLobato;

      constexpr
      H1(const Geometry::Mesh<Context>& mesh,
          size_t vdim,
          FiniteElementOrder order = FiniteElementOrder(1),
          Basis basis = DefaultBasis)
        : Parent(mesh, vdim, order, basis)
      {
        m_fe.resize(mesh.getDimension() + 1);
        for (size_t i = 0; i < mesh.getDimension() + 1; i++)
          m_fe[i].resize(mesh.getSimplexCount(i));
      }

      constexpr
      H1(const H1& other)
        : Parent(other)
      {}

      constexpr
      H1(H1&& other)
        : Parent(std::move(other))
      {}

      constexpr
      H1& operator=(H1&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      inline
      const FiniteElement<H1<Math::Vector, Context>>&
      getFiniteElement(const Geometry::Simplex& simplex) const
      {
        assert(m_fe.size() > simplex.getDimension());
        assert(m_fe[simplex.getDimension()].size() > simplex.getIndex());
        auto& fe = m_fe[simplex.getDimension()][simplex.getIndex()];
        if (fe.has_value())
        {
          return fe.value();
        }
        else
        {
          if (simplex.getDimension() == this->getMesh().getDimension())
          {
            fe.emplace(this->getVectorDimension(), simplex, this->getHandle().GetFE(simplex.getIndex()));
            return fe.value();
          }
          else if (simplex.getDimension() == this->getMesh().getDimension() - 1)
          {
            fe.emplace(this->getVectorDimension(), simplex, this->getHandle().GetFaceElement(simplex.getIndex()));
            return fe.value();
          }
          else
          {
            assert(false);
            return FiniteElement<H1<Math::Vector, Context>>(0, simplex, nullptr);
          }
        }
      }

    private:
      mutable std::vector<std::vector<std::optional<FiniteElement<H1<Math::Vector, Context>>>>> m_fe;
  };

  template <class Context>
  H1(const Geometry::Mesh<Context>& mesh,
      size_t vdim,
      FiniteElementOrder order = FiniteElementOrder(1),
      typename H1Base<H1<Math::Vector, Context>, Context>::Basis basis =
        H1Base<H1<Math::Vector, Context>, Context>::DefaultBasis) -> H1<Math::Vector, Context>;
}

#endif
