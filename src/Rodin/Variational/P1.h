/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_H
#define RODIN_VARIATIONAL_P1_H

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "ForwardDecls.h"
#include "P1Element.h"
#include "GridFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  template <>
  class P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpaceBase
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using RangeType = Scalar;
      using Context = Context::Serial;
      using Element = P1Element<RangeType>;
      using Parent = FiniteElementSpaceBase;

      P1(const Geometry::Mesh<Context>& mesh);

      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_elements(other.m_elements)
      {}

      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh),
          m_elements(other.m_elements)
      {}

      P1& operator=(P1&& other) = default;

      inline
      const ScalarP1Element& getFiniteElement(size_t d, Index i) const
      {
        return m_elements[d][i];
      }

      inline
      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount();
      }

      inline
      size_t getVectorDimension() const override
      {
        return 1;
      }

      inline
      const Geometry::Mesh<Context>& getMesh() const override
      {
        return m_mesh.get();
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(local < static_cast<size_t>(p.size()));
        return p(local);
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<Context>> m_mesh;
      std::vector<std::vector<P1Element<Scalar>>> m_elements;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Scalar, Context, Geometry::Mesh<Context>>;

  template <>
  class P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpaceBase
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using RangeType = Math::Vector;
      using Context = Context::Serial;
      using Element = P1Element<RangeType>;
      using Parent = FiniteElementSpaceBase;

      P1(const Geometry::Mesh<Context>& mesh, size_t vdim);

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      inline
      const VectorP1Element& getFiniteElement(size_t d, Index i) const
      {
        return m_elements[d][i];
      }

      inline
      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount() * m_vdim;
      }

      inline
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      inline
      const Geometry::Mesh<Context>& getMesh() const override
      {
        return m_mesh.get();
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        const size_t q = local / m_vdim;
        const size_t r = local % m_vdim;
        assert(q < static_cast<size_t>(p.size()));
        return p(q) + r * m_mesh.get().getVertexCount();
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<Context>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<P1Element<Math::Vector>>> m_elements;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&, size_t) -> P1<Math::Vector, Context, Geometry::Mesh<Context>>;

  template <class ... Ts>
  class GridFunction<P1<Ts...>> final : public GridFunctionBase<GridFunction<P1<Ts...>>, P1<Ts...>>
  {
    public:
      using FES = P1<Ts...>;
      using Parent = GridFunctionBase<GridFunction<P1<Ts...>>, P1<Ts...>>;

      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;

      /**
       * @brief Constructs a grid function on a finite element space.
       * @param[in] fes Finite element space to which the function belongs
       * to.
       */
      GridFunction(const FES& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copies the grid function.
       * @param[in] other Other grid function to copy.
       */
      GridFunction(const GridFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructs the grid function.
       * @param[in] other Other grid function to move.
       */
      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Move assignment operator.
       */
      inline
      constexpr
      GridFunction& operator=(GridFunction&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      GridFunction& operator=(const GridFunction&)  = delete;
  };

  template <class ... Ts>
  GridFunction(const P1<Ts...>&) -> GridFunction<P1<Ts...>>;
}



#endif