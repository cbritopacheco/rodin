/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0_H
#define RODIN_VARIATIONAL_P0_P0_H

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P0Element.h"

namespace Rodin::FormLanguage
{
  template <class Number, class Mesh>
  struct Traits<Variational::P0<Number, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Number;
    using RangeType = ScalarType;
    using ContextType = typename MeshType::Context;
    using ElementType = Variational::P0Element<RangeType>;
  };

  template <class Number, class Mesh>
  struct Traits<Variational::P0<Math::Vector<Number>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Number;
    using RangeType = Math::Vector<ScalarType>;
    using ContextType = typename MeshType::Context;
    using ElementType = Variational::P0Element<RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P0Specializations P0 Template Specializations
   * @brief Template specializations of the P0 class.
   * @see P0
   */

  template <class Range, class Mesh = Geometry::Mesh<Context::Local>>
  class P0;

  /**
   * @ingroup P0Specializations
   * @brief Real valued Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) : v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Real type.
   */
  template <>
  class P0<Real, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, P0<Real, Geometry::Mesh<Context::Local>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Real;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P0 space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P0Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P0<RangeType, MeshType>>;

      template <class Callable>
      class Mapping :
        public FiniteElementSpaceMappingBase<Mapping<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Mapping(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Mapping(const Mapping&) = default;

          decltype(auto) operator()(const Math::SpatialPoint& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      template <class Function>
      Mapping(const Geometry::Polytope&, Function&&)
        -> Mapping<Function>;

      template <class Callable>
      class InverseMapping :
        public FiniteElementSpaceInverseMappingBase<InverseMapping<Callable>>
      {
        public:
          using CallableType = Callable;

          /**
           * @param[in] polytope Reference to polytope on the mesh.
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          template <class Function>
          InverseMapping(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          InverseMapping(const InverseMapping&) = default;

          constexpr
          decltype(auto) operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      template <class Function>
      InverseMapping(Function&&) -> InverseMapping<Function>;

      P0(const MeshType& mesh)
        : m_mesh(mesh)
      {
        const size_t n = mesh.getCellCount();
        m_dofs.reserve(n);
        for (size_t i = 0; i < n; ++i)
          m_dofs.push_back(IndexArray{{i}});
      }

      P0(const P0& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      P0(P0&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh)
      {}

      P0& operator=(P0&& other) = default;

      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return s_elements[getMesh().getGeometry(d, i)];
      }

      size_t getSize() const override
      {
        return m_mesh.get().getCellCount();
      }

      size_t getVectorDimension() const override
      {
        return 1;
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        assert(d == getMesh().getDimension());
        return m_dofs.at(i);
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        assert(d == getMesh().getDimension());
        return i;
      }

      /**
       * @brief Returns the mapping of the function from the physical element
       * to the reference element.
       * @param[in] idx Index of the element in the mesh
       * @param[in] v Function defined on an element of the mesh
       *
       * For all @f$ \tau \in \mathcal{T}_h @f$ in the mesh, the finite
       * element space is generated by the bijective mapping:
       * @f[
       *  \psi_\tau : \mathbb{P}_1(\tau) \rightarrow \mathbb{P}_1(R)
       * @f]
       * taking a function @f$ v \in V(\tau) @f$ from the global element @f$
       * \tau @f$ element @f$ R @f$.
       */
      template <class Function>
      auto getMapping(const std::pair<size_t, Index>& idx, const Function& v) const
      {
        const auto [d, i] = idx;
        const auto& mesh = getMesh();
        return Mapping(*mesh.getPolytope(d, i), v);
      }

      template <class Function>
      auto getMapping(const Geometry::Polytope& polytope, const Function& v) const
      {
        return Mapping(polytope, v);
      }

      /**
       * @brief Returns the inverse mapping of the function from the physical
       * element to the reference element.
       * @param[in] idx Index of the element in the mesh.
       * @param[in] v Callable type
       */
      template <class Function>
      auto getInverseMapping(const std::pair<size_t, Index>& idx, const Function& v) const
      {
        return InverseMapping(v);
      }

      template <class Function>
      auto getInverseMapping(const Geometry::Polytope& polytope, const Function& v) const
      {
        return InverseMapping(v);
      }

    private:
      static const Geometry::GeometryIndexed<ElementType> s_elements;

      std::vector<IndexArray> m_dofs;
      std::reference_wrapper<const MeshType> m_mesh;
  };

  template <class Context>
  P0(const Geometry::Mesh<Context>&) -> P0<Real, Geometry::Mesh<Context>>;

  template <class Mesh>
  using RealP0 = P0<Real, Mesh>;

  template <class Mesh>
  using ComplexP0 = P0<Complex, Mesh>;
}

#endif

