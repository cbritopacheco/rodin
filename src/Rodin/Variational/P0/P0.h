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
    using NumberType = Number;
    using RangeType = NumberType;
    using ContextType = typename MeshType::Context;
    using ElementType = Variational::P0Element<RangeType>;
  };

  template <class Number, class Mesh>
  struct Traits<Variational::P0<Math::Vector<Number>, Mesh>>
  {
    using MeshType = Mesh;
    using NumberType = Number;
    using RangeType = Math::Vector<NumberType>;
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

  template <class Range, class Mesh = Geometry::Mesh<Context::Sequential>>
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
  template <class Number>
  class P0<Number, Geometry::Mesh<Context::Sequential>> final
    : public FiniteElementSpace<P0<Number, Geometry::Mesh<Context::Sequential>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using NumberType = Number;

      /// Range type of value
      using RangeType = NumberType;

      /// Represents the Context of the P0 space
      using ContextType = Context::Sequential;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P0Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<P0<RangeType, MeshType>>;

      /**
       * @brief Mapping for the scalar P0 space.
       */
      template <class FunctionDerived>
      class Mapping : public FiniteElementSpaceMappingBase<Mapping<FunctionDerived>>
      {
        public:
          using FunctionType = FunctionBase<FunctionDerived>;

          Mapping(const Geometry::Polytope& polytope, const FunctionType& v)
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

          inline
          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          inline
          constexpr
          const FunctionType& getFunction() const
          {
            assert(m_v);
            return *m_v;
          }

        private:
          Geometry::Polytope m_polytope;
          std::reference_wrapper<const Geometry::PolytopeTransformation> m_trans;
          std::unique_ptr<FunctionType> m_v;
      };

      /**
       * @brief Inverse mapping for the scalar P0 space.
       */
      template <class CallableType>
      class InverseMapping
        : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
      {
        public:
          using FunctionType = CallableType;

          /**
           * @param[in] polytope Reference to polytope on the mesh.
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          InverseMapping(const FunctionType& v)
            : m_v(v)
          {}

          InverseMapping(const InverseMapping&) = default;

          inline
          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          inline
          constexpr
          const FunctionType& getFunction() const
          {
            return m_v.get();
          }

        private:
          std::reference_wrapper<const FunctionType> m_v;
      };

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

      inline
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return s_elements[getMesh().getGeometry(d, i)];
      }

      inline
      size_t getSize() const override
      {
        return m_mesh.get().getCellCount();
      }

      inline
      size_t getVectorDimension() const override
      {
        return 1;
      }

      inline
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      inline
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        assert(d == getMesh().getDimension());
        return m_dofs.at(i);
      }

      inline
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
      template <class FunctionDerived>
      inline
      auto getMapping(const std::pair<size_t, Index>& idx, const FunctionBase<FunctionDerived>& v) const
      {
        const auto [d, i] = idx;
        const auto& mesh = getMesh();
        return Mapping<FunctionDerived>(*mesh.getPolytope(d, i), v);
      }

      template <class FunctionDerived>
      inline
      auto getMapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v) const
      {
        return Mapping<FunctionDerived>(polytope, v);
      }

      /**
       * @brief Returns the inverse mapping of the function from the physical
       * element to the reference element.
       * @param[in] idx Index of the element in the mesh.
       * @param[in] v Callable type
       */
      template <class CallableType>
      inline
      auto getInverseMapping(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return InverseMapping<CallableType>(v);
      }

      template <class CallableType>
      inline
      auto getInverseMapping(const Geometry::Polytope& polytope, const CallableType& v) const
      {
        return InverseMapping<CallableType>(v);
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

  template <class NumberType>
  const Geometry::GeometryIndexed<P0Element<NumberType>>
  P0<NumberType, Geometry::Mesh<Context::Sequential>>::s_elements =
  {
    { Geometry::Polytope::Type::Point, P0Element<NumberType>(Geometry::Polytope::Type::Point) },
    { Geometry::Polytope::Type::Segment, P0Element<NumberType>(Geometry::Polytope::Type::Segment) },
    { Geometry::Polytope::Type::Triangle, P0Element<NumberType>(Geometry::Polytope::Type::Triangle) },
    { Geometry::Polytope::Type::Quadrilateral, P0Element<NumberType>(Geometry::Polytope::Type::Quadrilateral) },
    { Geometry::Polytope::Type::Tetrahedron, P0Element<NumberType>(Geometry::Polytope::Type::Tetrahedron) }
  };
}

#endif

