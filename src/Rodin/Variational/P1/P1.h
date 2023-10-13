/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1_H
#define RODIN_VARIATIONAL_P1_P1_H

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P1Element.h"

namespace Rodin::FormLanguage
{
  template <class ContextType, class MeshType>
  struct Traits<Variational::P1<Scalar, ContextType, MeshType>>
  {
    using RangeType = Scalar;
    using Context = ContextType;
    using Mesh = MeshType;
    using Element = Variational::ScalarP1Element;
  };

  template <class ContextType, class MeshType>
  struct Traits<Variational::P1<Math::Vector, ContextType, MeshType>>
  {
    using RangeType = Math::Vector;
    using Context = ContextType;
    using Mesh = MeshType;
    using Element = Variational::VectorP1Element;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P1Specializations P1 Template Specializations
   * @brief Template specializations of the P1 class.
   * @see P1
   */

  template <class Range, class Context, class Mesh = Geometry::Mesh<Context>>
  class P1;

  /**
   * @ingroup P1Specializations
   * @brief Scalar valued Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) \mid v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Scalar type.
   */
  template <>
  class P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpace<P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Serial>;

      /// Range type of value
      using RangeType = Scalar;

      /// Represents the Context of the P1 space
      using Context = Context::Serial;

      //// Type of finite element
      using Element = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<P1<Scalar, Context, Geometry::Mesh<Context>>>;

      /**
       * @brief Mapping for the scalar P1 space.
       */
      template <class FunctionDerived>
      class Mapping : public FiniteElementSpaceMappingBase<Mapping<FunctionDerived>>
      {
        public:
          using Function = FunctionBase<FunctionDerived>;

          Mapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v)
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

          inline
          auto operator()(const Math::SpatialVector& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          inline
          constexpr
          const Function& getFunction() const
          {
            assert(m_v);
            return *m_v;
          }

        private:
          Geometry::Polytope m_polytope;
          std::reference_wrapper<const Geometry::PolytopeTransformation> m_trans;
          std::unique_ptr<Function> m_v;
      };

      /**
       * @brief Inverse mapping for the scalar P1 space.
       */
      template <class CallableType>
      class InverseMapping : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
      {
        public:
          using Function = CallableType;

          /**
           * @param[in] polytope Reference to polytope on the mesh.
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          InverseMapping(const Geometry::Polytope& polytope, const CallableType& v)
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          InverseMapping(const InverseMapping&) = default;

          inline
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          inline
          constexpr
          const Function& getFunction() const
          {
            return m_v.get();
          }

        private:
          Geometry::Polytope m_polytope;
          std::reference_wrapper<const Geometry::PolytopeTransformation> m_trans;
          std::reference_wrapper<Function> m_v;
      };

      P1(const Geometry::Mesh<Context>& mesh);

      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh)
      {}

      P1& operator=(P1&& other) = default;

      inline
      const ScalarP1Element& getFiniteElement(size_t d, Index i) const
      {
        return s_elements[getMesh().getGeometry(d, i)];
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
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return getMesh().getConnectivity().getPolytope(d, i);
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(local < static_cast<size_t>(p.size()));
        return p(local);
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
        return Mapping(*mesh.getPolytope(d, i), v);
      }

      template <class FunctionDerived>
      inline
      auto getMapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v) const
      {
        return Mapping(polytope, v);
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
        const auto [d, i] = idx;
        const auto& mesh = getMesh();
        return InverseMapping(*mesh.getPolytope(d, i), v);
      }

      template <class CallableType>
      inline
      auto getInverseMapping(const Geometry::Polytope& polytope, const CallableType& v) const
      {
        return InverseMapping(polytope, v);
      }

    private:
      static const Geometry::GeometryIndexed<ScalarP1Element> s_elements;

      std::reference_wrapper<const Geometry::Mesh<Context>> m_mesh;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Scalar, Context, Geometry::Mesh<Context>>;

  /// Alias for a scalar valued P1 finite element space
  template <class Context>
  using ScalarP1 = P1<Scalar, Context, Geometry::Mesh<Context>>;

  /**
   * @ingroup P1Specializations
   * @brief Vector valued Lagrange finite element space
   *
   * Represents the finite element space composed of @f$ d @f$ dimensional
   * vector valued, continuous, piecewise linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h)^d = \{ v \in C^0(\mathcal{T}_h)^d \mid v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is vector valued, i.e. evaluations of the function are of
   * Math::Vector type.
   */
  template <>
  class P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpaceBase
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      /// Range value type
      using RangeType = Math::Vector;

      /// Context type
      using Context = Context::Serial;

      /// Type of finite element
      using Element = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpaceBase;

      template <class FunctionDerived>
      class Mapping
      {
        public:
          using Function = FunctionBase<FunctionDerived>;

          Mapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v)
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

          inline
          auto operator()(const Math::Vector& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          inline
          constexpr
          const Function& getFunction() const
          {
            assert(m_v);
            return *m_v;
          }

        private:
          Geometry::Polytope m_polytope;
          std::reference_wrapper<const Geometry::PolytopeTransformation> m_trans;
          std::unique_ptr<Function> m_v;
      };

      P1(const Geometry::Mesh<Context>& mesh, size_t vdim);

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      inline
      const VectorP1Element& getFiniteElement(size_t d, Index i) const
      {
        return s_elements[m_vdim][getMesh().getGeometry(d, i)];
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
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_dofs[d][i];
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

      template <class FunctionDerived>
      inline
      auto getMapping(const std::pair<size_t, Index>& idx, const FunctionBase<FunctionDerived>& v) const
      {
        const auto [d, i] = idx;
        const auto& mesh = getMesh();
        return Mapping(*mesh.getPolytope(d, i), v);
      }

    private:
      static const std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION> s_elements;

      std::reference_wrapper<const Geometry::Mesh<Context>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<IndexArray>> m_dofs;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&, size_t) -> P1<Math::Vector, Context, Geometry::Mesh<Context>>;

  /// Alias for a vector valued P1 finite element space
  template <class Context>
  using VectorP1 = P1<Math::Vector, Context, Geometry::Mesh<Context>>;

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Elements();
  }
}

#endif
