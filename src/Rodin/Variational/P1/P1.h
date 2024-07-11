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
  template <class Scalar, class Mesh>
  struct Traits<Variational::P1<Scalar, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = ScalarType;
    using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;
    using ElementType = Variational::P1Element<RangeType>;
  };

  template <class Scalar, class Mesh>
  struct Traits<Variational::P1<Math::Vector<Scalar>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<ScalarType>;
    using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;
    using ElementType = Variational::P1Element<RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P1Specializations P1 Template Specializations
   * @brief Template specializations of the P1 class.
   * @see P1
   */

  template <class Range, class Mesh = Geometry::Mesh<Context::Sequential>>
  class P1;

  /**
   * @ingroup P1Specializations
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
  class P1<Real, Geometry::Mesh<Context::Sequential>> final
    : public FiniteElementSpace<P1<Real, Geometry::Mesh<Context::Sequential>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Real;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P1 space
      using ContextType = Context::Sequential;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<P1<RangeType, MeshType>>;

      /**
       * @brief Mapping for the scalar/complex P1 space.
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

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(res, p);
          }

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
       * @brief Inverse mapping for the scalar/complex P1 space.
       */
      template <class CallableType>
      class InverseMapping : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
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

          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return getFunction()(res, p.getReferenceCoordinates());
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

      P1(const MeshType& mesh)
        : m_mesh(mesh)
      {}

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
      const ElementType& getFiniteElement(size_t d, Index i) const
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
      const MeshType& getMesh() const override
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
      static const Geometry::GeometryIndexed<P1Element<RangeType>> s_elements;

      std::reference_wrapper<const MeshType> m_mesh;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Real, Geometry::Mesh<Context>>;

  /// Alias for a scalar valued P1 finite element space
  template <class Mesh>
  using RealP1 = P1<Real, Mesh>;

  /**
   * @ingroup P1Specializations
   * @brief Real valued Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) : v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Complex type.
   */
  template <>
  class P1<Complex, Geometry::Mesh<Context::Sequential>> final
    : public FiniteElementSpace<P1<Complex, Geometry::Mesh<Context::Sequential>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Complex;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P1 space
      using ContextType = Context::Sequential;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<P1<RangeType, MeshType>>;

      /**
       * @brief Mapping for the scalar/complex P1 space.
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

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(res, p);
          }

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
       * @brief Inverse mapping for the scalar/complex P1 space.
       */
      template <class CallableType>
      class InverseMapping : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
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

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          constexpr
          const FunctionType& getFunction() const
          {
            return m_v.get();
          }

        private:
          std::reference_wrapper<const FunctionType> m_v;
      };

      P1(const MeshType& mesh)
        : m_mesh(mesh)
      {}

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
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return s_elements[getMesh().getGeometry(d, i)];
      }

      inline
      size_t getSize() const override
      {
        return 2 * m_mesh.get().getVertexCount();
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
        return getMesh().getConnectivity().getPolytope(d, i);
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        const size_t q = local / 2;
        const size_t r = local % 2;
        assert(q < static_cast<size_t>(p.size()));
        return p(q) + r * m_mesh.get().getVertexCount();
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
      static const Geometry::GeometryIndexed<P1Element<RangeType>> s_elements;

      std::reference_wrapper<const MeshType> m_mesh;
  };

  template <class Mesh>
  using ComplexP1 = P1<Complex, Mesh>;

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
   * Math::Vector<Real> type.
   */
  template <>
  class P1<Math::Vector<Real>, Geometry::Mesh<Context::Sequential>> final
    : public FiniteElementSpace<P1<Math::Vector<Real>, Geometry::Mesh<Context::Sequential>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Real;

      /// Range type of value
      using RangeType = Math::Vector<ScalarType>;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Sequential>;

      /// Represents the Context of the P1 space
      using ContextType = Context::Sequential;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<P1<RangeType, MeshType>>;

      template <class FunctionDerived>
      class Mapping
      {
        public:
          using FunctionType = FunctionBase<FunctionDerived>;

          Mapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v)
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(res, p);
          }

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

      template <class CallableType>
      class InverseMapping : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
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

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return getFunction()(res, p.getReferenceCoordinates());
          }

          constexpr
          const FunctionType& getFunction() const
          {
            return m_v.get();
          }

        private:
          std::reference_wrapper<const FunctionType> m_v;
      };

      P1(const Geometry::Mesh<ContextType>& mesh, size_t vdim);

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
      const Geometry::Mesh<ContextType>& getMesh() const override
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

      template <class CallableType>
      inline
      auto getInverseMapping(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return InverseMapping(v);
      }

      template <class CallableType>
      inline
      auto getInverseMapping(const Geometry::Polytope& polytope, const CallableType& v) const
      {
        return InverseMapping(v);
      }

    private:
      static const std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION> s_elements;

      std::reference_wrapper<const Geometry::Mesh<ContextType>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<IndexArray>> m_dofs;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&, size_t)
    -> P1<Math::Vector<Real>, Geometry::Mesh<Context>>;

  /// Alias for a vector valued P1 finite element space
  template <class Mesh>
  using VectorP1 = P1<Math::Vector<Real>, Mesh>;

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Elements();
  }
}

#endif
