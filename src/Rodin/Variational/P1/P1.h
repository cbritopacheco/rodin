/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1_H
#define RODIN_VARIATIONAL_P1_P1_H

#include <boost/multi_array.hpp>
#include <functional>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

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
    using ElementType = Variational::P1Element<RangeType>;
  };

  template <class Scalar, class Mesh>
  struct Traits<Variational::P1<Math::Vector<Scalar>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<ScalarType>;
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

  template <class Range, class Mesh = Geometry::Mesh<Context::Local>>
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
  template <class Scalar>
  class P1<Scalar, Geometry::Mesh<Context::Local>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, P1<Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P1 space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

      /**
       * @brief Mapping for the scalar/complex P1 space.
       */
      template <class Callable>
      class Mapping : public FiniteElementSpaceMappingBase<Mapping<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Mapping(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Mapping(const Mapping&) = default;

          decltype(auto) operator()(const Math::SpatialVector<Real>& r) const
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

      /**
       * @brief Inverse mapping for the scalar/complex P1 space.
       */
      template <class Callable>
      class InverseMapping
        : public FiniteElementSpaceInverseMappingBase<InverseMapping<Callable>>
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

      P1(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      virtual ~P1() = default;

      P1& operator=(P1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
        }
        return *this;
      }

      P1& operator=(const P1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
        }
        return *this;
      }

      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Point);
            return s_element;
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Segment);
            return s_element;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Triangle);
            return s_element;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Quadrilateral);
            return s_element;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Tetrahedron);
            return s_element;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Wedge);
            return s_element;
          }
        }
        assert(false);
        static thread_local constexpr ElementType s_null;
        return s_null;
      }

      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount();
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
        return getMesh().getConnectivity().getPolytope(d, i);
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(p.size() >= 0);
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
      template <class Callable>
      auto getMapping(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Mapping(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /**
       * @brief Returns the inverse mapping of the function from the physical
       * element to the reference element.
       * @param[in] idx Index of the element in the mesh.
       * @param[in] v Callable type
       */
      template <class Callable>
      auto getInverseMapping(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return InverseMapping(std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Real, Geometry::Mesh<Context>>;

  /// Alias for a scalar valued P1 finite element space
  template <class Mesh>
  using RealP1 = P1<Real, Mesh>;

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
   * Math::Vector<Scalar> type.
   */
  template <class Scalar>
  class P1<Math::Vector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        P1<Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = Math::Vector<ScalarType>;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Local>;

      /// Represents the Context of the P1 space
      using ContextType = Context::Local;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

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

      P1(const Geometry::Mesh<ContextType>& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        const size_t vn = mesh.getVertexCount();
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = 0; d <= mesh.getDimension(); d++)
        {
          const size_t n = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(n);
          for (size_t i = 0; i < n; i++)
          {
            const auto& polytope = mesh.getConnectivity().getPolytope(d, i);
            const size_t count = polytope.size();
            auto& dofs = m_dofs[d].emplace_back(count * vdim);
            for (size_t local = 0; local < count * vdim; local++)
            {
              const size_t q = local / vdim;
              const size_t r = local % vdim;
              assert(q < count);
              dofs.coeffRef(local) = polytope(q) + r * vn;
            }
          }
        }
      }

      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_dofs(other.m_dofs)
      {}

      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim)),
          m_dofs(std::move(other.m_dofs))
      {}

      virtual ~P1() = default;

      P1& operator=(P1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = std::move(other.m_vdim);
          m_dofs = std::move(other.m_dofs);
        }
        return *this;
      }

      P1& operator=(const P1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
          m_dofs = other.m_dofs;
        }
        return *this;
      }

      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto& g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Point),
              ElementType(1, Geometry::Polytope::Type::Point),
              ElementType(2, Geometry::Polytope::Type::Point),
              ElementType(3, Geometry::Polytope::Type::Point)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Segment),
              ElementType(1, Geometry::Polytope::Type::Segment),
              ElementType(2, Geometry::Polytope::Type::Segment),
              ElementType(3, Geometry::Polytope::Type::Segment)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Triangle),
              ElementType(1, Geometry::Polytope::Type::Triangle),
              ElementType(2, Geometry::Polytope::Type::Triangle),
              ElementType(3, Geometry::Polytope::Type::Triangle)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Quadrilateral),
              ElementType(1, Geometry::Polytope::Type::Quadrilateral),
              ElementType(2, Geometry::Polytope::Type::Quadrilateral),
              ElementType(3, Geometry::Polytope::Type::Quadrilateral)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Tetrahedron),
              ElementType(1, Geometry::Polytope::Type::Tetrahedron),
              ElementType(2, Geometry::Polytope::Type::Tetrahedron),
              ElementType(3, Geometry::Polytope::Type::Tetrahedron)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local constexpr std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(0, Geometry::Polytope::Type::Wedge),
              ElementType(1, Geometry::Polytope::Type::Wedge),
              ElementType(2, Geometry::Polytope::Type::Wedge),
              ElementType(3, Geometry::Polytope::Type::Wedge)
            };
            return s_elements[m_vdim];
          }
        }
        assert(false);
        static thread_local constexpr ElementType s_null(0, Geometry::Polytope::Type::Point);
        return s_null;
      }

      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount() * m_vdim;
      }

      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_dofs[d][i];
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        const size_t q = local / m_vdim;
        const size_t r = local % m_vdim;
        assert(q < static_cast<size_t>(p.size()));
        return p(q) + r * m_mesh.get().getVertexCount();
      }

      template <class Callable>
      auto getMapping(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Mapping(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      template <class Callable>
      auto getInverseMapping(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return InverseMapping(std::forward<Callable>(v));
      }

    private:
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
}

#endif
