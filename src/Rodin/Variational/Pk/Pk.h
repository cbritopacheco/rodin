/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PK_PK_H
#define RODIN_VARIATIONAL_PK_PK_H

#include <boost/multi_array.hpp>
#include <functional>
#include <map>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "PkElement.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Scalar, class Mesh>
  struct Traits<Variational::Pk<K, Scalar, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = ScalarType;
    using ElementType = Variational::PkElement<K, RangeType>;
  };

  template <size_t K, class Scalar, class Mesh>
  struct Traits<Variational::Pk<K, Math::Vector<Scalar>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<ScalarType>;
    using ElementType = Variational::PkElement<K, RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup PkSpecializations Pk Template Specializations
   * @brief Template specializations of the Pk class.
   * @see Pk
   */

  template <size_t K, class Range, class Mesh = Geometry::Mesh<Context::Local>>
  class Pk;

  /**
   * @ingroup PkSpecializations
   * @brief Polynomial degree K Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise polynomial functions of degree K:
   * @f[
   *  \mathbb{P}_K (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) : v|_{\tau} \in \mathbb{P}_K(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Real type or Complex type.
   *
   * @tparam K Polynomial degree (0, 1, 2, 3, ...)
   */
  template <size_t K, class Scalar>
  class Pk<K, Scalar, Geometry::Mesh<Context::Local>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, Pk<K, Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the Pk space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = PkElement<K, RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, Pk<K, RangeType, MeshType>>;

      /**
       * @brief Pullback for the scalar/complex Pk space.
       */
      template <class Callable>
      class Pullback : public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Pullback(const Pullback&) = default;

          decltype(auto) operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      /**
       * @brief Inverse Pullback for the scalar/complex Pk space.
       */
      template <class Callable>
      class Pushforward
        : public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          using CallableType = Callable;

          /**
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          template <class Function>
          Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          decltype(auto) operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      /**
       * @brief Constructs a Pk finite element space on the given mesh.
       * @param[in] mesh Mesh on which to build the finite element space
       *
       * Creates the Pk space. The DOF structure depends on the polynomial
       * degree K and the mesh topology.
       */
      Pk(const MeshType& mesh);

      /**
       * @brief Copy constructor.
       * @param[in] other Pk space to copy
       */
      Pk(const Pk& other);

      /**
       * @brief Move constructor.
       * @param[in] other Pk space to move from
       */
      Pk(Pk&& other);

      virtual ~Pk() = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other Pk space to move from
       * @return Reference to this Pk space
       */
      Pk& operator=(Pk&& other);

      /**
       * @brief Copy assignment operator.
       * @param[in] other Pk space to copy
       * @return Reference to this Pk space
       */
      Pk& operator=(const Pk& other);

      /**
       * @brief Gets the finite element associated with a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Reference to the Pk element for this polytope type
       *
       * Returns the appropriate Pk element based on the polytope geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const;

      /**
       * @brief Gets the total number of degrees of freedom.
       * @return Number of DOFs
       *
       * The number of DOFs depends on the polynomial degree K and the
       * mesh structure.
       */
      size_t getSize() const override;

      /**
       * @brief Gets the vector dimension of the space.
       * @return Vector dimension (1 for scalar Pk)
       *
       * Returns the number of components per DOF. For scalar Pk, this is 1.
       */
      size_t getVectorDimension() const override;

      /**
       * @brief Gets the underlying mesh.
       * @return Reference to the mesh
       */
      const MeshType& getMesh() const override;

      /**
       * @brief Gets the global DOF indices for a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Array of global DOF indices
       */
      const IndexArray& getDOFs(size_t d, Index i) const override;

      /**
       * @brief Converts local to global DOF index.
       * @param[in] idx Pair of (dimension, polytope index)
       * @param[in] local Local DOF index within the polytope
       * @return Global DOF index
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override;

      /**
       * @brief Returns the Pullback of the function from the physical element
       * to the reference element.
       * @param[in] idx Index of the element in the mesh
       * @param[in] v Function defined on an element of the mesh
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /**
       * @brief Returns the inverse Pullback of the function from the physical
       * element to the reference element.
       * @param[in] idx Index of the element in the mesh.
       * @param[in] v Callable type
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      void build();

      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<std::vector<IndexArray>> m_dofs;
      size_t m_totalDOFs;
  };

  /// Alias for a scalar valued Pk finite element space
  template <size_t K, class Mesh>
  using RealPk = Pk<K, Real, Mesh>;

  template <size_t K, class Mesh>
  using ComplexPk = Pk<K, Complex, Mesh>;

  /**
   * @ingroup PkSpecializations
   * @brief Vector valued polynomial degree K Lagrange finite element space
   *
   * Represents the finite element space composed of d-dimensional
   * vector valued, continuous, piecewise polynomial functions of degree K:
   * @f[
   *  \mathbb{P}_K (\mathcal{T}_h)^d = \{ v \in C^0(\mathcal{T}_h)^d \mid v|_{\tau} \in \mathbb{P}_K(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is vector valued, i.e. evaluations of the function are of
   * Math::Vector<Scalar> type.
   *
   * @tparam K Polynomial degree (0, 1, 2, 3, ...)
   */
  template <size_t K, class Scalar>
  class Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        Pk<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>>
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

      /// Represents the Context of the Pk space
      using ContextType = Context::Local;

      /// Type of finite element
      using ElementType = PkElement<K, RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, Pk<K, RangeType, MeshType>>;

      template <class Callable>
      class Pullback :
        public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Pullback(const Pullback&) = default;

          decltype(auto) operator()(const Math::SpatialPoint& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      template <class Callable>
      class Pushforward :
        public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          decltype(auto) operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      Pk(const Geometry::Mesh<ContextType>& mesh, size_t vdim);

      Pk(const Pk& other);

      Pk(Pk&& other);

      virtual ~Pk() = default;

      Pk& operator=(Pk&& other);

      Pk& operator=(const Pk& other);

      const ElementType& getFiniteElement(size_t d, Index i) const;

      size_t getSize() const override;

      size_t getVectorDimension() const override;

      const MeshType& getMesh() const override;

      const IndexArray& getDOFs(size_t d, Index i) const override;

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override;

      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      void build();

      std::reference_wrapper<const Geometry::Mesh<ContextType>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<IndexArray>> m_dofs;
      size_t m_totalDOFs;
  };

  /// Alias for a vector valued Pk finite element space
  template <size_t K, class Mesh>
  using VectorPk = Pk<K, Math::Vector<Real>, Mesh>;
}

#include "Pk.hpp"

#endif
