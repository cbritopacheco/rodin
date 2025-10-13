/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Base class for finite element spaces.
   *
   * FiniteElementSpaceBase provides the foundation for defining finite element
   * spaces @f$ V_h \subset V @f$ where @f$ V @f$ is a function space and
   * @f$ V_h @f$ is its finite-dimensional approximation. The finite element
   * space defines the basis functions, degrees of freedom, and mesh association
   * required for finite element computations.
   *
   * ## Mathematical Foundation
   * A finite element space is characterized by:
   * - **Mesh**: Geometric discretization @f$ \mathcal{T}_h @f$
   * - **Element**: Local finite element @f$ (K, P, \Sigma) @f$ where:
   *   - @f$ K @f$ is the reference element geometry
   *   - @f$ P @f$ is the polynomial space
   *   - @f$ \Sigma @f$ is the set of degrees of freedom
   * - **Global Space**: @f$ V_h = \{v \in V : v|_K \in P_K \text{ for all } K \in \mathcal{T}_h\} @f$
   *
   * ## Key Features
   * - **DOF Management**: Global and local degree of freedom indexing
   * - **Basis Functions**: Access to shape functions and their derivatives
   * - **Mesh Association**: Strong coupling with underlying mesh structure
   * - **Conformity**: Support for @f$ H^1 @f$, @f$ H(\text{div}) @f$, @f$ H(\text{curl}) @f$ spaces
   */
  class FiniteElementSpaceBase
  {
    public:
      /**
       * @brief Local degree of freedom indexing structure.
       *
       * This structure provides the mapping between local element degrees of
       * freedom and their global counterparts, essential for assembly operations
       * that need to map local element contributions to the global system.
       */
      struct LocalIndex
      {
        std::pair<size_t, Index> p;  ///< Global index information (dimension, index)
        Index local;                 ///< Local element index
      };

      constexpr
      FiniteElementSpaceBase() = default;

      constexpr
      FiniteElementSpaceBase(const FiniteElementSpaceBase&) = default;

      constexpr
      FiniteElementSpaceBase(FiniteElementSpaceBase&&) = default;

      constexpr
      FiniteElementSpaceBase& operator=(FiniteElementSpaceBase&&) = default;

      constexpr
      FiniteElementSpaceBase& operator=(const FiniteElementSpaceBase&) = default;

      virtual ~FiniteElementSpaceBase() = default;

      constexpr
      bool operator==(const FiniteElementSpaceBase& other) const
      {
        return this == &other;
      }

      constexpr
      bool operator!=(const FiniteElementSpaceBase& other) const
      {
        return this != &other;
      }

      /**
       * @brief Gets the total number of degrees of freedom.
       * @returns Size of the finite element space
       */
      virtual size_t getSize() const = 0;

      /**
       * @brief Gets the dimension of the range space, i.e. the number of
       * components of each basis function.
       * @returns Vector dimension of the finite element space.
       */
      virtual size_t getVectorDimension() const = 0;

      /**
       * @brief Gets the constant reference to the mesh upon which the finite
       * element space is built on.
       * @returns Constant reference to mesh
       */
      virtual const Geometry::MeshBase& getMesh() const = 0;

      /**
       * @brief Gets a set of global degree of freedom indices associated to
       * the polytope of dimension @f$ d @f$ and index @f$ i @f$.
       * @returns Set of indices associated to the @f$ (d, i) @f$-polytope.
       */
      virtual const IndexArray& getDOFs(size_t d, Index i) const = 0;

      /**
       * @brief Gets the global index for the local degree of freedom on the
       * @f$ (d, i) @f$-polytope.
       * @param[in] idx Pair representing the @f$ (d, i) @f$-polytope.
       * @param[in] local Local degree of freedom index.
       */
      virtual Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const
      {
        const auto [d, i] = idx;
        return getDOFs(d, i).coeff(local);
      }
  };

  /**
   * @brief Represernts a finite element space.
   */
  template <class Mesh, class Derived>
  class FiniteElementSpace : public FiniteElementSpaceBase
  {
    public:
      /// Parent class
      using Parent = FiniteElementSpaceBase;

      using MeshType = Mesh;

      constexpr
      FiniteElementSpace() = default;

      constexpr
      FiniteElementSpace(const FiniteElementSpace&) = default;

      constexpr
      FiniteElementSpace(FiniteElementSpace&&) = default;

      constexpr
      FiniteElementSpace& operator=(FiniteElementSpace&&) = default;

      constexpr
      FiniteElementSpace& operator=(const FiniteElementSpace&) = default;

      virtual ~FiniteElementSpace() = default;

      const Mesh& getMesh() const override
      {
        return static_cast<const Derived&>(*this).getMesh();
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      decltype(auto) getFiniteElement(size_t d, Index i) const
      {
        return static_cast<const Derived&>(*this).getFiniteElement(d, i);
      }

      /**
       * @brief Returns the mapping of the function from the physical element
       * to the reference element.
       * @tparam T Callable type
       * @param[in] p Index of the element in the mesh
       * @param[in] v Function defined on an element of the mesh
       *
       * @note CRTP function to be overriden in Derived class.
       */
      template <class Callable>
      decltype(auto) getPullback(const std::pair<size_t, Index>& p, const Callable& v) const
      {
        return static_cast<const Derived&>(*this).getPullback(p, v);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      template <class CallableType>
      decltype(auto) getPushforward(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return static_cast<const Derived&>(*this).getPushforward(idx, v);
      }
  };

  /**
   * @brief Base class for mappings taking functions defined on physical elements
   * to reference elements.
   *
   * For all @f$ \tau \in \mathcal{T}_h @f$ the mapping
   * @f[
   *  \psi : V(\tau) \rightarrow V(K)
   * @f]
   * takes functions defined on the physical Banach space @f$ V(\tau) @f$ to
   * the reference Banach space @f$ V(K) @f$. Here @f$ \tau = x(K) @f$ is the
   * physical element, @f$ K @f$ is the reference element, and @f$ x : K
   * \rightarrow \tau @f$ is the
   * polytope transformation.
   *
   * @see FiniteElementSpaceInverseMappingBase
   */
  template <class Derived>
  class FiniteElementSpacePullbackBase
  {
    public:
      constexpr
      FiniteElementSpacePullbackBase() = default;

      constexpr
      FiniteElementSpacePullbackBase(const FiniteElementSpacePullbackBase&) = default;

      constexpr
      FiniteElementSpacePullbackBase(FiniteElementSpacePullbackBase&&) = default;

      constexpr
      FiniteElementSpacePullbackBase& operator=(FiniteElementSpacePullbackBase&&) = default;

      constexpr
      FiniteElementSpacePullbackBase& operator=(const FiniteElementSpacePullbackBase&) = default;

      virtual ~FiniteElementSpacePullbackBase() = default;

      /**
       * @brief Evaluates the mapped function on the reference coordinates.
       *
       * For the given function @f$ v \in V(\tau) @f$, performs the following
       * evaluation:
       * @f[
       *  \psi(v)(r)
       * @f]
       * on the reference coordinates @f$ r \in K @f$.
       *
       * @note CRTP function to be overriden in Derived class.
       */
      decltype(auto) operator()(const Math::SpatialVector<Real>& r) const
      {
        return static_cast<const Derived&>(*this).operator()(r);
      }
  };

  /**
   * @brief Base class for inverse mappings taking functions defined on
   * reference elements to physical elements.
   *
   * For all @f$ \tau \in \mathcal{T}_h @f$ the inverse mapping
   * @f[
   *  \psi^{-1} : V(K) \rightarrow V(\tau)
   * @f]
   * takes functions defined on the reference Banach space @f$ V(K) @f$ to the
   * physical Banach space @f$ V(K) @f$. Here @f$ \tau = x(K) @f$ is the
   * physical element, @f$ K @f$ is the reference element, and @f$ x : K
   * \rightarrow \tau @f$ is the
   * polytope transformation.
   *
   * @see FiniteElementSpaceMappingBase
   */
  template <class Derived>
  class FiniteElementSpacePushforwardBase
  {
    public:
      constexpr
      FiniteElementSpacePushforwardBase() = default;

      constexpr
      FiniteElementSpacePushforwardBase(const FiniteElementSpacePushforwardBase&) = default;

      constexpr
      FiniteElementSpacePushforwardBase(FiniteElementSpacePushforwardBase&&) = default;

      constexpr
      FiniteElementSpacePushforwardBase& operator=(FiniteElementSpacePushforwardBase&&) = default;

      constexpr
      FiniteElementSpacePushforwardBase& operator=(const FiniteElementSpacePushforwardBase&) = default;

      virtual ~FiniteElementSpacePushforwardBase() = default;

      /**
       * @brief Evaluates the mapped function on the physical coordinates.
       *
       * For the given function @f$ v \in V(K) @f$, performs the following
       * evaluation:
       * @f[
       *  \psi^{-1}(v)(p)
       * @f]
       * on the physical coordinates @f$ p \in \tau @f$.
       *
       * @note CRTP function to be overriden in Derived class.
       */
      decltype(auto) operator()(const Geometry::Point& pc) const
      {
        return static_cast<const Derived&>(*this).operator()(pc);
      }
  };
}

#endif
