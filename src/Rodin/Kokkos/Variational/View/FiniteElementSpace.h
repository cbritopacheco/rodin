/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_VIEW_FINITEELEMENTSPACE_H
#define RODIN_KOKKOS_VARIATIONAL_VIEW_FINITEELEMENTSPACE_H

/**
 * @file
 * @brief Physical basis evaluation and numbering of finite-element spaces.
 */

#include "Common.h"

namespace Rodin::Kokkos
{
  /**
   * @brief Physical element evaluation and coefficient numbering of a Kokkos space.
   * @tparam FES Real scalar or vector finite-element space exposing getDOFsView().
   *
   * ## Architecture
   *
   * The view retains the space's prepared device DOF handles and its physical value
   * dimension. Scalar reference bases are replicated by component, preserving
   * the P1, H1 and P0g numbering and sparse component optimizations. An intrinsic
   * vector basis has one coefficient per basis function, independently of its
   * physical value dimension. Its DOF view supplies physical basis values and
   * the supported differential operators, including any Piola mapping and
   * cell orientation transformations. No nodal interpolation is assumed here.
   *
   * getCount() and getGlobalIndex() describe complete local basis functions.
   * getComponentCount(), getStride() and getCapacity() retain the grouped
   * representation used by existing assembly kernels; intrinsic vector bases
   * have one component group and zero component stride. All data are prepared
   * by the space before constructing this view; evaluation is const.
   */
  template <Variational::Space FES>
    requires(Variational::ScalarSpace<FES> || Variational::VectorSpace<FES>)
  struct View<FES>
  {
    /// @brief Device DOF and family evaluation tables.
    using DOFsType = typename FES::DOFsType;

    /// @brief Physical range of the element.
    using RangeType = typename FES::RangeType;

    /// @brief A scalar generator or an intrinsically vector-valued basis.
    using BasisValue = std::decay_t<decltype(std::declval<const DOFsType&>().getBasis(
      std::declval<const Variational::IntegrationPointView&>(), Index(0)))>;

    /// @brief Whether independent scalar generators are replicated by component.
    static constexpr bool Componentwise = std::is_same_v<BasisValue, Real>;

    /// @brief Maximum number of independent copies of a basis generator.
    static constexpr Index ComponentCapacity =
      Componentwise && !std::is_same_v<RangeType, Real> ? RODIN_MAXIMAL_SPACE_DIMENSION
                                                        : 1;

    /// @param fes Space whose prepared device tables are borrowed.
    explicit View(const FES& fes)
      : m_dofs(fes.getDOFsView()),
        m_dimension(fes.getVectorDimension())
    {}

    /// @returns Prepared device tables, without changing their ownership.
    KOKKOS_INLINE_FUNCTION const DOFsType& getDOFs() const
    {
      return m_dofs;
    }

    /// @returns Physical value dimension, independent of coefficient multiplicity.
    KOKKOS_INLINE_FUNCTION Index getVectorDimension() const
    {
      if constexpr (std::is_same_v<RangeType, Real>)
        return 1;
      else
        return m_dimension;
    }

    /// @returns Number of independently numbered component copies.
    KOKKOS_INLINE_FUNCTION Index getComponentCount() const
    {
      if constexpr (!Componentwise || std::is_same_v<RangeType, Real>)
        return 1;
      else
        return m_dimension;
    }

    /// @returns Offset between independent component copies, zero for intrinsic bases.
    KOKKOS_INLINE_FUNCTION Index getStride() const
    {
      if constexpr (Componentwise)
        return m_dofs.getStride();
      else
        return 0;
    }

    /// @param vertices Vertex capacity. @returns Maximum number of basis generators.
    KOKKOS_INLINE_FUNCTION static constexpr Index getCapacity(Index vertices)
    {
      return DOFsType::getCapacity(vertices);
    }

    /// @param polytope Bound entity. @returns Its generator indices in local order.
    KOKKOS_INLINE_FUNCTION auto bind(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope) const
    {
      return m_dofs.bind(polytope);
    }

    /// @param polytope Bound entity. @returns Complete local basis count.
    KOKKOS_INLINE_FUNCTION Index getCount(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope) const
    {
      return bind(polytope).size() * getComponentCount();
    }

    /**
     * @param polytope Bound entity.
     * @param local Complete local basis index.
     * @returns Coefficient index of the basis function, in mesh-local numbering.
     */
    KOKKOS_INLINE_FUNCTION Index getGlobalIndex(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Polytope>& polytope,
      Index local) const
    {
      const Index components = getComponentCount();
      const auto indices = bind(polytope);
      KOKKOS_ASSERT(local < indices.size() * components);
      return indices[local / components] + (local % components) * getStride();
    }

    /// @param ip Integration point. @param dof Generator index. @returns Generator value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getGenerator(const Point& ip, Index dof) const
    {
      if constexpr (Componentwise)
        return m_dofs.getBasis(ip, dof);
      else
      {
        const Geometry::Vector result = m_dofs.getBasis(ip, dof);
        return result;
      }
    }

    /// @param ip Integration point. @param local Complete local basis index. @returns Physical basis value.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getBasis(const Point& ip, Index local) const
    {
      const Index components = getComponentCount();
      return getBasis(ip, local / components, local % components);
    }

    /**
     * @param ip Integration point. @param dof Generator index.
     * @param component Independent copy index; zero for intrinsic bases.
     * @returns Complete physical basis value.
     */
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getBasis(
      const Point& ip, Index dof, Index component) const
    {
      KOKKOS_ASSERT(component < getComponentCount());
      if constexpr (Componentwise && !std::is_same_v<RangeType, Real>)
      {
        return Variational::UnitVector{
          m_dimension, component, m_dofs.getBasis(ip, dof)};
      }
      else
        return getGenerator(ip, dof);
    }

    /// @param ip Integration point. @param dof Scalar generator index. @returns Physical gradient.
    template <class Point>
    KOKKOS_INLINE_FUNCTION decltype(auto) getGradient(const Point& ip, Index dof) const
      requires requires(const DOFsType& dofs) { dofs.getGradient(ip, dof); }
    {
      return m_dofs.getGradient(ip, dof);
    }

    /**
     * @param ip Integration point. @param dof Generator index.
     * @param component Independent copy index; zero for intrinsic bases.
     * @returns Complete physical basis Jacobian, when supplied by the family.
     */
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getJacobian(
      const Point& ip, Index dof, Index component) const
      requires((Componentwise && requires(const DOFsType& dofs) {
        dofs.getGradient(ip, dof);
      }) || (!Componentwise && requires(const DOFsType& dofs) {
        dofs.getJacobian(ip, dof);
      }))
    {
      KOKKOS_ASSERT(component < getComponentCount());
      if constexpr (Componentwise)
      {
        return Variational::RowMatrix{
          m_dimension, component, m_dofs.getGradient(ip, dof)};
      }
      else
      {
        const Geometry::Matrix result = m_dofs.getJacobian(ip, dof);
        return result;
      }
    }

    /**
     * @param ip Integration point. @param dof Generator index.
     * @param component Independent copy index; zero for intrinsic bases.
     * @returns Physical divergence, evaluated directly by an intrinsic family.
     */
    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getDivergence(
      const Point& ip, Index dof, Index component) const
      requires Variational::VectorSpace<FES> &&
      ((Componentwise && requires(const DOFsType& dofs) {
        dofs.getGradient(ip, dof);
      }) || (!Componentwise && requires(const DOFsType& dofs) {
        dofs.getDivergence(ip, dof);
      }))
    {
      KOKKOS_ASSERT(component < getComponentCount());
      if constexpr (Componentwise)
        return m_dofs.getGradient(ip, dof)(component);
      else
        return m_dofs.getDivergence(ip, dof);
    }

    /// @param ip Integration point. @param dof Generator index. @returns Family-provided physical curl.
    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getCurl(const Point& ip, Index dof) const
      requires requires(const DOFsType& dofs) { dofs.getCurl(ip, dof); }
    {
      return m_dofs.getCurl(ip, dof);
    }

    private:
      DOFsType m_dofs;
      Index m_dimension;
  };
}

#endif
