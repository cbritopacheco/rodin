/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_PROJECTION_H
#define RODIN_KOKKOS_VARIATIONAL_PROJECTION_H

/**
 * @file
 * @brief Device projection of functions on grid functions of device spaces.
 */

#include <cstdint>
#include <string>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Variational/GridFunction.h"

#include "DeviceView.h"
#include "Integral.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Constant value of a projection.
   * @tparam Value Real or device vector.
   */
  template <class Value>
  struct ConstantView
  {
    /// @param p Device point. @return The constant.
    KOKKOS_INLINE_FUNCTION
    Value getValue([[maybe_unused]] const Geometry::PointView& p) const
    {
      return value;
    }

    /// @brief Constant value.
    Value value;
  };

  /**
   * @brief Callable taking a device point, as the source of a projection.
   * @tparam F Callable returning a real or a device vector.
   */
  template <class F>
  struct CallableView
  {
    /// @param p Device point. @return Value of the callable.
    KOKKOS_INLINE_FUNCTION
    auto getValue(const Geometry::PointView& p) const
    {
      return f(p);
    }

    /// @brief Device copy of the callable.
    F f;
  };

  /// @brief Device predicate selecting every polytope of a projection.
  struct SelectAll
  {
    /// @param p Device polytope. @return true.
    KOKKOS_INLINE_FUNCTION
    bool operator()([[maybe_unused]] const Geometry::PolytopeView& p) const
    {
      return true;
    }
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Projection on grid functions of device spaces of a resident
   * Kokkos mesh.
   *
   * The source is evaluated on the device at the node of every degree of
   * freedom of the selected polytopes, as the host evaluates the linear
   * forms of the elements. When polytopes share a degree of freedom, the one
   * of largest index sets it, as the last one visited on the host. Degrees of
   * freedom outside the selection keep their values.
   *
   * The source is a function expression with a device lowering, a callable
   * taking a const Rodin::Kokkos::Geometry::PointView&, or a constant. A
   * device predicate, callable on a const Rodin::Kokkos::Geometry::PolytopeView&
   * and returning bool, further restricts the selected polytopes.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  struct DeviceProjection<Context::Kokkos<Context::Local, Backend>>
  {
    /// @brief Mesh context.
    using ContextType = Context::Kokkos<Context::Local, Backend>;

    /// @brief Resident mesh type.
    using MeshType = Geometry::Mesh<ContextType>;

    /// @brief Execution space of the mesh context.
    using ExecutionSpace = typename MeshType::ExecutionSpace;

    /// @brief Memory space of the mesh context.
    using MemorySpace = typename MeshType::MemorySpace;

    /// @brief Tells whether a callable is a device predicate on polytopes.
    template <class Pred>
    static constexpr bool IsPredicate = std::is_invocable_r_v<bool, const Pred&,
      const ::Rodin::Kokkos::Geometry::PolytopeView&>;

    /**
     * @brief Projects a source on a grid function.
     * @param gf Grid function of a device space.
     * @param region Region on which to project.
     * @param fn Source to evaluate.
     * @param attrs Mesh attributes selecting the region; empty selects it
     * all.
     * @param select Device predicate selecting the polytopes.
     */
    template <class GridFunction, class Function,
      class Select = ::Rodin::Kokkos::Variational::SelectAll>
    static void project(GridFunction& gf, Geometry::Region region, const Function& fn,
      const FlatSet<Geometry::Attribute>& attrs, const Select& select = {})
    {
      const auto& fes = gf.getFiniteElementSpace();
      const auto& mesh = fes.getMesh();
      const size_t dimension = mesh.getDimension();
      if ((region == Geometry::Region::Boundary ||
            region == Geometry::Region::Interface) &&
        !mesh.getConnectivity().hasIncidence(dimension - 1, dimension))
      {
        Alert::MemberFunctionException(gf, __func__)
          << "Compute face-to-cell incidence before projecting on the boundary or an "
             "interface."
          << Alert::Raise;
      }
      const size_t d =
        ::Rodin::Kokkos::Variational::Internal::getRegionDimension(dimension, region);
      const auto& execution = mesh.getContext().getExecutionSpace();
      const auto dofs = fes.getDeviceDOFs();
      const Index stride = dofs.getStride();
      const Index vdim = fes.getVectorDimension();
      const auto attributes =
        ::Rodin::Kokkos::Variational::Internal::uploadAttributes(execution, attrs);
      // Zero marks the degrees of freedom outside the selection.
      ::Kokkos::View<Index*, MemorySpace> owners(
        ::Kokkos::view_alloc(execution, std::string("Rodin projection owners")),
        stride);
      ::Kokkos::View<Real*, MemorySpace> values(
        ::Kokkos::view_alloc(execution, std::string("Rodin projection values")),
        stride * vdim);
      evaluate(execution, mesh.getDeviceView(), region, d, attributes, select, owners,
        values, lower(fn, vdim), dofs, mesh.getPolytopeCount(d), vdim);

      const auto hostOwners = ::Kokkos::create_mirror_view(owners);
      const auto hostValues = ::Kokkos::create_mirror_view(values);
      ::Kokkos::deep_copy(execution, hostOwners, owners);
      ::Kokkos::deep_copy(execution, hostValues, values);
      execution.fence("Rodin projection values");
      if constexpr (requires { gf.flush(); })
        static_cast<const GridFunction&>(gf).flush();
      for (Index r = 0; r < vdim; ++r)
      {
        for (Index i = 0; i < stride; ++i)
        {
          if (hostOwners(i) != 0)
            gf[i + r * stride] = hostValues(i + r * stride);
        }
      }
      if constexpr (requires { gf.flush(); })
        gf.flush();
    }

    /**
     * @brief Lowers a source to a device view.
     * @param fn Source.
     * @param vdim Vector dimension of the space.
     * @return Device view providing getValue(const PointView&).
     */
    template <class Function>
    static auto lower(const Function& fn, [[maybe_unused]] Index vdim)
    {
      if constexpr (decltype(isFunction(static_cast<const Function*>(nullptr)))::value)
      {
        return lowerFunction(fn);
      }
      else if constexpr (std::is_arithmetic_v<Function>)
      {
        return ::Rodin::Kokkos::Variational::ConstantView<Real>{static_cast<Real>(fn)};
      }
      else if constexpr (FormLanguage::IsVectorRange<Function>::Value)
      {
        ::Rodin::Kokkos::Geometry::Vector value(vdim);
        for (Index r = 0; r < vdim; ++r)
          value(r) = fn(static_cast<std::uint8_t>(r));
        return ::Rodin::Kokkos::Variational::ConstantView<
          ::Rodin::Kokkos::Geometry::Vector>{value};
      }
      else
      {
        static_assert(std::is_invocable_v<const Function&,
                        const ::Rodin::Kokkos::Geometry::PointView&>,
          "A function projected on a Kokkos mesh must be a KOKKOS_LAMBDA taking a "
          "const Rodin::Kokkos::Geometry::PointView&.");
        return ::Rodin::Kokkos::Variational::CallableView<Function>{fn};
      }
    }

    /**
     * @brief Lowers a function expression to its device view.
     * @param fn Function expression.
     * @return Device view of the expression.
     */
    template <class Derived>
    static auto lowerFunction(const FunctionBase<Derived>& fn)
    {
      return ::Rodin::Kokkos::Variational::DeviceView<FunctionBase<Derived>>(fn);
    }

    /**
     * @brief Tells that a source is a function expression.
     * @param source Pointer to the source.
     * @return std::true_type, used only in unevaluated context.
     */
    template <class Derived>
    static std::true_type isFunction(const FunctionBase<Derived>* source);

    /**
     * @brief Tells that a source is not a function expression.
     * @param source Pointer to the source.
     * @return std::false_type, used only in unevaluated context.
     */
    static std::false_type isFunction(const void* source);

    /**
     * @brief Marks every degree of freedom of a selected polytope with one plus
     * the polytope index, keeping the largest.
     */
    template <class AttributeView, class Select, class DOFs>
    struct Owners
    {
      /// @brief Device mesh.
      ::Rodin::Kokkos::Geometry::MeshView view;

      /// @brief Region on which to project.
      Geometry::Region region;

      /// @brief Dimension of the polytopes of the region.
      Index d;

      /// @brief Selected attributes; empty selects the whole region.
      AttributeView attributes;

      /// @brief Device predicate selecting the polytopes.
      Select select;

      /**
       * @brief One plus the index of the polytope owning each scalar degree of
       * freedom, zero when none does.
       */
      ::Kokkos::View<Index*, MemorySpace> owners;

      /// @brief Device degrees of freedom.
      DOFs dofs;

      /**
       * @brief Marks the degrees of freedom of a polytope if it is selected.
       * @param i Index of the polytope of dimension d.
       */
      KOKKOS_INLINE_FUNCTION void operator()(Index i) const
      {
        const auto polytope = view.getPolytope(d, i);
        if (!::Rodin::Kokkos::Variational::isIntegrated(
              view, region, polytope, attributes) ||
          !select(polytope))
          return;
        const auto indices = dofs.bind(polytope);
        for (Index q = 0; q < indices.size(); ++q)
          ::Kokkos::atomic_max(&owners(indices[q]), i + 1);
      }
    };

    /**
     * @brief Evaluates the source at the nodes of the degrees of freedom a
     * selected polytope owns.
     */
    template <class AttributeView, class Select, class ValueView, class DOFs>
    struct Values
    {
      /// @brief Device mesh.
      ::Rodin::Kokkos::Geometry::MeshView view;

      /// @brief Region on which to project.
      Geometry::Region region;

      /// @brief Dimension of the polytopes of the region.
      Index d;

      /// @brief Selected attributes; empty selects the whole region.
      AttributeView attributes;

      /// @brief Device predicate selecting the polytopes.
      Select select;

      /**
       * @brief One plus the index of the polytope owning each scalar degree of
       * freedom, zero when none does.
       */
      ::Kokkos::View<Index*, MemorySpace> owners;

      /**
       * @brief Projected values, component r of the degree of freedom i at
       * i + r * stride.
       */
      ::Kokkos::View<Real*, MemorySpace> values;

      /// @brief Device view of the source.
      ValueView value;

      /// @brief Device degrees of freedom.
      DOFs dofs;

      /// @brief Offset between components.
      Index stride;

      /// @brief Vector dimension of the space.
      Index vdim;

      /**
       * @brief Evaluates the source at the nodes of the degrees of freedom
       * owned by a selected polytope.
       * @param i Index of the polytope of dimension d.
       */
      KOKKOS_INLINE_FUNCTION void operator()(Index i) const
      {
        const auto polytope = view.getPolytope(d, i);
        if (!::Rodin::Kokkos::Variational::isIntegrated(
              view, region, polytope, attributes) ||
          !select(polytope))
          return;
        const auto indices = dofs.bind(polytope);
        for (Index q = 0; q < indices.size(); ++q)
        {
          const Index dof = indices[q];
          if (owners(dof) != i + 1)
            continue;
          const auto p = polytope.getPoint(dofs.getNode(polytope, q));
          store(values, dof, stride, vdim, value.getValue(p));
        }
      }
    };

    /**
     * @brief Selects the degrees of freedom and evaluates the source.
     * @param execution Execution space instance.
     * @param view Device mesh.
     * @param region Region on which to project.
     * @param d Dimension of the polytopes of the region.
     * @param attributes Selected attributes; empty selects the whole region.
     * @param select Device predicate selecting the polytopes.
     * @param owners One plus the largest selected polytope of every scalar
     * degree of freedom, zero when none is.
     * @param values Value of component r of the degree of freedom i at
     * i + r * stride.
     * @param value Device view of the source.
     * @param dofs Device degrees of freedom.
     * @param count Number of polytopes of dimension d.
     * @param vdim Vector dimension of the space.
     *
     * The count is given by the host, since the topology of the device mesh
     * is not readable on the host.
     */
    template <class AttributeView, class Select, class ValueView, class DOFs>
    static void evaluate(const ExecutionSpace& execution,
      const ::Rodin::Kokkos::Geometry::MeshView& view, Geometry::Region region, Index d,
      const AttributeView& attributes, const Select& select,
      const ::Kokkos::View<Index*, MemorySpace>& owners,
      const ::Kokkos::View<Real*, MemorySpace>& values, const ValueView& value,
      const DOFs& dofs, Index count, Index vdim)
    {
      const Index stride = dofs.getStride();
      const auto policy =
        ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
          execution, 0, count);
      ::Kokkos::parallel_for("Rodin projection owners", policy,
        Owners<AttributeView, Select, DOFs>{
          view, region, d, attributes, select, owners, dofs});
      ::Kokkos::parallel_for("Rodin projection values", policy,
        Values<AttributeView, Select, ValueView, DOFs>{view, region, d, attributes,
          select, owners, values, value, dofs, stride, vdim});
    }

    /**
     * @brief Stores a scalar value of a degree of freedom.
     * @param values Projected values.
     * @param dof Scalar degree of freedom.
     * @param stride Offset between components.
     * @param vdim Vector dimension of the space.
     * @param v Value at the node.
     */
    KOKKOS_INLINE_FUNCTION static void store(
      const ::Kokkos::View<Real*, MemorySpace>& values, Index dof,
      [[maybe_unused]] Index stride, [[maybe_unused]] Index vdim, Real v)
    {
      values(dof) = v;
    }

    /**
     * @brief Stores the components of a vector value of a degree of freedom.
     * @param values Projected values, component r of the degree of freedom i
     * at i + r * stride.
     * @param dof Scalar degree of freedom.
     * @param stride Offset between components.
     * @param vdim Vector dimension of the space.
     * @param v Value at the node.
     */
    KOKKOS_INLINE_FUNCTION static void store(
      const ::Kokkos::View<Real*, MemorySpace>& values, Index dof, Index stride,
      Index vdim, const ::Rodin::Kokkos::Geometry::Vector& v)
    {
      for (Index r = 0; r < vdim; ++r)
        values(dof + r * stride) = v(r);
    }
  };
}

#endif
