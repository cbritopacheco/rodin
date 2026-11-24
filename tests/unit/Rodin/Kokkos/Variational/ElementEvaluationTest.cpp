/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

#include <array>
#include <cstdlib>
#include <functional>
#include <type_traits>

namespace Rodin::Tests::Unit
{
  template <Index Capacity>
  class CoupledSpace;

  /// @brief Metadata for a synthetic coupled vector basis, without nodal functionals.
  struct CoupledElement
  {
    Optional<size_t> getDegree() const
    {
      return 2;
    }
    size_t getCount() const
    {
      return 3;
    }
  };
}

namespace Rodin::FormLanguage
{
  template <Index Capacity>
  struct Traits<Tests::Unit::CoupledSpace<Capacity>>
  {
    using ScalarType = Real;
    using RangeType = Math::SpatialVector<Real>;
    using MeshType = Tests::Unit::DeviceMesh;
    using ElementType = Tests::Unit::CoupledElement;
  };
}

namespace Rodin::Tests::Unit
{
  /**
   * @brief Synthetic intrinsic vector generators with optional orientation mixing.
   * No scalar gradient, component stride or nodal interpolation is provided.
   * @tparam Capacity Bound forcing either thread or team integration in tests.
   */
  template <Index Capacity>
  struct CoupledDOFs
  {
    bool covariant;
    bool oriented;
    ::Kokkos::Array<Index, 3> indices{2, 0, 1};

    KOKKOS_INLINE_FUNCTION static constexpr Index getCapacity(Index)
    {
      return Capacity;
    }

    KOKKOS_INLINE_FUNCTION auto bind(
      const ::Rodin::Kokkos::View<Geometry::Polytope>&) const
    {
      return ::Rodin::Kokkos::View<IndexVector>{indices.data(), 3};
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getReference(const Point& ip, Index q) const
    {
      const auto& r = ip.getPoint().getReferenceCoordinates();
      DeviceGeometry::Vector result(2);
      if (q == 0)
      {
        result(0) = 1;
        result(1) = r(0);
      }
      else if (q == 1)
      {
        result(0) = r(1) - (oriented ? 1 : 0);
        result(1) = 1 - (oriented ? r(0) : 0);
      }
      else
      {
        const Real sign = oriented ? -1 : 1;
        result(0) = sign * r(0) * r(1);
        result(1) = sign * (r(0) - r(1));
      }
      return result;
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getMap(const Point& ip) const
    {
      const auto point = ip.getPoint();
      const auto J = point.getJacobian();
      const auto inverse = point.getJacobianInverse();
      const Real det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
      DeviceGeometry::Matrix map(2, 2);
      for (Index i = 0; i < 2; ++i)
      {
        for (Index j = 0; j < 2; ++j)
          map(i, j) = covariant ? inverse(j, i) : J(i, j) / det;
      }
      return map;
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getBasis(const Point& ip, Index q) const
    {
      const auto ref = getReference(ip, q);
      const auto map = getMap(ip);
      DeviceGeometry::Vector result(2);
      for (Index i = 0; i < 2; ++i)
      {
        for (Index j = 0; j < 2; ++j)
          result(i) += map(i, j) * ref(j);
      }
      return result;
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION auto getJacobian(const Point& ip, Index q) const
    {
      const auto& r = ip.getPoint().getReferenceCoordinates();
      DeviceGeometry::Matrix reference(2, 2);
      if (q == 0)
        reference(1, 0) = 1;
      else if (q == 1)
      {
        reference(0, 1) = 1;
        reference(1, 0) = oriented ? -1 : 0;
      }
      else
      {
        const Real sign = oriented ? -1 : 1;
        reference(0, 0) = sign * r(1);
        reference(0, 1) = sign * r(0);
        reference(1, 0) = sign;
        reference(1, 1) = -sign;
      }
      const auto map = getMap(ip);
      const auto inverse = ip.getPoint().getJacobianInverse();
      DeviceGeometry::Matrix result(2, 2);
      for (Index i = 0; i < 2; ++i)
      {
        for (Index j = 0; j < 2; ++j)
        {
          for (Index a = 0; a < 2; ++a)
          {
            for (Index b = 0; b < 2; ++b)
              result(i, j) += map(i, a) * reference(a, b) * inverse(b, j);
          }
        }
      }
      return result;
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getDivergence(const Point& ip, Index q) const
    {
      const auto J = getJacobian(ip, q);
      return J(0, 0) + J(1, 1);
    }

    template <class Point>
    KOKKOS_INLINE_FUNCTION Real getCurl(const Point& ip, Index q) const
    {
      const auto J = getJacobian(ip, q);
      return J(1, 0) - J(0, 1);
    }
  };

  /// @brief An intrinsic family may expose divergence without gradients or a Jacobian.
  struct DivergenceOnlyDOFs : CoupledDOFs<3>
  {
    template <class Point>
    void getJacobian(const Point&, Index) const = delete;

    template <class Point>
    void getCurl(const Point&, Index) const = delete;
  };

  struct DivergenceOnlySpace
  {
    using RangeType = Math::SpatialVector<Real>;
    using DOFsType = DivergenceOnlyDOFs;
    DOFsType getDOFsView() const
    {
      return {};
    }
    size_t getVectorDimension() const
    {
      return 2;
    }
  };

  template <class Element>
  concept HasElementGradient = requires(const Element& element,
    const ::Rodin::Kokkos::Variational::IntegrationPointView& ip) {
    element.getGradient(ip, 0);
  };

  template <class Element>
  concept HasElementJacobian = requires(const Element& element,
    const ::Rodin::Kokkos::Variational::IntegrationPointView& ip) {
    element.getJacobian(ip, 0, 0);
  };

  template <class Element>
  concept HasElementDivergence = requires(const Element& element,
    const ::Rodin::Kokkos::Variational::IntegrationPointView& ip) {
    element.getDivergence(ip, 0, 0);
  };

  template <class Element>
  concept HasElementCurl = requires(const Element& element,
    const ::Rodin::Kokkos::Variational::IntegrationPointView& ip) {
    element.getCurl(ip, 0);
  };

  using DivergenceElement = ::Rodin::Kokkos::View<DivergenceOnlySpace>;
  static_assert(HasElementDivergence<DivergenceElement>);
  static_assert(!HasElementGradient<DivergenceElement>);
  static_assert(!HasElementJacobian<DivergenceElement>);
  static_assert(!HasElementCurl<DivergenceElement>);
  using ScalarElement = ::Rodin::Kokkos::View<Variational::P1<Real, DeviceMesh>>;
  static_assert(HasElementGradient<ScalarElement>);
  static_assert(!HasElementDivergence<ScalarElement>);
  static_assert(!HasElementCurl<ScalarElement>);

  /// @brief Test-only space of three coupled modes, not an RT or Nedelec family.
  template <Index Capacity>
  class CoupledSpace final
    : public Variational::FiniteElementSpace<DeviceMesh, CoupledSpace<Capacity>>
  {
    public:
      using MeshType = DeviceMesh;
      using ScalarType = Real;
      using RangeType = Math::SpatialVector<Real>;
      using ElementType = CoupledElement;
      using DOFsType = CoupledDOFs<Capacity>;

      CoupledSpace(const DeviceMesh& mesh, bool covariant, bool oriented)
        : m_mesh(mesh),
          m_dofs{covariant, oriented}
      {}

      size_t getSize() const override
      {
        return 3;
      }
      size_t getVectorDimension() const override
      {
        return 2;
      }
      const DeviceMesh& getMesh() const override
      {
        return m_mesh.get();
      }
      const ElementType& getFiniteElement(size_t, Index) const
      {
        return m_element;
      }
      DOFsType getDOFsView() const
      {
        return m_dofs;
      }
      const IndexArray& getDOFs(size_t, Index) const override
      {
        static const IndexArray indices{{2, 0, 1}};
        return indices;
      }

    private:
      std::reference_wrapper<const DeviceMesh> m_mesh;
      DOFsType m_dofs;
      CoupledElement m_element;
  };
}

namespace Rodin::Variational
{
  /// @brief Stateless shape leaf for the synthetic coupled-vector test space.
  template <class Derived, Index Capacity, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, Tests::Unit::CoupledSpace<Capacity>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, Tests::Unit::CoupledSpace<Capacity>, Space>,
        Tests::Unit::CoupledSpace<Capacity>, Space>
  {
    public:
      using FESType = Tests::Unit::CoupledSpace<Capacity>;
      using Parent = ShapeFunctionBase<ShapeFunction, FESType, Space>;
      explicit ShapeFunction(const FESType& fes)
        : Parent(fes)
      {}
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 2;
      }
      size_t getDOFs(const Geometry::Polytope&) const
      {
        return 3;
      }
      Math::SpatialVector<Real> getBasis(size_t) const
      {
        Alert::Exception() << "The synthetic fixture supplies device evaluation only."
                           << Alert::Raise;
        std::abort();
      }
      ShapeFunction& setPolytope(const Geometry::Polytope&)
      {
        return *this;
      }
      ShapeFunction& setIntegrationPoint(const IntegrationPoint&)
      {
        return *this;
      }
      const IntegrationPoint& getIntegrationPoint() const
      {
        Alert::Exception() << "The synthetic fixture has no host point binding."
                           << Alert::Raise;
        std::abort();
      }
      ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }
  };
}

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /**
   * @brief Evaluates and assembles complete coupled modes with one coefficient each.
   * @tparam Capacity Forces the existing thread or team kernel without changing the basis.
   */
  template <Index Capacity>
  void checkCoupledEvaluation()
  {
    const auto source = HostMesh::Builder()
                          .initialize(2)
                          .nodes(3)
                          .vertex({0, 0})
                          .vertex({2, 0})
                          .vertex({0.5, 3})
                          .polytope(Polytope::Type::Triangle, {0, 1, 2})
                          .finalize();
    const DeviceMesh mesh(TestContext{}, source);
    const auto& formula = QF::QuadratureFormula::get(4, Polytope::Type::Triangle);
    for (const bool covariant : {false, true})
    {
      for (const bool oriented : {false, true})
      {
        CoupledSpace<Capacity> fes(mesh, covariant, oriented);
        using Element = ::Rodin::Kokkos::View<decltype(fes)>;
        static_assert(!Element::Componentwise);
        static_assert(Element::ComponentCapacity == 1);
        const Element element(fes);
        EXPECT_EQ(element.getVectorDimension(), 2u);
        EXPECT_EQ(element.getComponentCount(), 1u);
        GridFunction field(fes);
        Math::Vector<Real> coefficients(3);
        coefficients << 2, -3, 5;
        field.getData().upload(coefficients);
        const ::Rodin::Kokkos::Variational::GridFunctionView<decltype(fes)> value(
          fes, field.getData());
        const auto meshView = mesh.getView();
        ::Kokkos::View<Real**, TestMemorySpace> output("coupled evaluation", 3, 18);
        ::Kokkos::parallel_for(
          "evaluate coupled basis", ::Kokkos::RangePolicy<TestExecutionSpace>(0, 3),
          KOKKOS_LAMBDA(Index q) {
            DeviceGeometry::Vector r(2);
            r(0) = Real(0.2);
            r(1) = Real(0.3);
            const auto polytope = meshView.getPolytope(2, 0);
            const ::Rodin::Kokkos::Variational::IntegrationPointView ip(
              polytope.getPoint(r));
            const auto basis = element.getBasis(ip, q);
            const auto jacobian = element.getJacobian(ip, q, 0);
            const auto fieldValue = value.interpolateVector(ip);
            output(q, 0) = element.getCount(polytope);
            output(q, 1) = element.getGlobalIndex(polytope, q);
            output(q, 2) = basis(0);
            output(q, 3) = basis(1);
            output(q, 4) = element.getDivergence(ip, q, 0);
            output(q, 5) = element.getCurl(ip, q);
            output(q, 6) = fieldValue(0);
            output(q, 7) = fieldValue(1);
            output(q, 8) = value.interpolateDivergence(ip);
            const auto fieldJacobian = value.interpolateJacobian(ip);
            for (Index i = 0; i < 2; ++i)
            {
              for (Index j = 0; j < 2; ++j)
              {
                output(q, 9 + 2 * i + j) = jacobian(i, j);
                output(q, 14 + 2 * i + j) = fieldJacobian(i, j);
              }
            }
            output(q, 13) = value.interpolate(ip, 1);
          });
        const auto actual =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, output);
        Math::Matrix<Real> J(2, 2);
        J << 2, 0.5, 0, 3;
        const Math::Matrix<Real> inverse = J.inverse();
        Math::Matrix<Real> map(2, 2);
        if (covariant)
          map = inverse.transpose();
        else
          map = J / J.determinant();
        const auto reference = [&](Real x, Real y) {
          Math::Matrix<Real> basis(2, 3);
          basis << 1, y, x * y, x, 1, x - y;
          if (oriented)
          {
            basis.col(1) -= basis.col(0);
            basis.col(2) *= -1;
          }
          return (map * basis).eval();
        };
        Math::Matrix<Real> expected = reference(Real(0.2), Real(0.3));
        const std::array<Index, 3> indices{2, 0, 1};
        Math::Vector<Real> local(3);
        for (Index q = 0; q < 3; ++q)
          local(q) = coefficients(indices[q]);
        const Math::Vector<Real> fieldValue = expected * local;
        Real divergence = 0;
        Math::Matrix<Real> fieldJacobian = Math::Matrix<Real>::Zero(2, 2);
        for (Index q = 0; q < 3; ++q)
        {
          Math::Matrix<Real> derivative(2, 2);
          if (q == 0)
            derivative << 0, 0, 1, 0;
          else if (q == 1)
            derivative << 0, 1, oriented ? -1 : 0, 0;
          else
          {
            derivative << Real(0.3), Real(0.2), 1, -1;
            if (oriented)
              derivative *= -1;
          }
          const Math::Matrix<Real> physical = map * derivative * inverse;
          EXPECT_EQ(actual(q, 0), 3);
          EXPECT_EQ(actual(q, 1), indices[q]);
          EXPECT_NEAR(actual(q, 2), expected(0, q), 1e-13);
          EXPECT_NEAR(actual(q, 3), expected(1, q), 1e-13);
          EXPECT_NEAR(actual(q, 4), physical.trace(), 1e-13);
          EXPECT_NEAR(actual(q, 5), physical(1, 0) - physical(0, 1), 1e-13);
          EXPECT_NEAR(actual(q, 6), fieldValue(0), 1e-13);
          EXPECT_NEAR(actual(q, 7), fieldValue(1), 1e-13);
          for (Index i = 0; i < 2; ++i)
          {
            for (Index j = 0; j < 2; ++j)
              EXPECT_NEAR(actual(q, 9 + 2 * i + j), physical(i, j), 1e-13);
          }
          EXPECT_NEAR(actual(q, 13), fieldValue(1), 1e-13);
          divergence += physical.trace() * local(q);
          fieldJacobian += physical * local(q);
        }
        EXPECT_NEAR(actual(0, 8), divergence, 1e-13);
        for (Index i = 0; i < 2; ++i)
        {
          for (Index j = 0; j < 2; ++j)
            EXPECT_NEAR(actual(0, 14 + 2 * i + j), fieldJacobian(i, j), 1e-13);
        }
        TrialFunction u(fes);
        TestFunction v(fes);
        BilinearForm mass(u, v);
        auto integral = Integral(u, v);
        EXPECT_FALSE(integral.isBlockDiagonal());
        CoupledSpace<Capacity> other(mesh, covariant, oriented);
        TrialFunction otherU(other);
        TestFunction otherV(other);
        const auto otherIntegral = Integral(otherU, otherV);
        EXPECT_FALSE(integral.getLayout() == otherIntegral.getLayout());
        integral.setOrder(4);
        mass = integral;
        mass.assemble();
        Math::Matrix<Real> referenceMass = Math::Matrix<Real>::Zero(3, 3);
        for (Index q = 0; q < formula.getSize(); ++q)
        {
          const auto& r = formula.getPoint(q);
          const auto basis = reference(r(0), r(1));
          const Math::Matrix<Real> contribution =
            formula.getWeight(q) * J.determinant() * basis.transpose() * basis;
          for (Index i = 0; i < 3; ++i)
          {
            for (Index j = 0; j < 3; ++j)
              referenceMass(indices[i], indices[j]) += contribution(i, j);
          }
        }
        const auto assembled = mass.getOperator().download();
        ASSERT_EQ(assembled.rows(), 3);
        ASSERT_EQ(assembled.cols(), 3);
        EXPECT_LE((Math::Matrix<Real>(assembled) - referenceMass).norm(), 1e-12);
        LinearForm rhs(v);
        auto load = Integral(field, v);
        load.setOrder(4);
        rhs = load;
        rhs.assemble();
        EXPECT_LE(
          (rhs.getVector().download() - referenceMass * coefficients).norm(), 1e-12);
      }
    }
  }

  /// @brief Intrinsic bases preserve value dimension, mapped derivatives and one coefficient per mode.
  TEST(Kokkos_Variational_Element, CoupledVectorEvaluationAndThreadAssembly)
  {
    checkCoupledEvaluation<3>();
  }

  /// @brief The team kernel consumes the same intrinsic element contract as the thread kernel.
  TEST(Kokkos_Variational_Element, CoupledVectorTeamAssembly)
  {
    checkCoupledEvaluation<100>();
  }

  /// @brief Complete coefficient indices agree with the existing componentwise space API.
  template <class FES>
  void checkComponentwiseElement(const FES& fes)
  {
    using Element = ::Rodin::Kokkos::View<FES>;
    static_assert(Element::Componentwise);
    const Element element(fes);
    const auto& expected = fes.getDOFs(2, 0);
    const auto mesh = fes.getMesh().getView();
    const Index count = expected.size();
    ::Kokkos::View<Index*, TestMemorySpace> indices(
      "complete element indices", count + 1);
    ::Kokkos::parallel_for(
      "read complete element indices",
      ::Kokkos::RangePolicy<TestExecutionSpace>(0, count), KOKKOS_LAMBDA(Index i) {
        const auto polytope = mesh.getPolytope(2, 0);
        indices(i) = element.getGlobalIndex(polytope, i);
        if (i == 0)
          indices(count) = element.getCount(polytope);
      });
    const auto actual =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, indices);
    EXPECT_EQ(actual(count), count);
    for (Index i = 0; i < count; ++i)
      EXPECT_EQ(actual(i), expected(i));
  }

  /// @brief Scalar and vector P1/H1 spaces retain their existing coefficient numbering.
  TEST(Kokkos_Variational_Element, ComponentwiseNumberingIsPreserved)
  {
    DeviceMesh mesh(TestContext{}, makeTriangleMesh());
    mesh.getConnectivity().compute(2, 1);
    P1 scalar(mesh);
    P1 vector(mesh, 2);
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    H1 cubic(std::integral_constant<size_t, 3>{}, mesh, 2);
    checkComponentwiseElement(scalar);
    checkComponentwiseElement(vector);
    checkComponentwiseElement(quadratic);
    checkComponentwiseElement(cubic);
  }

}
