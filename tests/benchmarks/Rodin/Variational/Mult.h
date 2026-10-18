/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Public scalar-weighted basis-value production without pair contraction.
 */
#ifndef RODIN_RODIN_VARIATIONAL_MULT_H
#define RODIN_RODIN_VARIATIONAL_MULT_H

#include <benchmark/benchmark.h>
#include <string>
#include <type_traits>
#include <vector>
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0g.h"
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  inline bool basisScalingFailed = false;

  /// @brief Space constructors selected by the basis-production fixture.
  enum class BasisFamily
  {
    H1,
    P1,
    P0,
    P0g
  };

  /// @brief Public node whose complete quadrature/basis sweep is measured.
  enum class BasisScalingOperation
  {
    Leaf,
    Left,
    Right
  };

  /**
   * @brief Measures leaf and weighted basis production through public expression nodes.
   *
   * ## Mathematical foundation
   *
   * For each mapped sample @f$x_q@f$ and local basis index @f$a@f$, stores
   * @f$\phi_a(x_q)@f$, @f$\rho(x_q)\phi_a(x_q)@f$ or
   * @f$\phi_a(x_q)\rho(x_q)@f$, respectively. The scalar coefficient is
   * @f$\rho=s(2+x_0)@f$ over the reals and
   * @f$\rho=s(2+x_0+i/2)@f$ over the complex numbers. Point uses @f$x_0=0@f$.
   * Component values are stored in quadrature/basis/component order, with
   * row-major logical components for matrices. This is an output layout,
   * not an assumption about global DOF numbering or nodal interpolation.
   *
   * ## Measurement boundary
   *
   * One sweep includes integration-point rebinding, coefficient-cache binding,
   * public getBasis calls and active-value stores. Mesh/space/node construction,
   * mapped quadrature and public-binding priming are outside timing. Binding
   * priming alone does not guarantee reference-table readiness for matrix
   * leaves. Prepared-reference scopes explicitly construct the scalar table
   * before priming and require its presence immediately before and after the
   * measured loop. Preparation is not included in their reported consumption
   * time. No basis
   * pairs, Frobenius contractions, weighted integration, scatter or solving are
   * included. The leaf control is a separate operation, not a timing to subtract
   * to obtain coefficient or allocation cost. All node borrows refer to one
   * owned integration-point context that outlives the nodes.
   *
   * ## Numerical qualification
   *
   * Returned logical vector and matrix dimensions are checked outside timing.
   * Every output is compared exactly with the public unweighted basis and the
   * corresponding scalar multiplication in its original operand order. This
   * isolates materialisation, not independent basis accuracy. Changed output
   * must be rejected. Coefficient scales @f$1\to2\to1@f$ exercise rebinding
   * without changing logical sample identities. Only one flattened active-value
   * output is retained; a second full basis-value table is not constructed.
   *
   * @tparam K H1 degree, one for P1, or zero for the constant families.
   * @tparam Range Scalar, spatial-vector or spatial-matrix field range.
   * @tparam Family Space constructor family.
   */
  template <size_t K, class Range, BasisFamily Family>
  class BasisScaling
  {
    public:
      using Scalar = typename FormLanguage::Traits<Range>::ScalarType;

      /**
       * @brief Measures one qualified sweep of a selected public node.
       * @tparam PreparedReference Whether an H1 scalar reference table is prepared.
       * @param state Logical rows, columns, quadrature order and iteration state.
       * @param geometry Reference-cell family, including Point.
       * @param operation Leaf or coefficient operand order.
       */
      template <bool PreparedReference = false>
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        BasisScalingOperation operation)
      {
        static_assert(!PreparedReference || Family == BasisFamily::H1,
          "Prepared reference scopes require an H1 element");
        using namespace Geometry;
        using namespace Variational;
        const size_t rows = state.range(0), columns = state.range(1);
        const size_t components = rows * columns, order = state.range(2);
        auto mesh = [&]() {
          if (geometry == Polytope::Type::Point)
            return LocalMesh::Builder()
              .initialize(3)
              .nodes(1)
              .vertex({0, 0, 0})
              .finalize();
          return Convergence::UniformGrid(geometry).makeMesh(2);
        }();
        for (size_t d = 1; d <= mesh.getDimension(); ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        using Mesh = decltype(mesh);
        using Space = std::conditional_t<Family == BasisFamily::H1, H1<K, Range, Mesh>,
          std::conditional_t<Family == BasisFamily::P1, P1<Range, Mesh>,
            std::conditional_t<Family == BasisFamily::P0, P0<Range, Mesh>,
              P0g<Range, Mesh>>>>;
        auto space = [&]() {
          if constexpr (Family == BasisFamily::H1)
          {
            if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              return Space(std::integral_constant<size_t, K>{}, mesh, rows, columns);
            else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              return Space(std::integral_constant<size_t, K>{}, mesh, rows);
            else
              return Space(std::integral_constant<size_t, K>{}, mesh);
          }
          else if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
            return Space(mesh, rows, columns);
          else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return Space(mesh, rows);
          else
            return Space(mesh);
        }();
        const auto cell = *mesh.getCell();
        const auto& formula = QF::PolytopeQuadratureFormula::get(order, geometry);
        const PolytopeQuadrature quadrature(cell, formula);
        IntegrationPoint ip(quadrature, 0);
        TrialFunction u(space);
        auto& leaf = u.getDerived();
        Real scale = 1;
        const auto coefficient = [&]() {
          const auto evaluate = [&scale](const Point& point) {
            const Real x = point.getDimension() == 0 ? Real(0) : point.x();
            if constexpr (std::is_same_v<Scalar, Complex>)
              return scale * Complex(Real(2) + x, Real(0.5));
            else
              return scale * (Real(2) + x);
          };
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(evaluate);
          else
            return RealFunction(evaluate);
        }();
        auto left = coefficient * u;
        auto right = u * coefficient;
        const size_t count = leaf.getDOFs(cell), samples = quadrature.getSize();
        if (components == 0 || count == 0 || count % components != 0 || samples == 0)
        {
          basisScalingFailed = true;
          state.SkipWithError("Basis/sample dimensions do not match the selected range");
          return;
        }
        std::vector<Scalar> output(samples * count * components);
        const auto component = [&](const auto& value, size_t c) -> Scalar {
          if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
            return value(c / columns, c % columns);
          else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return value(c);
          else
            return value;
        };
        const auto sweepNode = [&](auto& node) {
          for (size_t qp = 0; qp < samples; ++qp)
          {
            ip = IntegrationPoint(quadrature, qp);
            node.setIntegrationPoint(ip);
            for (size_t a = 0; a < count; ++a)
            {
              decltype(auto) value = node.getBasis(a);
              for (size_t c = 0; c < components; ++c)
                output[(qp * count + a) * components + c] = component(value, c);
            }
          }
        };
        const auto sweep = [&]() {
          switch (operation)
          {
            case BasisScalingOperation::Leaf:
              sweepNode(leaf);
              break;
            case BasisScalingOperation::Left:
              sweepNode(left);
              break;
            case BasisScalingOperation::Right:
              sweepNode(right);
              break;
          }
        };
        const auto validNode = [&](auto& node) {
          for (size_t qp = 0; qp < samples; ++qp)
          {
            ip = IntegrationPoint(quadrature, qp);
            leaf.setIntegrationPoint(ip);
            node.setIntegrationPoint(ip);
            const Scalar weight = coefficient.getValue(ip);
            for (size_t a = 0; a < count; ++a)
            {
              decltype(auto) actual = node.getBasis(a);
              if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              {
                if (actual.rows() != rows || actual.cols() != columns)
                  return false;
              }
              else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              {
                if (actual.size() != rows)
                  return false;
              }
              decltype(auto) basis = leaf.getBasis(a);
              for (size_t c = 0; c < components; ++c)
              {
                const Scalar value = component(basis, c);
                const Scalar expected = operation == BasisScalingOperation::Leaf ? value
                  : operation == BasisScalingOperation::Left ? weight * value
                                                             : value * weight;
                if (output[(qp * count + a) * components + c] != expected ||
                  component(actual, c) != expected)
                  return false;
              }
            }
          }
          return true;
        };
        const auto valid = [&]() {
          switch (operation)
          {
            case BasisScalingOperation::Leaf:
              return validNode(leaf);
            case BasisScalingOperation::Left:
              return validNode(left);
            case BasisScalingOperation::Right:
              return validNode(right);
          }
          return false;
        };
        const auto fail = [&](const char* message) {
          basisScalingFailed = true;
          state.SkipWithError(message);
        };
        if constexpr (PreparedReference)
        {
          const H1Element<K, Scalar> element(geometry);
          const auto& table = element.getTabulation(formula);
          if (table.nqp != samples || table.ndof != count / components ||
            table.dim != Polytope::Traits(geometry).getDimension())
          {
            fail("Prepared reference dimensions differ from the selected space");
            return;
          }
        }
        const auto referencePresent = [&]() {
          if constexpr (Family == BasisFamily::H1)
          {
            const H1Element<K, Scalar> element(geometry);
            return bool(element.findTabulation(formula.getCacheIdentity()));
          }
          else
            return false;
        };
        sweep();
        if (!valid())
        {
          fail("Public basis values differ from ordered component multiplication");
          return;
        }
        output.front() += Scalar(1);
        if (valid())
        {
          fail("Basis-value oracle accepted changed output");
          return;
        }
        for (Real value : {Real(2), Real(1)})
        {
          scale = value;
          sweep();
          if (!valid())
          {
            fail("Coefficient rebinding changed ordered basis values");
            return;
          }
        }
        const bool readyBefore = referencePresent();
        if constexpr (PreparedReference)
        {
          if (!readyBefore)
          {
            fail("Prepared reference table is absent before timing");
            return;
          }
        }
        for (auto _ : state)
        {
          sweep();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        const bool readyAfter = referencePresent();
        if constexpr (PreparedReference)
        {
          if (!readyAfter)
          {
            fail("Prepared reference table is absent after timing");
            return;
          }
        }
        if (!valid())
          fail("Timed basis sweep changed its qualified output");
        state.counters["reference_prepared_protocol"] = PreparedReference;
        state.counters["reference_table_present_before"] = readyBefore;
        state.counters["reference_table_present_after"] = readyAfter;
        state.counters["cells"] = 1;
        state.counters["fixture_cells"] = mesh.getPolytopeCount(mesh.getDimension());
        state.counters["dimension"] = mesh.getDimension();
        state.counters["ambient_dimension"] = mesh.getSpaceDimension();
        state.counters["global_dofs"] = space.getSize();
        state.counters["local_basis"] = count;
        state.counters["scalar_basis"] = count / components;
        state.counters["field_rows"] = rows;
        state.counters["field_columns"] = columns;
        state.counters["components"] = components;
        state.counters["degree"] = K;
        state.counters["quadrature_order"] = order;
        state.counters["quadrature_points"] = samples;
        state.counters["basis_values"] = samples * count;
        state.counters["stored_scalar_values"] = output.size();
        state.counters["retained_output_bytes"] = output.size() * sizeof(Scalar);
      }

      /// @brief Registers ranges and sample orders without a spatial refinement claim.
      template <bool PreparedReference = false>
      static void registerCases()
      {
        static_assert(!PreparedReference || Family == BasisFamily::H1,
          "Prepared reference scopes require an H1 element");
        using Geometry::Polytope;
        const std::string family = Family == BasisFamily::H1 ? "H1/P" + std::to_string(K)
          : Family == BasisFamily::P1                        ? "P1"
          : Family == BasisFamily::P0                        ? "P0"
                                                             : "P0g";
        const std::string range = FormLanguage::IsMatrixRange<Range>::Value ? "Matrix"
          : FormLanguage::IsVectorRange<Range>::Value                       ? "Vector"
                                                                            : "Scalar";
        for (auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (auto operation : {BasisScalingOperation::Leaf, BasisScalingOperation::Left,
                 BasisScalingOperation::Right})
          {
            const std::string name =
              std::string(
                PreparedReference ? "BasisScaling/PreparedReference/" : "BasisScaling/") +
              (operation == BasisScalingOperation::Leaf      ? "Leaf/"
                  : operation == BasisScalingOperation::Left ? "Left/"
                                                             : "Right/") +
              family + "/" + (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
              range + "/" +
              (geometry == Polytope::Type::Point
                  ? "Point"
                  : std::string(Convergence::UniformGrid::getGeometryName(geometry)));
            auto* entry = benchmark::RegisterBenchmark(
              name.c_str(), [geometry, operation](auto& state) {
                BasisScaling::template run<PreparedReference>(state, geometry, operation);
              });
            for (int64_t order : {2, 4, 8})
            {
              if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              {
                entry->Args({1, 2, order});
                entry->Args({2, 3, order});
              }
              else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              {
                for (int64_t rows : {1, 2, 3})
                  entry->Args({rows, 1, order});
              }
              else
                entry->Args({1, 1, order});
            }
            entry->UseRealTime();
          }
        }
      }
  };

  template <size_t K, BasisFamily Family, bool PreparedReference = false>
  void registerBasisScaling()
  {
    BasisScaling<K, Real, Family>::template registerCases<PreparedReference>();
    BasisScaling<K, Complex, Family>::template registerCases<PreparedReference>();
    BasisScaling<K, Math::SpatialVector<Real>,
      Family>::template registerCases<PreparedReference>();
    BasisScaling<K, Math::SpatialVector<Complex>,
      Family>::template registerCases<PreparedReference>();
    BasisScaling<K, Math::SpatialMatrix<Real>,
      Family>::template registerCases<PreparedReference>();
    BasisScaling<K, Math::SpatialMatrix<Complex>,
      Family>::template registerCases<PreparedReference>();
  }

  /// @brief Registers H1 degrees one through three in one compilation unit.
  void registerBasisScalingH1();
  /// @brief Registers H1 degrees four through six in one compilation unit.
  void registerBasisScalingHigherOrder();
  /// @brief Registers P1, P0 and P0g in one compilation unit.
  void registerBasisScalingFixedDegree();
}

#endif // RODIN_RODIN_VARIATIONAL_MULT_H
