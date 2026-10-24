/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  class ModalVectorElement;
  class ModalVectorSpace;
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Tests::Unit::ModalVectorElement>
  {
    using ScalarType = Real;
    using RangeType = Math::SpatialVector<Real>;
  };

  template <>
  struct Traits<Tests::Unit::ModalVectorSpace>
  {
    using MeshType = Geometry::LocalMesh;
    using ScalarType = Real;
    using RangeType = Math::SpatialVector<Real>;
    using ElementType = Tests::Unit::ModalVectorElement;
  };
}

namespace Rodin::Tests::Unit
{
  /// @brief Coupled vector polynomial modes without nodes or a component DOF layout.
  class ModalVectorElement final
    : public Variational::FiniteElementBase<ModalVectorElement>
  {
    public:
      ModalVectorElement()
        : FiniteElementBase(Geometry::Polytope::Type::Triangle)
      {}

      struct Basis
      {
        size_t local;

        Math::SpatialVector<Real> operator()(const Math::SpatialPoint& point) const
        {
          const Real x = point(0), y = point(1);
          switch (local)
          {
            case 0:
              return {1, x};
            case 1:
              return {y, 1};
            default:
              return {x * y, x - y};
          }
        }
      };

      size_t getCount() const
      {
        return 3;
      }

      Basis getBasis(size_t local) const
      {
        return {local};
      }
  };

  /// @brief Permuted DOFs with contravariant/covariant maps and optional resampling.
  class ModalVectorSpace final
    : public Variational::FiniteElementSpace<Geometry::LocalMesh, ModalVectorSpace>
  {
    public:
      template <class Function>
      class Pushforward
      {
        public:
          Pushforward(Function function, bool covariant, bool resample)
            : m_function(function),
              m_covariant(covariant),
              m_resample(resample)
          {}

          Math::SpatialVector<Real> operator()(const Geometry::Point& point) const
          {
            auto reference = point.getReferenceCoordinates();
            if (m_resample)
              reference(0) += Real(0.125);
            const auto value = m_function(reference);
            if (m_covariant)
              return point.getJacobianInverse().transpose() * value;
            const auto& jacobian = point.getJacobian();
            const Real determinant =
              jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
            return (jacobian * value) / determinant;
          }

        private:
          Function m_function;
          bool m_covariant, m_resample;
      };

      ModalVectorSpace(const Geometry::LocalMesh& mesh, bool covariant, bool resample)
        : m_mesh(mesh),
          m_covariant(covariant),
          m_resample(resample)
      {}

      size_t getSize() const override
      {
        return 3;
      }

      size_t getVectorDimension() const override
      {
        return 2;
      }

      const Geometry::LocalMesh& getMesh() const override
      {
        return m_mesh.get();
      }

      const ModalVectorElement& getFiniteElement(size_t, Index) const
      {
        return m_element;
      }

      const IndexArray& getDOFs(size_t, Index) const override
      {
        static const IndexArray indices{{2, 0, 1}};
        return indices;
      }

      template <class Function>
      auto getPushforward(const std::pair<size_t, Index>&, Function&& function) const
      {
        return Pushforward<std::decay_t<Function>>(function, m_covariant, m_resample);
      }

    private:
      std::reference_wrapper<const Geometry::LocalMesh> m_mesh;
      ModalVectorElement m_element;
      bool m_covariant, m_resample;
  };

  /// @brief Field evaluation respects mode order, live coefficients and both Piola maps.
  TEST(Rodin_Variational_GridFunction, NonNodalVectorExpansionUsesMappingAndDOFOrder)
  {
    for (bool covariant : {false, true})
    {
      for (bool resample : {false, true})
      {
        SCOPED_TRACE(covariant);
        SCOPED_TRACE(resample);
        auto mesh = Geometry::LocalMesh::Builder()
                      .initialize(2)
                      .nodes(3)
                      .vertex({0, 0})
                      .vertex({2, 0})
                      .vertex({0.5, 3})
                      .polytope(Geometry::Polytope::Type::Triangle, {0, 1, 2})
                      .finalize();
        ModalVectorSpace space(mesh, covariant, resample);
        Variational::GridFunction field(space);
        field.getData() << 2, -3, 5;
        const auto cell = *mesh.getCell();
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
        for (size_t mutation = 0; mutation < 2; ++mutation)
        {
          const Geometry::PolytopeQuadrature quadrature(cell, formula);
          for (size_t sample = 0; sample < quadrature.getSize(); ++sample)
          {
            const Real x = formula.getPoint(sample)(0) + (resample ? Real(0.125) : 0);
            const Real y = formula.getPoint(sample)(1);
            // Local coefficients are global entries 2, 0, 1, not component pairs.
            const Real first = Real(5 + mutation) + 2 * y - 3 * x * y;
            const Real second = Real(5 + mutation) * x + 2 - 3 * (x - y);
            const Real scale = mutation == 0 ? 1 : 2;
            Math::SpatialVector<Real> expected(2);
            if (covariant)
            {
              expected(0) = first / (2 * scale);
              expected(1) = (-first / 12 + second / 3) / scale;
            }
            else
            {
              expected(0) = (2 * first + Real(0.5) * second) / (6 * scale);
              expected(1) = second / (2 * scale);
            }
            const Variational::IntegrationPoint certified(quadrature, sample);
            const Variational::IntegrationPoint provenance(
              quadrature.getPoint(sample), &formula, sample);
            EXPECT_NEAR((field.getValue(certified) - expected).norm(), 0, 1e-13);
            EXPECT_NEAR((field.getValue(provenance) - expected).norm(), 0, 1e-13);
            EXPECT_NEAR(
              (field.getValue(quadrature.getPoint(sample)) - expected).norm(), 0, 1e-13);
          }
          mesh.scale(2);
          field.getData()(2) += 1;
        }
      }
    }
  }

  /// @brief Coefficient callbacks may evict tables before the next reference read.
  TEST(Rodin_Variational_GridFunction, ReferenceExpansionSurvivesCoefficientCacheEviction)
  {
    const auto geometry = Geometry::Polytope::Type::Segment;
    const auto& formula = QF::PolytopeQuadratureFormula::get(3, geometry);
    const Variational::H1Element<2, Real> scalar(geometry);
    const auto check = [&](const auto& element) {
      using Range = typename std::decay_t<decltype(element)>::RangeType;
      const auto coefficient = [](size_t local) { return Real(1 + 2 * local) / 7; };
      Range expected, actual;
      element.evaluate(expected, coefficient, formula.getPoint(0));
      element.evaluate(
        actual,
        [&](size_t local) {
          // Nine distinct formulas exceed the bounded reference cache capacity.
          for (size_t order = 10; order < 19; ++order)
            scalar.getTabulation(QF::PolytopeQuadratureFormula::get(order, geometry));
          return coefficient(local);
        },
        formula, 0);
      if constexpr (std::is_same_v<Range, Real>)
        EXPECT_NEAR(actual, expected, 1e-13);
      else
        EXPECT_NEAR((actual - expected).norm(), 0, 1e-13);
    };
    check(scalar);
    check(Variational::H1Element<2, Math::SpatialVector<Real>>(geometry, 2));
    check(Variational::H1Element<2, Math::SpatialMatrix<Real>>(geometry, 2, 2));
  }

  /// @brief Certified scalar, vector and matrix contractions agree with direct bases.
  TEST(Rodin_Variational_GridFunction, ReferenceExpansionMatchesDirectBasis)
  {
    const auto check = []<size_t Degree, class Scalar>() {
      for (const auto geometry : {Geometry::Polytope::Type::Segment,
             Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
             Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Hexahedron,
             Geometry::Polytope::Type::Pyramid, Geometry::Polytope::Type::Wedge})
      {
        SCOPED_TRACE(static_cast<int>(geometry));
        SCOPED_TRACE(Degree);
        const auto& formula = QF::PolytopeQuadratureFormula::get(2, geometry);
        const auto coefficient = [](size_t local) {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(Real(1 + local) / 7, -Real(2 + local) / 11);
          else
            return Scalar(Real(1 + local) / 7);
        };
        const auto compare = [&](const auto& element) {
          using Range = typename std::decay_t<decltype(element)>::RangeType;
          for (size_t sample = 0; sample < formula.getSize(); ++sample)
          {
            Range expected, actual;
            element.evaluate(expected, coefficient, formula.getPoint(sample));
            element.evaluate(actual, coefficient, formula, sample);
            if constexpr (std::is_same_v<Range, Scalar>)
              EXPECT_NEAR(std::abs(actual - expected), 0, 1e-11);
            else
              EXPECT_NEAR((actual - expected).norm(), 0, 1e-11);
          }
        };
        compare(Variational::H1Element<Degree, Scalar>(geometry));
        compare(Variational::H1Element<Degree, Math::SpatialVector<Scalar>>(geometry, 2));
        compare(
          Variational::H1Element<Degree, Math::SpatialMatrix<Scalar>>(geometry, 2, 2));
      }
    };
    check.template operator()<2, Real>();
    check.template operator()<3, Complex>();
  }

}
