/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Complete active element tangents agree with explicit component contraction.
#include <gtest/gtest.h>
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  namespace
  {
    template <size_t Degree>
    void checkActiveTangent(Geometry::Polytope::Type geometry, size_t dimension)
    {
      using namespace Variational;
      auto mesh =
        Geometry::LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      for (size_t i = 1; i <= dimension; ++i)
      {
        for (size_t j = 0; j <= dimension; ++j)
        {
          if (i != j)
          {
            mesh.getConnectivity().compute(i, j);
          }
        }
      }
      using Space = std::conditional_t<Degree == 1, P1<Math::SpatialVector<Real>>,
        H1<Degree, Math::SpatialVector<Real>>>;
      auto space = [&] {
        if constexpr (Degree == 1)
        {
          return Space(mesh, dimension);
        }
        else
        {
          return Space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
        }
      }();
      TrialFunction trial(space);
      TestFunction test(space);
      GridFunction displacement(space);
      Solid::ActiveFiberLaw::Parameters parameters;
      parameters.stiffness = 2;
      parameters.damping = 0.3;
      parameters.contractility = 0.8;
      parameters.crossBridgeStiffness = 0.5;
      parameters.initial.stiffness = 1.2;
      parameters.initial.stress = 0.4;
      const Solid::ActiveContraction prototype(
        Solid::NeoHookean(2, 3), Solid::ActiveFiberLaw(parameters));
      const auto cell = *mesh.getCell();
      const auto& element = space.getFiniteElement(cell.getDimension(), cell.getIndex());
      const size_t count = element.getCount();
      const auto& formula = QF::QuadratureFormula::get(
        2 * element.getDegree().value_or(2 * element.getOrder()), geometry);
      const auto& quadrature = cell.getQuadrature(formula);
      for (bool dynamic : {false, true})
      {
        const auto input = [dynamic, dimension](Solid::ConstitutivePoint& point) {
          Math::SpatialVector<Real> direction(dimension);
          for (size_t i = 0; i < dimension; ++i)
          {
            direction(i) = i + 1;
          }
          point.set<Solid::Tags::FiberDirection>(direction);
          point.set<Solid::Tags::ActiveExtension>(-0.01);
          if (dynamic)
          {
            point.set<Solid::Tags::TimeStep>(0.02);
            point.set<Solid::Tags::PreviousActiveExtension>(-0.015);
            point.set<Solid::Tags::PreviousActiveGamma>(1.2);
            point.set<Solid::Tags::PreviousActiveBeta>(0.8);
            point.set<Solid::Tags::ElectricalActivation>(0.4);
            point.set<Solid::Tags::ActiveGamma>(1.1);
            point.set<Solid::Tags::ActiveBeta>(0.7);
            point.set<Solid::Tags::PreviousFiberStrain>(0.015);
          }
        };
        Solid::InternalVirtualWork work(prototype, displacement);
        work.setInput(input);
        auto tangent = work.Tangent(trial, test);
        // Reuse the integrator after changing the displacement coefficients.
        for (Real amplitude : {0.001, 0.002})
        {
          for (size_t i = 0; i < space.getSize(); ++i)
          {
            displacement.getData()(i) = amplitude * (i % 11);
          }
          tangent.setPolytope(cell);
          Math::Matrix<Real> expected = Math::Matrix<Real>::Zero(count, count);
          auto law = prototype;
          auto stateGradient = Jacobian(displacement);
          auto trialGradient = Jacobian(trial);
          auto testGradient = Jacobian(test);
          trialGradient.setPolytope(cell);
          testGradient.setPolytope(cell);
          for (size_t q = 0; q < quadrature.getSize(); ++q)
          {
            const IntegrationPoint ip(quadrature, q);
            trialGradient.setIntegrationPoint(ip);
            testGradient.setIntegrationPoint(ip);
            Solid::KinematicState state(dimension);
            state.setDisplacementGradient(stateGradient.getValue(ip));
            Solid::ConstitutivePoint point(quadrature.getPoint(q), state);
            input(point);
            law.setConstitutivePoint(point);
            const Real weight = formula.getWeight(q) * ip.getPoint().getDistortion();
            for (size_t tr = 0; tr < count; ++tr)
            {
              Math::SpatialMatrix<Real> increment;
              law.getMaterialTangent(increment, trialGradient.getBasis(tr));
              for (size_t te = 0; te < count; ++te)
              {
                const auto gradient = testGradient.getBasis(te);
                Real value = 0;
                for (size_t row = 0; row < dimension; ++row)
                {
                  for (size_t col = 0; col < dimension; ++col)
                  {
                    value += increment(row, col) * gradient(row, col);
                  }
                }
                expected(te, tr) += weight * value;
              }
            }
          }
          Math::Matrix<Real> actual(count, count);
          for (size_t tr = 0; tr < count; ++tr)
          {
            for (size_t te = 0; te < count; ++te)
            {
              actual(te, tr) = tangent.integrate(tr, te);
            }
          }
          SCOPED_TRACE(::testing::Message()
            << "degree=" << Degree << " dimension=" << dimension
            << " geometry=" << static_cast<int>(geometry) << " dynamic=" << dynamic
            << " amplitude=" << amplitude);
          ASSERT_TRUE(actual.allFinite());
          ASSERT_GT(expected.norm(), 1);
          EXPECT_LE((actual - expected).norm(), 1e-12 * expected.norm());
        }
      }
    }

    template <size_t Degree>
    void checkGeometries()
    {
      for (auto geometry :
        {Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral})
      {
        checkActiveTangent<Degree>(geometry, 2);
      }
      for (auto geometry :
        {Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Hexahedron})
      {
        checkActiveTangent<Degree>(geometry, 3);
      }
    }

    TEST(SolidInternalVirtualWorkTangent, RejectsUnsupportedSpatialDimension)
    {
      using namespace Variational;
      auto mesh =
        Geometry::LocalMesh::UniformGrid(Geometry::Polytope::Type::Segment, {2});
      P1 space(mesh, size_t{1});
      TrialFunction trial(space);
      TestFunction test(space);
      GridFunction displacement(space);
      displacement.getData().setZero();
      Solid::InternalVirtualWork work(Solid::NeoHookean(2, 3), displacement);
      auto tangent = work.Tangent(trial, test);
      EXPECT_THROW(tangent.setPolytope(*mesh.getCell()), Alert::Exception);
    }

    TEST(SolidInternalVirtualWorkTangent, P1ComponentContraction)
    {
      checkGeometries<1>();
    }

    TEST(SolidInternalVirtualWorkTangent, P2ComponentContraction)
    {
      checkGeometries<2>();
    }

    TEST(SolidInternalVirtualWorkTangent, P3ComponentContraction)
    {
      checkGeometries<3>();
    }
  }
}
