/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_COMMON_H
#define KELVIN_BALL_COMMON_H

#include <array>
#include <initializer_list>

#include <Rodin/Geometry.h>
#ifdef RODIN_USE_UMFPACK
#include <Rodin/Solver/UMFPack.h>
#else
#include <Rodin/Solver/SparseLU.h>
#endif
#include <Rodin/Variational.h>

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  inline constexpr Attribute Obstacle = 2;
  inline constexpr Attribute Fluid = 3;
  inline constexpr Attribute Outer = 7;
  inline constexpr Attribute Gamma = 13;
  inline constexpr Attribute SigmaPlus = 31;
  inline constexpr Attribute SigmaMinus = 32;
  inline constexpr Attribute SigmaXYPlus = 33;
  inline constexpr Attribute SigmaXYMinus = 35;

  inline constexpr Real Mu = 1;
  inline constexpr Real DefaultNitschePenalty = 320;
  inline constexpr Real DefaultStabilizationFactor = 0.05;
  inline constexpr Real LinearResidualTolerance = 1e-8;
  inline constexpr size_t ChamberMultiplicity = 24;

#ifdef RODIN_USE_UMFPACK
  inline constexpr const char* DirectSolverName = "UMFPACK";
#else
  inline constexpr const char* DirectSolverName = "Eigen SparseLU";
#endif

  using Mesh = Geometry::Mesh<Context::Local>;
  using VelocitySpace = P1<Math::SpatialVector<Real>, Mesh>;
  using PressureSpace = P1<Real, Mesh>;

  struct RotationPair
  {
      RotationPair(
        Attribute slave, Attribute master, std::initializer_list<Real> coefficients);

      Attribute slave;
      Attribute master;
      Math::SpatialMatrix<Real> rotation;
  };

  extern const std::array<RotationPair, 2> RotationPairs;
  extern const FlatSet<Attribute> MasterCuts;

  Math::SpatialPoint centroid(const Mesh& mesh, const Polytope& face);

  void splitSelfPairedCut(Mesh& mesh);

  void prepare(Mesh& mesh);

  template <class Problem>
  void solveDirect(Problem& problem)
  {
#ifdef RODIN_USE_UMFPACK
    Solver::UMFPack(problem).solve();
#else
    Solver::SparseLU(problem).solve();
#endif
  }
}

#endif
