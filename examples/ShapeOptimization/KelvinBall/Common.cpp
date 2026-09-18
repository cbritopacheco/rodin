/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace KelvinBall
{
  const std::array<RotationPair, 2> RotationPairs{
    {{SigmaXYMinus, SigmaXYPlus, {0, 1, 0, 1, 0, 0, 0, 0, -1}},
      {SigmaMinus, SigmaPlus, {1, 0, 0, 0, 0, -1, 0, 1, 0}}}};

  const FlatSet<Attribute> MasterCuts{SigmaPlus, SigmaXYPlus};

  RotationPair::RotationPair(
    Attribute slave, Attribute master, std::initializer_list<Real> coefficients)
    : slave(slave),
      master(master),
      rotation(3, 3)
  {
    auto coefficient = coefficients.begin();
    for (size_t row = 0; row < 3; ++row)
      for (size_t column = 0; column < 3; ++column)
        rotation(row, column) = *coefficient++;
  }

  Math::SpatialPoint centroid(const Mesh& mesh, const Polytope& face)
  {
    Math::SpatialPoint result(3);
    result.setZero();
    for (const Index vertex : face.getVertices())
      result += mesh.getVertexCoordinates(vertex);
    return result / static_cast<Real>(face.getVertices().size());
  }

  void splitSelfPairedCut(Mesh& mesh)
  {
    const size_t faceDimension = mesh.getDimension() - 1;
    for (auto face = mesh.getPolytope(faceDimension); face; ++face)
    {
      if (face->getAttribute() == SigmaXYPlus && centroid(mesh, *face).z() < 0)
        mesh.setAttribute({faceDimension, face->getIndex()}, SigmaXYMinus);
    }
  }

  void prepare(Mesh& mesh)
  {
    splitSelfPairedCut(mesh);
    auto& connectivity = mesh.getConnectivity();
    connectivity.discover(3, 2);
    connectivity.discover(3, 1);
    connectivity.restrict(1, 0);
    connectivity.restrict(2, 0);
    connectivity.restrict(2, 3);
  }
}
