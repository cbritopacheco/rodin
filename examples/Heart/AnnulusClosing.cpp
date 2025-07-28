/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <stdexcept>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <Eigen/Eigenvalues>

#include <Rodin/Geometry.h>
#include <Rodin/Math.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static constexpr Attribute AnnulusLabel = 54;

using Vector2 = Rodin::Math::Vector2<Real>;
using Vector3 = Rodin::Math::Vector3<Real>;

namespace
{
  struct TriangleLocal { int i0, i1, i2; };

  bool isCCW(const Vector2& a, const Vector2& b, const Vector2& c)
  {
    return ((b.x() - a.x()) * (c.y() - a.y())
           - (b.y() - a.y()) * (c.x() - a.x())) > 1e-12;
  }

  bool pointInTriangle(const Vector2& p, const Vector2& a,
                       const Vector2& b, const Vector2& c)
  {
    auto sign = [](const Vector2& p1, const Vector2& p2, const Vector2& p3) {
      return (p1.x() - p3.x()) * (p2.y() - p3.y())
           - (p2.x() - p3.x()) * (p1.y() - p3.y());
    };
    bool b1 = sign(p, a, b) < 0.0;
    bool b2 = sign(p, b, c) < 0.0;
    bool b3 = sign(p, c, a) < 0.0;
    return (b1 == b2) && (b2 == b3);
  }

  // Ear clipping triangulation using local indices [0..N-1]
  std::vector<TriangleLocal> triangulateEarClipping(
    const std::vector<Vector2>& pts2D,
    std::vector<int> V)
  {
    std::vector<TriangleLocal> result;
    size_t n = V.size();
    while (n > 3)
    {
      bool earFound = false;
      for (size_t i = 0; i < n; ++i)
      {
        int iPrev = V[(i + n - 1) % n];
        int iCurr = V[i];
        int iNext = V[(i + 1) % n];
        const Vector2& prev = pts2D[iPrev];
        const Vector2& curr = pts2D[iCurr];
        const Vector2& next = pts2D[iNext];
        if (!isCCW(prev, curr, next)) continue;
        bool hasInside = false;
        for (size_t j = 0; j < n; ++j)
        {
          int vi = V[j];
          if (vi == iPrev || vi == iCurr || vi == iNext) continue;
          if (pointInTriangle(pts2D[vi], prev, curr, next))
          {
            hasInside = true;
            break;
          }
        }
        if (!hasInside)
        {
          result.push_back({iPrev, iCurr, iNext});
          V.erase(V.begin() + i);
          n = V.size();
          earFound = true;
          break;
        }
      }
      if (!earFound)
        throw std::runtime_error("Ear clipping failed: polygon may be non-simple or degenerate");
    }
    if (V.size() == 3)
      result.push_back({V[0], V[1], V[2]});
    else
      throw std::runtime_error("Triangulation ended with invalid vertex count: " + std::to_string(V.size()));
    return result;
  }

  void projectToBestFitPlane(
    const std::vector<Vector3>& points3D,
    std::vector<Vector2>& projected2D,
    Vector3& origin)
  {
    origin = Vector3::Zero();
    for (const auto& p : points3D) origin += p;
    origin /= points3D.size();

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& p : points3D)
    {
      Eigen::Vector3d d = p - origin;
      cov += d * d.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
    Vector3 axis1 = eig.eigenvectors().col(1);
    Vector3 axis2 = eig.eigenvectors().col(2);
    projected2D.reserve(points3D.size());
    for (const auto& p : points3D)
    {
      Vector3 d = p - origin;
      projected2D.emplace_back(d.dot(axis1), d.dot(axis2));
    }
  }
}

int main()
{
  Mesh mesh;
  mesh.load("LeftAtriumTetra.mesh", IO::FileFormat::MEDIT);

  std::cout << "---------------------" << AnnulusLabel << "\n";
  std::cout << "Closing annuli in mesh with label " << AnnulusLabel << "...\n";
  std::cout << "---------------------" << AnnulusLabel << "\n";

  std::unordered_map<int, std::vector<int>> adj;
  for (size_t i = 0; i < mesh.getPolytopeCount(1); ++i)
  {
    auto e = mesh.getPolytope(1, i);
    if (e->getAttribute() == AnnulusLabel)
    {
      auto vs = e->getVertices();
      adj[vs[0]].push_back(vs[1]);
      adj[vs[1]].push_back(vs[0]);
    }
  }
  for (const auto& [v, nbr] : adj)
    if (nbr.size() != 2)
      throw std::runtime_error("Loop requirement not satisfied!");

  std::vector<int> loop;
  int current = adj.begin()->first;
  int prev = -1;
  do {
    loop.push_back(current);
    const auto& nbr = adj[current];
    int next = (nbr[0] == prev ? nbr[1] : nbr[0]);
    prev = current;
    current = next;
  } while (current != loop.front());

  std::cout << "Found a loop with " << loop.size() << " vertices.\n";

  std::vector<Vector3> points3D;
  for (int v : loop)
    points3D.push_back(mesh.getVertex(v)->getCoordinates());

  std::vector<Vector2> points2D;
  Vector3 origin;
  projectToBestFitPlane(points3D, points2D, origin);

  // Ensure CCW winding of projected2D
  double area = 0;
  for (size_t i = 0; i < points2D.size(); ++i)
  {
    const auto& p0 = points2D[i];
    const auto& p1 = points2D[(i + 1) % points2D.size()];
    area += (p0.x() * p1.y() - p1.x() * p0.y());
  }
  if (area < 0)
  {
    std::reverse(points2D.begin(), points2D.end());
    std::reverse(loop.begin(), loop.end());
  }

  // Prepare local indices [0..N-1]
  std::vector<int> localIdx(loop.size());
  std::iota(localIdx.begin(), localIdx.end(), 0);

  // Triangulate using ear clipping
  std::vector<TriangleLocal> localTris;
  try {
    localTris = triangulateEarClipping(points2D, localIdx);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\nFalling back to fan triangulation...\n";
    for (size_t i = 1; i + 1 < loop.size(); ++i)
      localTris.push_back({0, (int)i, (int)(i + 1)});
  }

  // Map back to global vertex indices and output
  std::cout << "Triangles\n" << localTris.size() << "\n";
  for (const auto& lt : localTris)
  {
    int g0 = loop[lt.i0];
    int g1 = loop[lt.i1];
    int g2 = loop[lt.i2];
    std::cout << g0 + 1 << " " << g1 + 1 << " " << g2 + 1 << " " << AnnulusLabel << "\n";
  }

  return 0;
}
