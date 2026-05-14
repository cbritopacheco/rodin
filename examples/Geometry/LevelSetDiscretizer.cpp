/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetInterfaceGraph.h>
#include <Rodin/IO.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static LocalMesh makeInterfaceMesh(const InterfaceGraph& graph)
{
  LocalMesh::Builder builder;
  builder
    .initialize(2)
    .nodes(graph.vertices.size())
    .reserve(1, graph.edges.size());

  for (const auto& vertex : graph.vertices)
    builder.vertex(vertex.x);

  for (Index i = 0; i < graph.edges.size(); ++i)
  {
    const auto& edge = graph.edges[i];
    builder
      .polytope(Polytope::Type::Segment, {edge.v0, edge.v1})
      .attribute({1, i}, edge.interfaceAttribute);
  }

  return builder.finalize();
}

int main(int, char**)
{
  static constexpr Attribute Gamma = 10;
  static constexpr size_t n = 20;

  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
  mesh.scale(1.0 / static_cast<Real>(n - 1));
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  mesh.save("LevelSetInterfaceGraph.background.mesh", IO::FileFormat::MEDIT);

  P1 phiSpace(mesh);
  GridFunction phi(phiSpace);

  struct Circle
  {
    Real cx;
    Real cy;
    Real r;
  };

  const std::array<Circle, 3> circles = {{
    {0.30, 0.35, 0.16},
    {0.70, 0.35, 0.14},
    {0.50, 0.72, 0.12}
  }};
  const Real openLineX = 0.08;

  for (Index i = 0; i < mesh.getVertexCount(); ++i)
  {
    const auto x = mesh.getVertexCoordinates(i);
    Real value = x[0] - openLineX;
    for (const auto& circle : circles)
    {
      const Real dx = x[0] - circle.cx;
      const Real dy = x[1] - circle.cy;
      value *= dx * dx + dy * dy - circle.r * circle.r;
    }
    phi[i] = value;
  }

  const auto graph = LevelSetInterfaceGraph(phi)
    .setSignTolerance(1e-12)
    .setInterfaceAttribute(Gamma)
    .extract();

  const auto interfaceMesh = makeInterfaceMesh(graph);
  interfaceMesh.save("LevelSetInterfaceGraph.interface.mesh", IO::FileFormat::MEDIT);

  IO::XDMF xdmf("LevelSetInterfaceGraph");
  xdmf.grid("background").setMesh(mesh);
  xdmf.grid("interface").setMesh(interfaceMesh);
  xdmf.write();
  xdmf.close();

  std::ofstream points("LevelSetInterfaceGraph.points.csv");
  points << "id,x,y,kind,original_vertex,parent_edge,t\n";
  for (Index i = 0; i < graph.vertices.size(); ++i)
  {
    const auto& v = graph.vertices[i];
    points << i << ','
           << v.x[0] << ','
           << v.x[1] << ','
           << (v.kind == InterfaceVertexKind::OriginalVertex ? "original" : "edge_cut")
           << ',';
    if (v.originalVertex)
      points << *v.originalVertex;
    points << ',';
    if (v.parentEdge)
      points << *v.parentEdge;
    points << ',';
    if (v.t)
      points << *v.t;
    points << '\n';
  }

  std::ofstream segments("LevelSetInterfaceGraph.segments.csv");
  segments << "id,v0,v1,attribute,provenance_count\n";
  for (Index i = 0; i < graph.edges.size(); ++i)
  {
    const auto& e = graph.edges[i];
    segments << i << ','
             << e.v0 << ','
             << e.v1 << ',';
    if (e.interfaceAttribute)
      segments << *e.interfaceAttribute;
    segments << ','
             << e.provenance.size() << '\n';
  }

  std::ofstream provenance("LevelSetInterfaceGraph.edge-provenance.csv");
  provenance << "edge,entry,parent_cell,parent_edge0,parent_edge1,parent_cell_attribute\n";
  for (Index i = 0; i < graph.edges.size(); ++i)
  {
    const auto& edge = graph.edges[i];
    for (Index j = 0; j < edge.provenance.size(); ++j)
    {
      const auto& p = edge.provenance[j];
      provenance << i << ','
                 << j << ','
                 << p.parentCell << ','
                 << p.parentEdges[0] << ','
                 << p.parentEdges[1] << ',';
      if (p.parentCellAttribute)
        provenance << *p.parentCellAttribute;
      provenance << '\n';
    }
  }

  std::ofstream chains("LevelSetInterfaceGraph.chains.csv");
  chains << "chain,closed,position,vertex,edge\n";
  for (Index i = 0; i < graph.chains.size(); ++i)
  {
    const auto& chain = graph.chains[i];
    for (Index j = 0; j < chain.vertices.size(); ++j)
    {
      chains << i << ','
             << (chain.closed ? 1 : 0) << ','
             << j << ','
             << chain.vertices[j] << ',';
      if (j < chain.edges.size())
        chains << chain.edges[j];
      chains << '\n';
    }
  }

  std::cout << "Wrote " << graph.vertices.size()
            << " interface vertices and " << graph.edges.size()
            << " interface segments." << std::endl;
  std::cout << "Wrote MEDIT meshes:"
            << " LevelSetInterfaceGraph.background.mesh and"
            << " LevelSetInterfaceGraph.interface.mesh." << std::endl;
  std::cout << "Wrote XDMF overlay: LevelSetInterfaceGraph.xdmf." << std::endl;
  std::cout << "Detected " << graph.chains.size()
            << " ordered interface chains." << std::endl;
  for (Index i = 0; i < graph.chains.size(); ++i)
  {
    const auto& chain = graph.chains[i];
    std::cout << "  chain " << i << " ("
              << (chain.closed ? "closed" : "open") << "): "
              << chain.vertices.size() << " vertices, "
              << chain.edges.size() << " segments" << std::endl;
  }
  return 0;
}
