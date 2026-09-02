/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  struct Options
  {
      std::string input;
      std::string output;
      Real thickness = 0;
      size_t layers = 1;
      bool centered = false;
      bool zDirection = false;
      Optional<Attribute> wedgeAttribute;
      Attribute bottomAttribute = 1;
      Attribute topAttribute = 2;
      Attribute sideAttribute = 3;
      IO::FileFormat inputFormat = IO::FileFormat::MEDIT;
      IO::FileFormat outputFormat = IO::FileFormat::MEDIT;
  };

  struct TriangleData
  {
      std::array<Index, 3> vertices;
      Attribute attribute;
  };

  struct EdgeData
  {
      Index first;
      Index second;
      size_t count;
  };

  void usage(const char* exe)
  {
    std::cerr
      << "Usage: " << exe << " [options] <input.mesh> <output.mesh> <thickness>\n"
      << "\n"
      << "Options:\n"
      << "  --layers <n>             Wedge layers through thickness (default: 1)\n"
      << "  --centered               Extrude +/- thickness/2 around the input surface\n"
      << "  --direction <normal|z>   Extrude along vertex normals or +z (default: "
         "normal)\n"
      << "  --input-format <medit|mfem>\n"
      << "  --output-format <medit|mfem>\n"
      << "  --wedge-attribute <n>    Override wedge attributes\n"
      << "  --bottom-attribute <n>   Bottom surface attribute (default: 1)\n"
      << "  --top-attribute <n>      Top surface attribute (default: 2)\n"
      << "  --side-attribute <n>     Lateral boundary attribute (default: 3)\n";
  }

  Real parseReal(const std::string& value, const std::string& name)
  {
    char* end = nullptr;
    const Real out = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0')
      Alert::Exception() << "Invalid value for " << name << ": " << value << Alert::Raise;
    return out;
  }

  size_t parseSize(const std::string& value, const std::string& name)
  {
    char* end = nullptr;
    const unsigned long out = std::strtoul(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0')
      Alert::Exception() << "Invalid value for " << name << ": " << value << Alert::Raise;
    return static_cast<size_t>(out);
  }

  Attribute parseAttribute(const std::string& value, const std::string& name)
  {
    return static_cast<Attribute>(parseSize(value, name));
  }

  IO::FileFormat parseFormat(const std::string& value)
  {
    if (value == "medit")
      return IO::FileFormat::MEDIT;
    if (value == "mfem")
      return IO::FileFormat::MFEM;
    Alert::Exception() << "Unsupported mesh format: " << value
                       << ". Expected medit or mfem." << Alert::Raise;
    return IO::FileFormat::MEDIT;
  }

  std::string takeValue(int& i, int argc, char** argv, const std::string& name)
  {
    const std::string arg(argv[i]);
    const std::string prefix = name + "=";
    if (arg.rfind(prefix, 0) == 0)
      return arg.substr(prefix.size());
    if (i + 1 >= argc)
      Alert::Exception() << "Missing value for " << name << Alert::Raise;
    return argv[++i];
  }

  Options parseOptions(int argc, char** argv)
  {
    Options options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg == "--help" || arg == "-h")
      {
        usage(argv[0]);
        std::exit(0);
      }
      else if (arg == "--centered")
      {
        options.centered = true;
      }
      else if (arg == "--layers" || arg.rfind("--layers=", 0) == 0)
      {
        options.layers = parseSize(takeValue(i, argc, argv, "--layers"), "--layers");
      }
      else if (arg == "--direction" || arg.rfind("--direction=", 0) == 0)
      {
        const std::string value = takeValue(i, argc, argv, "--direction");
        if (value == "normal")
          options.zDirection = false;
        else if (value == "z")
          options.zDirection = true;
        else
          Alert::Exception() << "Invalid --direction: " << value
                             << ". Expected normal or z." << Alert::Raise;
      }
      else if (arg == "--input-format" || arg.rfind("--input-format=", 0) == 0)
      {
        options.inputFormat = parseFormat(takeValue(i, argc, argv, "--input-format"));
      }
      else if (arg == "--output-format" || arg.rfind("--output-format=", 0) == 0)
      {
        options.outputFormat = parseFormat(takeValue(i, argc, argv, "--output-format"));
      }
      else if (arg == "--wedge-attribute" || arg.rfind("--wedge-attribute=", 0) == 0)
      {
        options.wedgeAttribute = parseAttribute(
          takeValue(i, argc, argv, "--wedge-attribute"), "--wedge-attribute");
      }
      else if (arg == "--bottom-attribute" || arg.rfind("--bottom-attribute=", 0) == 0)
      {
        options.bottomAttribute = parseAttribute(
          takeValue(i, argc, argv, "--bottom-attribute"), "--bottom-attribute");
      }
      else if (arg == "--top-attribute" || arg.rfind("--top-attribute=", 0) == 0)
      {
        options.topAttribute =
          parseAttribute(takeValue(i, argc, argv, "--top-attribute"), "--top-attribute");
      }
      else if (arg == "--side-attribute" || arg.rfind("--side-attribute=", 0) == 0)
      {
        options.sideAttribute = parseAttribute(
          takeValue(i, argc, argv, "--side-attribute"), "--side-attribute");
      }
      else if (!arg.empty() && arg[0] == '-')
      {
        Alert::Exception() << "Unknown option: " << arg << Alert::Raise;
      }
      else
      {
        positional.push_back(arg);
      }
    }

    if (positional.size() != 3)
    {
      usage(argv[0]);
      Alert::Exception() << "Expected input mesh, output mesh, and thickness."
                         << Alert::Raise;
    }

    options.input = positional[0];
    options.output = positional[1];
    options.thickness = parseReal(positional[2], "thickness");

    if (options.thickness <= 0)
      Alert::Exception() << "Thickness must be positive." << Alert::Raise;
    if (options.layers == 0)
      Alert::Exception() << "--layers must be at least 1." << Alert::Raise;
    if (options.input == options.output)
      Alert::Exception() << "Input and output paths must be different." << Alert::Raise;
    return options;
  }

  Math::SpatialPoint point3(const LocalMesh& mesh, Index vertex)
  {
    Math::SpatialPoint out(3);
    out(0) = 0;
    out(1) = 0;
    out(2) = 0;

    const auto x = mesh.getVertexCoordinates(vertex);
    for (size_t i = 0; i < std::min<size_t>(3, x.size()); ++i)
      out(i) = x(i);
    return out;
  }

  std::vector<TriangleData> collectTriangles(const LocalMesh& mesh)
  {
    std::vector<TriangleData> triangles;
    triangles.reserve(mesh.getCellCount());

    for (auto it = mesh.getCell(); it; ++it)
    {
      if (it->getGeometry() != Polytope::Type::Triangle)
      {
        Alert::Exception() << "Expected a triangular surface mesh; encountered "
                           << it->getGeometry() << " cell." << Alert::Raise;
      }

      const auto& vs = it->getVertices();
      triangles.push_back(
        {{vs(0), vs(1), vs(2)}, it->getAttribute().value_or(Attribute{0})});
    }

    if (triangles.empty())
      Alert::Exception() << "Input mesh has no triangular cells." << Alert::Raise;
    return triangles;
  }

  std::vector<Math::SpatialPoint> computeVertexNormals(
    const LocalMesh& mesh, const std::vector<TriangleData>& triangles)
  {
    std::vector<Math::SpatialPoint> normals(
      mesh.getVertexCount(), Math::SpatialPoint::Zero(3));

    for (const auto& triangle : triangles)
    {
      const auto a = point3(mesh, triangle.vertices[0]);
      const auto b = point3(mesh, triangle.vertices[1]);
      const auto c = point3(mesh, triangle.vertices[2]);
      const auto normal = (b - a).cross(c - a);
      if (normal.norm() == 0)
      {
        Alert::Exception() << "Encountered a degenerate triangle while computing normals."
                           << Alert::Raise;
      }

      for (Index v : triangle.vertices)
        normals[v] += normal;
    }

    for (Index v = 0; v < static_cast<Index>(normals.size()); ++v)
    {
      if (normals[v].norm() == 0)
        Alert::Exception() << "Vertex " << v << " has no usable incident normal."
                           << Alert::Raise;
      normals[v] /= normals[v].norm();
    }
    return normals;
  }

  Index layeredVertex(Index vertex, size_t layer, size_t vertexCount)
  {
    return static_cast<Index>(layer * vertexCount + vertex);
  }

  std::vector<EdgeData> collectBoundaryEdges(const std::vector<TriangleData>& triangles)
  {
    std::map<std::pair<Index, Index>, EdgeData> edges;
    const auto addEdge = [&edges](Index a, Index b) {
      const auto key = std::minmax(a, b);
      auto it = edges.find(key);
      if (it == edges.end())
        edges.emplace(key, EdgeData{a, b, 1});
      else
        it->second.count++;
    };

    for (const auto& triangle : triangles)
    {
      addEdge(triangle.vertices[0], triangle.vertices[1]);
      addEdge(triangle.vertices[1], triangle.vertices[2]);
      addEdge(triangle.vertices[2], triangle.vertices[0]);
    }

    std::vector<EdgeData> boundary;
    for (const auto& [_, edge] : edges)
    {
      if (edge.count == 1)
        boundary.push_back(edge);
    }
    return boundary;
  }

  LocalMesh extrude(const LocalMesh& input, const Options& options)
  {
    if (input.getDimension() != 2)
    {
      Alert::Exception()
        << "Expected a 2D triangular surface mesh; input topological dimension is "
        << input.getDimension() << "." << Alert::Raise;
    }

    const auto triangles = collectTriangles(input);
    const auto boundaryEdges = collectBoundaryEdges(triangles);
    const std::vector<Math::SpatialPoint> directions = options.zDirection
      ? std::vector<Math::SpatialPoint>(
          input.getVertexCount(), Math::SpatialPoint({0, 0, 1}))
      : computeVertexNormals(input, triangles);

    LocalMesh::Builder builder;
    builder.initialize(3)
      .nodes((options.layers + 1) * input.getVertexCount())
      .reserve(0, (options.layers + 1) * input.getVertexCount())
      .reserve(2, 2 * triangles.size() + options.layers * boundaryEdges.size())
      .reserve(3, options.layers * triangles.size());

    const Real offset = options.centered ? -Real(0.5) * options.thickness : Real(0);
    for (size_t layer = 0; layer <= options.layers; ++layer)
    {
      const Real t = offset +
        options.thickness * static_cast<Real>(layer) / static_cast<Real>(options.layers);
      for (Index v = 0; v < static_cast<Index>(input.getVertexCount()); ++v)
      {
        const auto x = point3(input, v) + t * directions[v];
        builder.vertex(x);
        if (const auto attr = input.getVertex(v)->getAttribute())
          builder.attribute({0, layeredVertex(v, layer, input.getVertexCount())}, *attr);
      }
    }

    for (size_t layer = 0; layer < options.layers; ++layer)
    {
      for (const auto& triangle : triangles)
      {
        Index index;
        builder.polytope(Polytope::Type::Wedge,
          {layeredVertex(triangle.vertices[0], layer, input.getVertexCount()),
            layeredVertex(triangle.vertices[1], layer, input.getVertexCount()),
            layeredVertex(triangle.vertices[2], layer, input.getVertexCount()),
            layeredVertex(triangle.vertices[0], layer + 1, input.getVertexCount()),
            layeredVertex(triangle.vertices[1], layer + 1, input.getVertexCount()),
            layeredVertex(triangle.vertices[2], layer + 1, input.getVertexCount())},
          index);
        builder.attribute(
          {3, index}, options.wedgeAttribute.value_or(triangle.attribute));
      }
    }

    for (const auto& triangle : triangles)
    {
      Index bottom;
      builder.polytope(Polytope::Type::Triangle,
        {layeredVertex(triangle.vertices[2], 0, input.getVertexCount()),
          layeredVertex(triangle.vertices[1], 0, input.getVertexCount()),
          layeredVertex(triangle.vertices[0], 0, input.getVertexCount())},
        bottom);
      builder.attribute({2, bottom}, options.bottomAttribute);

      Index top;
      builder.polytope(Polytope::Type::Triangle,
        {layeredVertex(triangle.vertices[0], options.layers, input.getVertexCount()),
          layeredVertex(triangle.vertices[1], options.layers, input.getVertexCount()),
          layeredVertex(triangle.vertices[2], options.layers, input.getVertexCount())},
        top);
      builder.attribute({2, top}, options.topAttribute);
    }

    for (const auto& edge : boundaryEdges)
    {
      for (size_t layer = 0; layer < options.layers; ++layer)
      {
        Index side;
        builder.polytope(Polytope::Type::Quadrilateral,
          {layeredVertex(edge.first, layer, input.getVertexCount()),
            layeredVertex(edge.second, layer, input.getVertexCount()),
            layeredVertex(edge.second, layer + 1, input.getVertexCount()),
            layeredVertex(edge.first, layer + 1, input.getVertexCount())},
          side);
        builder.attribute({2, side}, options.sideAttribute);
      }
    }

    return builder.finalize();
  }
}

int main(int argc, char** argv)
{
  try
  {
    const Options options = parseOptions(argc, argv);

    LocalMesh input;
    input.load(options.input, options.inputFormat);

    LocalMesh output = extrude(input, options);
    output.save(options.output, options.outputFormat);

    std::cout << "Wrote " << options.output << ": " << output.getVertexCount()
              << " vertices, " << output.getPolytopeCount(Polytope::Type::Wedge)
              << " wedges, " << output.getPolytopeCount(Polytope::Type::Triangle)
              << " triangles, " << output.getPolytopeCount(Polytope::Type::Quadrilateral)
              << " quadrilaterals." << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
