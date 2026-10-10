/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <Rodin/Alert.h>
#include <Rodin/Geometry/Polytope.h>
#include <Rodin/Math/SpatialVector.h>

#include <boost/json.hpp>
#include "Ball.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>

namespace Examples
{
  using namespace Rodin;

  /**
   * Offline Euclidean covering search on Rodin reference polytopes.
   *
   * Architecture: Rodin traits define the bounded reference domain. Each
   * Voronoi cell is clipped by site-relative bisectors; its vertices give the
   * continuous covering radius in floating point. Checked smallest enclosing
   * balls recenter frozen cells, then the Voronoi partition is rebuilt.
   * JSON includes reference geometry, initial/final cells and search diagnostics;
   * visualizers need not reconstruct or solve the geometry.
   */
  class Witness
  {
    public:
      using Type = Geometry::Polytope::Type;
      using Point = Math::SpatialVector<Real>;
      using Points = std::vector<Point>;
      using Faces = std::vector<Points>;
      using Clock = std::chrono::steady_clock;

      // Numerical geometry tolerances in unit-reference coordinates.
      static constexpr Real GeometryTolerance = 1e-10;
      static constexpr Real MergeTolerance = 1e-12;

      struct Options
      {
        size_t count = 6, iterations = 100, resolution = 12, seed = 13, maxEvaluations = 0;
        Real stepTolerance = 1e-6;
        std::string geometry = "triangle", output, evaluate;
      };

      struct Plane
      {
        Point normal;
        Real bound;
      };

      struct Coverage
      {
        Real radius = std::numeric_limits<Real>::infinity();
        std::vector<Points> cells;
        Points holes;
      };

      static const std::vector<std::pair<std::string, Type>>& geometries()
      {
        static const std::vector<std::pair<std::string, Type>> values = {
          {"point", Type::Point}, {"segment", Type::Segment},
          {"triangle", Type::Triangle}, {"quadrilateral", Type::Quadrilateral},
          {"tetrahedron", Type::Tetrahedron}, {"pyramid", Type::Pyramid},
          {"hexahedron", Type::Hexahedron}, {"wedge", Type::Wedge}};
        return values;
      }

      Witness(Type type, Options options)
        : m_options(std::move(options)),
          m_dimension(Geometry::Polytope::Traits(type).getDimension()),
          m_coordinates(std::max(size_t(1), m_dimension)), m_center(m_coordinates)
      {
        const Geometry::Polytope::Traits traits(type);
        for (size_t i = 0; i < traits.getVertexCount(); ++i)
        {
          m_vertices.push_back(traits.getVertex(i));
          m_center += traits.getVertex(i);
        }
        m_center /= Real(m_vertices.size());
        const auto& halfSpace = traits.getHalfSpace();
        for (Eigen::Index i = 0; i < halfSpace.matrix.rows(); ++i)
          m_planes.push_back({Point(halfSpace.matrix.row(i).transpose()), halfSpace.vector(i)});
        if (m_dimension == 2)
          m_faces.push_back(m_vertices);
        else if (m_dimension == 3)
        {
          for (const auto& plane : m_planes)
          {
            Points face;
            for (const auto& vertex : m_vertices)
            {
              if (std::abs(plane.normal.dot(vertex)-plane.bound) <= GeometryTolerance)
                face.push_back(vertex);
            }
            order(face, plane.normal);
            m_faces.push_back(std::move(face));
          }
        }
      }

      Coverage evaluate(const Points& sites, bool record = true) const
      {
        Coverage result;
        if (sites.empty())
          return result;
        for (size_t i = 0; i < sites.size(); ++i)
        {
          if (!contains(sites[i]))
            return result;
          for (size_t j = 0; j < i; ++j)
          {
            if ((sites[i]-sites[j]).norm() <= GeometryTolerance)
              return result;
          }
        }
        Real radiusSquared = 0;
        Points candidates;
        for (size_t i = 0; i < sites.size(); ++i)
        {
          Points cell;
          if (m_dimension == 0)
            cell = m_vertices;
          else if (m_dimension == 1)
          {
            Real lower = 0, upper = 1;
            for (const auto& other : sites)
            {
              if (other(0) < sites[i](0))
                lower = std::max(lower, (other(0)+sites[i](0))/2);
              else if (other(0) > sites[i](0))
                upper = std::min(upper, (other(0)+sites[i](0))/2);
            }
            cell = {Point{lower}, Point{upper}};
          }
          else
          {
            Faces faces = m_faces;
            for (auto& face : faces)
            {
              for (auto& vertex : face)
                vertex -= sites[i];
            }
            for (size_t j = 0; j < sites.size(); ++j)
            {
              if (i == j)
                continue;
              const Point difference(sites[j]-sites[i]);
              const Real length = difference.norm();
              const Plane bisector{Point(difference/length), length/2};
              Faces clipped;
              Points cap;
              for (const auto& face : faces)
              {
                Points polygon;
                for (size_t k = 0; k < face.size(); ++k)
                {
                  const Point& a = face[k];
                  const Point& b = face[(k+1)%face.size()];
                  const Real da = bisector.normal.dot(a)-bisector.bound;
                  const Real db = bisector.normal.dot(b)-bisector.bound;
                  if (da <= 0)
                    append(polygon, a);
                  if ((da <= 0 && db > 0) || (db <= 0 && da > 0))
                  {
                    const Point intersection(a+(da/(da-db))*(b-a));
                    append(polygon, intersection);
                    append(cap, intersection);
                  }
                }
                if (polygon.size() >= 3)
                  clipped.push_back(std::move(polygon));
              }
              if (m_dimension == 3 && cap.size() >= 3)
              {
                order(cap, bisector.normal);
                clipped.push_back(std::move(cap));
              }
              faces = std::move(clipped);
              if (faces.empty())
                return result;
            }
            for (const auto& face : faces)
            {
              for (const auto& vertex : face)
                append(cell, Point(vertex+sites[i]));
            }
          }
          for (const auto& vertex : cell)
          {
            radiusSquared = std::max(radiusSquared, nearestSquared(vertex, sites));
            if (record)
              append(candidates, vertex);
          }
          if (record)
            result.cells.push_back(std::move(cell));
        }
        result.radius = std::sqrt(radiusSquared);
        if (record)
        {
          for (const auto& vertex : candidates)
          {
            if (std::abs(std::sqrt(nearestSquared(vertex, sites))-result.radius) <= GeometryTolerance)
              append(result.holes, vertex);
          }
        }
        return result;
      }

      boost::json::object run(const Points& supplied = {})
      {
        const auto start = Clock::now();
        if (!m_options.count || (m_dimension == 0 && m_options.count != 1))
          Alert::Exception() << "Invalid witness count for this geometry." << Alert::Raise;
        const size_t resolution = m_options.resolution;
        size_t initialResolution = 0, initialLatticeSize = 0;
        std::string initialization = "supplied";
        Points initial = supplied;
        if (initial.empty())
        {
          if (m_options.count == 1)
          {
            initial.push_back(m_center);
            initialization = "vertex_barycenter";
          }
          else
          {
            Points candidates;
            do
            {
              candidates = grid(++initialResolution);
            } while (candidates.size() < m_options.count);
            initialLatticeSize = candidates.size();
            if (candidates.size() == m_options.count)
            {
              initial = candidates;
              initialization = "complete_uniform_lattice";
            }
            else
            {
              initialization = "uniform_lattice_farthest_subset";
              const auto closest = std::min_element(candidates.begin(), candidates.end(), [&](const Point& a, const Point& b) {
                return (a-m_center).squaredNorm() < (b-m_center).squaredNorm();
              });
              initial.push_back(*closest);
              while (initial.size() < m_options.count)
              {
                const auto farthest = std::max_element(candidates.begin(), candidates.end(), [&](const Point& a, const Point& b) {
                  return nearestSquared(a, initial) < nearestSquared(b, initial);
                });
                initial.push_back(*farthest);
              }
            }
          }
        }
        const Coverage initialCoverage = evaluate(initial);
        if (!std::isfinite(initialCoverage.radius))
          Alert::Exception() << "Input witnesses must be distinct, finite and inside the reference geometry." << Alert::Raise;
        Points best = initial;
        Real bestRadius = initialCoverage.radius;
        boost::json::array history;
        size_t sweeps = 0, fallbacks = 0;
        Real localSeconds = 0, coverageSeconds = 0, motion = 0;
        std::string termination = supplied.empty() ? "iteration_limit" : "evaluation";
        history.push_back(boost::json::object{{"sweep", 0}, {"radius", bestRadius}});
        if (supplied.empty() && m_dimension == 1)
        {
          best.clear();
          for (size_t i = 0; i < m_options.count; ++i)
            best.push_back(Point{(Real(i)+0.5)/Real(m_options.count)});
          termination = "analytic";
          history.clear();
          history.push_back(boost::json::object{{"sweep", 0}, {"radius", Real(0.5)/m_options.count}});
        }
        else if (supplied.empty() && m_dimension > 1 && m_options.iterations)
        {
          Coverage coverage = initialCoverage;
          for (; sweeps < m_options.iterations && !exhausted(); ++sweeps)
          {
            const auto localStart = Clock::now();
            Points candidate;
            motion = 0;
            for (size_t i = 0; i < best.size(); ++i)
            {
              const Ball ball(coverage.cells[i], m_options.seed);
              fallbacks += ball.usedFallback();
              motion = std::max(motion, (ball.getCenter()-best[i]).norm());
              candidate.push_back(ball.getCenter());
            }
            localSeconds += std::chrono::duration<Real>(Clock::now()-localStart).count();
            const auto coverageStart = Clock::now();
            Coverage next = evaluate(candidate);
            ++m_evaluations;
            coverageSeconds += std::chrono::duration<Real>(Clock::now()-coverageStart).count();
            if (!std::isfinite(next.radius) || next.radius > coverage.radius+GeometryTolerance)
            {
              termination = "coverage_rejected";
              ++sweeps;
              break;
            }
            best = std::move(candidate);
            coverage = std::move(next);
            history.push_back(boost::json::object{{"sweep", sweeps+1}, {"radius", coverage.radius}, {"motion", motion}});
            if (motion <= m_options.stepTolerance)
            {
              termination = "motion_tolerance";
              ++sweeps;
              break;
            }
          }
          if (exhausted() && termination == "iteration_limit")
            termination = "evaluation_limit";
        }
        const Coverage finalCoverage = evaluate(best);
        if (!std::isfinite(finalCoverage.radius))
          Alert::Exception() << "Final Voronoi geometry is numerically unresolved." << Alert::Raise;
        const Real lowerBound = sampledRadius(best, grid(2*resolution));
        if (lowerBound > finalCoverage.radius+GeometryTolerance)
          Alert::Exception() << "Independent lattice check exceeds the Voronoi radius." << Alert::Raise;
        const Real seconds = std::chrono::duration<Real>(Clock::now()-start).count();
        boost::json::array edges;
        for (size_t i = 0; i < m_vertices.size(); ++i)
        {
          for (size_t j = i+1; j < m_vertices.size(); ++j)
          {
            Points shared;
            for (const auto& plane : m_planes)
            {
              if (std::abs(plane.normal.dot(m_vertices[i])-plane.bound) <= GeometryTolerance &&
                  std::abs(plane.normal.dot(m_vertices[j])-plane.bound) <= GeometryTolerance)
                shared.push_back(plane.normal);
            }
            if (m_dimension == 1 || (m_dimension == 2 && !shared.empty()) ||
                (m_dimension == 3 && shared.size() >= 2 &&
                 shared[0].getData().cross(shared[1].getData()).norm() > GeometryTolerance))
              edges.push_back(boost::json::array{i, j});
          }
        }
        return boost::json::object{
          {"schema_version", 2}, {"geometry", m_options.geometry}, {"dimension", m_dimension},
          {"count", best.size()}, {"points", json(best)}, {"covering_radius", finalCoverage.radius},
          {"worst_locations", json(finalCoverage.holes)}, {"cells", json(finalCoverage.cells)},
          {"initial_points", json(initial)}, {"initial_radius", initialCoverage.radius},
          {"initial_cells", json(initialCoverage.cells)},
          {"initial_worst_locations", json(initialCoverage.holes)},
          {"initialization", boost::json::object{
            {"method", initialization}, {"resolution", initialResolution},
            {"lattice_size", initialLatticeSize}}},
          {"reference_vertices", json(m_vertices)}, {"reference_edges", std::move(edges)},
          {"validation_lattice_lower_bound", lowerBound},
          {"optimality", m_dimension == 0 ? "unique_admissible_set" :
            (m_dimension == 1 && supplied.empty() ? "analytic_interval_solution" : "not_certified")},
          {"radius_evaluation", "clipped_voronoi_float64"},
          {"history", std::move(history)},
          {"search", boost::json::object{
            {"implementation", "checked_enclosing_ball"},
            {"enclosing_ball_backend", "CGAL::Min_sphere_of_spheres_d"},
            {"mode", supplied.empty() ? "search" : "evaluation"},
            {"objective", "exact"}, {"resolution", resolution},
            {"iterations", m_options.iterations}, {"seed", m_options.seed},
            {"step_tolerance", m_options.stepTolerance}, {"sweeps", sweeps},
            {"termination", termination}, {"motion", motion}, {"support_fallbacks", fallbacks},
            {"local_seconds", localSeconds}, {"coverage_seconds", coverageSeconds},
            {"max_evaluations", m_options.maxEvaluations}, {"evaluations", m_evaluations},
            {"wall_seconds", seconds}, {"threads", 1}}}};
      }

      Points readPoints(const boost::json::array& values) const
      {
        Points points;
        for (const auto& value : values)
        {
          const auto& coordinates = value.as_array();
          if (coordinates.size() != m_coordinates)
            Alert::Exception() << "Incorrect witness coordinate dimension." << Alert::Raise;
          Point point(m_coordinates);
          for (size_t axis = 0; axis < m_coordinates; ++axis)
            point(axis) = boost::json::value_to<Real>(coordinates[axis]);
          points.push_back(std::move(point));
        }
        return points;
      }

    private:
      bool contains(const Point& point) const
      {
        for (size_t axis = 0; axis < m_coordinates; ++axis)
        {
          if (!std::isfinite(point(axis)))
            return false;
        }
        if (!m_dimension)
          return point.norm() <= GeometryTolerance;
        for (const auto& plane : m_planes)
        {
          if (plane.normal.dot(point) > plane.bound+GeometryTolerance)
            return false;
        }
        return true;
      }

      void append(Points& points, const Point& candidate) const
      {
        for (const auto& point : points)
        {
          if ((point-candidate).norm() <= MergeTolerance)
            return;
        }
        points.push_back(candidate);
      }

      void order(Points& face, const Point& normal) const
      {
        Point center(m_coordinates);
        for (const auto& vertex : face)
          center += vertex;
        center /= Real(face.size());
        const Point u(Point(face.front()-center).normalized());
        const Point v(normal.getData().cross(u.getData()));
        std::sort(face.begin(), face.end(), [&](const Point& a, const Point& b) {
          const Point da(a-center), db(b-center);
          return std::atan2(da.dot(v), da.dot(u)) < std::atan2(db.dot(v), db.dot(u));
        });
      }

      Real nearestSquared(const Point& point, const Points& sites) const
      {
        Real distance = std::numeric_limits<Real>::infinity();
        for (const auto& site : sites)
          distance = std::min(distance, Real((point-site).squaredNorm()));
        return distance;
      }

      Real sampledRadius(const Points& sites, const Points& lattice) const
      {
        Real squared = 0;
        for (const auto& point : lattice)
          squared = std::max(squared, nearestSquared(point, sites));
        return std::sqrt(squared);
      }

      Points grid(size_t resolution) const
      {
        if (!m_dimension)
          return m_vertices;
        Points points;
        Point point(m_coordinates);
        auto visit = [&](auto&& self, size_t axis) -> void {
          if (axis == m_dimension)
          {
            if (contains(point))
              points.push_back(point);
            return;
          }
          for (size_t i = 0; i <= resolution; ++i)
          {
            point(axis) = Real(i)/Real(resolution);
            self(self, axis+1);
          }
        };
        visit(visit, 0);
        return points;
      }

      bool exhausted() const
      {
        return m_options.maxEvaluations && m_evaluations >= m_options.maxEvaluations;
      }

      boost::json::array json(const Points& points) const
      {
        boost::json::array array;
        for (const auto& point : points)
        {
          boost::json::array coordinates;
          for (size_t axis = 0; axis < m_coordinates; ++axis)
            coordinates.push_back(point(axis));
          array.push_back(std::move(coordinates));
        }
        return array;
      }

      boost::json::array json(const std::vector<Points>& cells) const
      {
        boost::json::array array;
        for (const auto& cell : cells)
          array.push_back(json(cell));
        return array;
      }

      const Options m_options;
      const size_t m_dimension, m_coordinates;
      Point m_center;
      Points m_vertices;
      Faces m_faces;
      std::vector<Plane> m_planes;
      size_t m_evaluations = 0;
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using Examples::Witness;
  try
  {
    Witness::Options options;
    bool countSet = false;
    for (int i = 1; i < argc; ++i)
    {
      std::string argument = argv[i], name, value;
      if (argument == "--help")
      {
        std::cout << "SWIFT_Witness [n] [--geometry triangle|...|all] [--output result.json]\n"
          "  --iterations 100 --resolution 12 --seed 13\n"
          "  --step-tolerance 1e-6 --max-evaluations 0 --evaluate input.json\n"
          "All witness positions are free. Search optimality is not certified.\n";
        return 0;
      }
      if (argument.starts_with("--"))
      {
        const auto equal = argument.find('=');
        name = argument.substr(2, equal == std::string::npos ? equal : equal-2);
        if (equal != std::string::npos)
          value = argument.substr(equal+1);
        else if (i+1 < argc)
          value = argv[++i];
        else
          Alert::Exception() << "Missing value for --" << name << Alert::Raise;
      }
      else
      {
        name = "n";
        value = argument;
      }
      if (name == "geometry")
        options.geometry = value;
      else if (name == "output")
        options.output = value;
      else if (name == "evaluate")
        options.evaluate = value;
      else if (name == "step-tolerance")
      {
        size_t consumed = 0;
        options.stepTolerance = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(options.stepTolerance) || options.stepTolerance <= 0)
          Alert::Exception() << "Invalid step tolerance." << Alert::Raise;
      }
      else
      {
        size_t number;
        const auto parsed = std::from_chars(value.data(), value.data()+value.size(), number);
        if (parsed.ec != std::errc() || parsed.ptr != value.data()+value.size())
          Alert::Exception() << "Invalid integer for --" << name << Alert::Raise;
        if (name == "n" && !countSet)
        {
          options.count = number;
          countSet = true;
        }
        else if (name == "iterations") options.iterations = number;
        else if (name == "resolution") options.resolution = number;
        else if (name == "seed") options.seed = number;
        else if (name == "max-evaluations") options.maxEvaluations = number;
        else Alert::Exception() << "Unknown or repeated option: " << name << Alert::Raise;
      }
    }
    if (!options.count || !options.resolution)
      Alert::Exception() << "Count and resolution must be positive." << Alert::Raise;
    boost::json::value input;
    if (!options.evaluate.empty())
    {
      std::ifstream stream(options.evaluate);
      if (!stream)
        Alert::Exception() << "Cannot read " << options.evaluate << Alert::Raise;
      input = boost::json::parse(std::string(std::istreambuf_iterator<char>(stream), {}));
      options.geometry = std::string(input.at("geometry").as_string());
      options.count = input.at("points").as_array().size();
    }
    boost::json::value output;
    boost::json::array all;
    bool found = false;
    for (const auto& [name, type] : Witness::geometries())
    {
      if (options.geometry != "all" && options.geometry != name)
        continue;
      auto settings = options;
      settings.geometry = name;
      if (name == "point" && options.geometry == "all")
        settings.count = 1;
      Witness generator(type, settings);
      const auto result = input.is_null() ? generator.run()
        : generator.run(generator.readPoints(input.at("points").as_array()));
      std::cout << name << ": n=" << result.at("count") << ", radius="
        << result.at("covering_radius") << ", seconds="
        << result.at("search").at("wall_seconds") << std::endl;
      if (options.geometry == "all")
        all.push_back(result);
      else
        output = result;
      found = true;
    }
    if (!found)
      Alert::Exception() << "Unknown reference geometry: " << options.geometry << Alert::Raise;
    if (options.geometry == "all")
      output = std::move(all);
    const std::filesystem::path path = options.output.empty()
      ? "witness-"+options.geometry+"-"+std::to_string(options.count)+".json" : options.output;
    if (path.has_parent_path())
      std::filesystem::create_directories(path.parent_path());
    std::ofstream stream(path);
    stream << boost::json::serialize(output) << '\n';
    if (!stream)
      Alert::Exception() << "Cannot write " << path.string() << Alert::Raise;
    std::cout << path << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
