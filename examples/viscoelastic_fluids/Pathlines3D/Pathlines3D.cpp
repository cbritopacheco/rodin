/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Pathlines3D.cpp
 * @brief Lagrangian pathlines through the velocity output of ArterialLesion3D.
 *
 * Reads the XDMF/HDF5 output of a 3D lesion case (every rank piece), rebuilds
 * the tetrahedral mesh, loads the velocity snapshots of the last cycle, finds
 * the inlet, seeds particles over it and integrates their trajectories in the
 * time-periodic velocity field with Rodin's characteristic tracer
 * (Variational::Flow, cell-by-cell, RK4 inside each cell). Writes the
 * pathlines as an XDMF polyline mesh with nodal fields, one row per
 * trajectory in a CSV, and the residence-time statistics per release phase.
 *
 *     Pathlines3D -pl_case <dir with case.xdmf> [options]
 *
 * Options (defaults):
 *   -pl_case      <dir>        the case folder (required)
 *   -pl_out       <stem>       output stem (default <dir>/pathlines)
 *   -pl_seeds     200          seeds over the inlet (Fibonacci disc)
 *   -pl_releases  8            release phases, equally spaced in the cycle
 *   -pl_substeps  4            tracing sub-steps per snapshot interval
 *   -pl_transits  3            stop after this many mean transit times L/U
 *   -pl_horizon   6            ... but never before this many cycles (washout)
 *   -pl_period    <s>          cycle period (default: read from run.log)
 *   -pl_cycles    1            cycles of snapshots to use (the last ones)
 *
 * What is found by the program itself:
 *   - the period, from the "[inflow] ... T=<s>" line of run.log (or -pl_period);
 *   - the axis of the vessel, as the principal direction of the vertex cloud;
 *   - the two caps (inlet, outlet): the boundary faces whose normal is
 *     aligned with the axis, split by their position along it; the inlet is
 *     the cap with negative mean flux over the cycle;
 *   - the integration length: a trajectory ends when it crosses the outlet
 *     cap ("outlet"), when it reaches another boundary ("wall", the no-slip
 *     wall is never crossed by the exact flow, so this flags a tracing
 *     limit), or after -pl_transits mean transit times L/U of the cycle
 *     ("trapped": a particle held in a recirculation), with L the axial
 *     extent and U the cycle-mean inlet speed, but no earlier than
 *     -pl_horizon cycles;
 *   - the lesion: the wall vertices whose distance to the axis departs from
 *     the parent radius R (the median over the wall) by more than 5 %; the
 *     lesion slab is their axial range. For a dilation (aneurysm) the sac is
 *     the part of the slab outside the parent cylinder, r > R.
 *
 * Residence in the lesion (the quantity that tells the fluids apart; the
 * time to the outlet is dominated by the length of the pipe): for each
 * trajectory, the time spent in the slab and in the sac, and the time from
 * the first entry to the last exit. Over the trajectories that entered, the
 * washout curve (the fraction still to leave for good, against the time
 * since entry, in cycles) and its median, p90 and half-life. Trajectories
 * stopped while still able to return (trapped, horizon) are censored: their
 * times are lower bounds, and the censored fraction is reported.
 *
 * Time: the snapshots of the last cycle are a periodic field, u(t + T) = u(t).
 * Each sub-step of length dtau freezes the field at the sub-step midpoint
 * (linear interpolation between the two bracketing snapshots), which is
 * second order in dtau; inside the sub-step the tracer integrates the frozen
 * field exactly up to its own tolerances.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <hdf5.h>

#include <Rodin/Alert.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/HDF5.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Location.h>
#include <Rodin/Math.h>
#include <Rodin/Math/RungeKutta/RK4.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/Flow.h>

namespace Rodin::Examples::ViscoelasticFluids
{
  using Rodin::Index;
  using Rodin::Real;
  using Rodin::Geometry::Point;
  using Rodin::Geometry::Polytope;
  using MeshType = Rodin::Geometry::Mesh<Rodin::Context::Local>;
  using VectorFESType = Rodin::Variational::VectorP1<MeshType>;
  using VectorGridFunctionType = Rodin::Variational::GridFunction<VectorFESType, Rodin::Math::Vector<Real>>;
  using ScalarFESType = Rodin::Variational::P1<Real, MeshType>;
  using ScalarGridFunctionType = Rodin::Variational::GridFunction<ScalarFESType, Rodin::Math::Vector<Real>>;

  /**
   * @brief Pathlines of the ArterialLesion3D output.
   */
  class Pathlines3D
  {
    public:
      struct Config
      {
        std::filesystem::path caseDir;
        std::filesystem::path outStem;
        size_t seeds = 200;
        size_t releases = 8;
        size_t substeps = 4;
        Real transits = 3.0;
        Real horizon = 6.0;     ///< minimum integration time, in cycles
        Real period = 0.0;      ///< 0: read from run.log
        size_t cycles = 1;
      };

      /// @brief XDMF id of a tetrahedron in a mixed topology stream.
      static constexpr std::uint64_t XDMFTetrahedron = 6;
      /// @brief Cosine of the largest angle between a cap normal and the
      ///        axis (10 degrees): the caps are plane cuts normal to the
      ///        vessel, the wall normals are nearly normal to the axis.
      static constexpr Real CapAlignment = 0.985;
      /// @brief Seeds are placed inside this fraction of the inlet radius: the
      ///        no-slip ring carries no flow and its particles would never leave.
      static constexpr Real SeedRadiusFraction = 0.97;
      /// @brief Inward offset of the seeds from the inlet plane, relative to
      ///        the inlet radius, so that the locator finds them in a cell.
      static constexpr Real SeedInwardOffset = 1.0e-3;
      /// @brief A trajectory that moves less than this fraction of the axial
      ///        extent during one period is stagnant and stops ("trapped").
      static constexpr Real StagnationFraction = 1.0e-4;
      /// @brief Relative departure of the wall radius from the parent radius
      ///        that marks the lesion (the same as the postprocessing).
      static constexpr Real LesionDeviation = 0.05;
      /// @brief Sampling of the washout curve, in cycles.
      static constexpr Real WashoutStep = 0.05;
      /// @brief Relative tolerance on coordinates when merging the vertices
      ///        shared by the rank pieces (they come from the same mesh, so
      ///        they match to rounding).
      static constexpr Real VertexMergeTolerance = 1.0e-9;

      explicit Pathlines3D(Config cfg)
        : m_cfg(std::move(cfg))
      {}

      void run()
      {
        readLayout();
        buildMesh();
        loadVelocity();
        findCaps();
        findLesion();
        seed();
        integrate();
        regions();
        write();
      }

    private:
      struct Piece
      {
        std::string meshFile, verticesPath, topologyPath;
        std::vector<std::pair<std::string, std::string>> velocity;   ///< per step: file, path
      };

      struct Cap
      {
        std::vector<Index> faces;
        Math::SpatialVector<Real> centroid;
        Math::SpatialVector<Real> normal;   ///< outward
        Real radius = 0.0;
        Real area = 0.0;
        Real meanFlux = 0.0;                ///< cycle mean of int u . n_out
      };

      struct Sample
      {
        Math::SpatialVector<Real> x;
        Real t;
        Real speed;
      };

      struct Trajectory
      {
        size_t seed = 0, release = 0;
        Real releaseTime = 0.0, residence = 0.0, length = 0.0, maxSpeed = 0.0;
        std::string exit;
        /// Per region (0 the lesion slab, 1 the sac): time inside, and the
        /// times of the first entry and last exit (< 0: never entered).
        std::array<Real, 2> inside{ { 0.0, 0.0 } }, firstIn{ { -1.0, -1.0 } }, lastIn{ { -1.0, -1.0 } };
        std::vector<Sample> samples;
      };

      // ---- the layout of the case: times, pieces, period --------------------
      void readLayout()
      {
        namespace pt = boost::property_tree;
        const auto xdmf = m_cfg.caseDir / "case.xdmf";
        if (!std::filesystem::exists(xdmf))
          throw std::runtime_error("no case.xdmf in " + m_cfg.caseDir.string());
        pt::ptree tree;
        pt::read_xml(xdmf.string(), tree);
        const auto& domain = tree.get_child("Xdmf.Domain");
        // Temporal collection of spatial collections (one per step), each
        // with one grid per rank.
        std::vector<Real> times;
        std::vector<std::vector<Piece>> steps;
        const auto splitRef = [](const std::string& ref) {
          const auto colon = ref.find(':');
          if (colon == std::string::npos)
            throw std::runtime_error("XDMF data item without file:path reference: " + ref);
          std::string file = ref.substr(0, colon), path = ref.substr(colon + 1);
          const auto trim = [](std::string& s) {
            const auto a = s.find_first_not_of(" \n\r\t"), b = s.find_last_not_of(" \n\r\t");
            s = (a == std::string::npos) ? std::string() : s.substr(a, b - a + 1);
          };
          trim(file); trim(path);
          return std::make_pair(file, path);
        };
        const auto visit = [&](const pt::ptree& grid, auto&& self) -> void {
          const std::string type = grid.get<std::string>("<xmlattr>.GridType", "Uniform");
          if (type == "Collection")
          {
            const std::string ct = grid.get<std::string>("<xmlattr>.CollectionType", "Spatial");
            if (ct == "Spatial")
            {
              const auto* t = grid.get_child_optional("Time") ? &grid.get_child("Time") : nullptr;
              times.push_back(t ? t->template get<Real>("<xmlattr>.Value") : Real(times.size()));
              steps.emplace_back();
            }
            for (const auto& [key, child] : grid)
              if (key == "Grid")
                self(child, self);
          }
          else
          {
            Piece piece;
            const auto geo = splitRef(grid.get<std::string>("Geometry.DataItem"));
            piece.meshFile = geo.first;
            piece.verticesPath = geo.second;
            piece.topologyPath = splitRef(grid.get<std::string>("Topology.DataItem")).second;
            const std::string topo = grid.get<std::string>("Topology.<xmlattr>.TopologyType", "");
            if (topo != "Mixed")
              throw std::runtime_error("expected a Mixed topology stream, got " + topo);
            for (const auto& [key, attr] : grid)
            {
              if (key != "Attribute" || attr.template get<std::string>("<xmlattr>.Name", "") != "velocity")
                continue;
              piece.velocity.push_back(splitRef(attr.template get<std::string>("DataItem")));
            }
            if (piece.velocity.empty())
              throw std::runtime_error("a grid of case.xdmf has no velocity attribute");
            if (steps.empty())
              steps.emplace_back(), times.push_back(0.0);
            steps.back().push_back(std::move(piece));
          }
        };
        for (const auto& [key, child] : domain)
          if (key == "Grid")
            visit(child, visit);
        if (steps.size() < 2)
          throw std::runtime_error("case.xdmf holds fewer than two snapshots");

        // The period, from run.log unless given.
        m_period = m_cfg.period;
        if (!(m_period > 0.0))
        {
          std::ifstream log(m_cfg.caseDir / "run.log");
          const std::regex re(R"(\bT=([0-9.eE+-]+) s)");
          for (std::string line; std::getline(log, line);)
          {
            std::smatch m;
            if (line.find("[inflow]") != std::string::npos && std::regex_search(line, m, re))
            {
              m_period = std::stod(m[1]);
              break;
            }
          }
          if (!(m_period > 0.0))
            throw std::runtime_error("period not found in run.log: give -pl_period <s>");
        }
        // Keep the snapshots of the last -pl_cycles cycles; the first kept
        // snapshot is the one closest to t_end - cycles T.
        const Real tEnd = times.back();
        const Real tStart = tEnd - m_cfg.cycles * m_period;
        size_t first = 0;
        for (size_t i = 0; i < times.size(); ++i)
          if (std::abs(times[i] - tStart) < std::abs(times[first] - tStart))
            first = i;
        const Real snapGap = std::abs(times[first] - tStart);
        for (size_t i = first; i < times.size(); ++i)
        {
          m_times.push_back(times[i]);
          m_pieces.push_back(std::move(steps[i]));
        }
        m_span = m_times.back() - m_times.front();
        const size_t ranks = m_pieces.front().size();
        for (const auto& s : m_pieces)
          if (s.size() != ranks)
            throw std::runtime_error("the snapshots do not all have the same number of rank pieces");
        Alert::Info() << "[layout] " << times.size() << " snapshots in case.xdmf, " << ranks
                      << " rank piece(s); using " << m_times.size() << " from t = " << m_times.front()
                      << " to " << tEnd << " s (period " << m_period << " s, span/period = "
                      << m_span / m_period << ", first snapshot " << snapGap
                      << " s from the cycle start)" << Alert::Raise;
        if (m_span < 0.9 * m_period)
          Alert::Warning() << "[layout] less than a cycle of snapshots: the field is extended "
                              "periodically over " << m_span << " s" << Alert::Raise;
      }

      static std::vector<Real> readReal(const std::filesystem::path& file, const std::string& path)
      {
        const auto h = IO::HDF5::File(H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!h)
          throw std::runtime_error("cannot open " + file.string());
        return IO::HDF5::readVectorDataset<Real>(h.get(), path);
      }

      static std::vector<IO::HDF5::U64> readUInt(const std::filesystem::path& file, const std::string& path)
      {
        const auto h = IO::HDF5::File(H5Fopen(file.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!h)
          throw std::runtime_error("cannot open " + file.string());
        return IO::HDF5::readVectorDataset<IO::HDF5::U64>(h.get(), path);
      }

      // ---- the mesh: the union of the owned cells of every piece ------------
      void buildMesh()
      {
        const auto& pieces = m_pieces.front();
        // Global vertex of each (piece, local vertex), merging the shared ones.
        std::map<std::array<std::int64_t, 3>, Index> lookup;
        std::vector<Math::SpatialVector<Real>> vertices;
        std::vector<std::array<Index, 4>> cells;
        m_pieceVertex.assign(pieces.size(), {});
        Real scale = 0.0;
        std::vector<std::vector<Real>> coords(pieces.size());
        for (size_t r = 0; r < pieces.size(); ++r)
        {
          coords[r] = readReal(m_cfg.caseDir / pieces[r].meshFile, pieces[r].verticesPath);
          for (Real c : coords[r])
            scale = std::max(scale, std::abs(c));
        }
        const Real quantum = VertexMergeTolerance * scale;
        for (size_t r = 0; r < pieces.size(); ++r)
        {
          const auto& xs = coords[r];
          if (xs.size() % 3 != 0)
            throw std::runtime_error("vertex array of piece " + std::to_string(r) + " is not 3D");
          const size_t n = xs.size() / 3;
          auto& map = m_pieceVertex[r];
          map.resize(n);
          for (size_t i = 0; i < n; ++i)
          {
            std::array<std::int64_t, 3> key;
            for (size_t c = 0; c < 3; ++c)
              key[c] = static_cast<std::int64_t>(std::llround(xs[3 * i + c] / quantum));
            const auto it = lookup.find(key);
            if (it == lookup.end())
            {
              lookup.emplace(key, vertices.size());
              map[i] = vertices.size();
              vertices.push_back(Math::SpatialVector<Real>{{ xs[3 * i], xs[3 * i + 1], xs[3 * i + 2] }});
            }
            else
              map[i] = it->second;
          }
          const auto topo = readUInt(m_cfg.caseDir / pieces[r].meshFile, pieces[r].topologyPath);
          for (size_t k = 0; k < topo.size();)
          {
            if (topo[k] != XDMFTetrahedron)
              throw std::runtime_error("only tetrahedra are supported (XDMF type "
                                       + std::to_string(topo[k]) + " found)");
            if (k + 4 >= topo.size())
              throw std::runtime_error("truncated topology stream");
            cells.push_back({ map[topo[k + 1]], map[topo[k + 2]], map[topo[k + 3]], map[topo[k + 4]] });
            k += 5;
          }
        }
        MeshType::Builder builder;
        builder.initialize(3).nodes(vertices.size()).reserve(0, vertices.size()).reserve(3, cells.size());
        for (const auto& x : vertices)
          builder.vertex(x);
        for (const auto& c : cells)
          builder.polytope(Polytope::Type::Tetrahedron, { c[0], c[1], c[2], c[3] });
        m_mesh = builder.finalize();
        const size_t D = m_mesh.getDimension();
        m_mesh.getConnectivity().compute(D, D);
        m_mesh.getConnectivity().compute(D, 0);
        m_mesh.getConnectivity().compute(D, D - 1);
        m_mesh.getConnectivity().compute(D - 1, D);
        m_mesh.getConnectivity().compute(D - 1, 0);
        m_mesh.getConnectivity().compute(D - 1, 1);
        m_mesh.getConnectivity().compute(1, 0);
        // The axis and extent of the vessel: principal direction of the
        // vertex cloud (the pipe is long against its diameter).
        Math::SpatialVector<Real> mean = Math::SpatialVector<Real>::Zero(3);
        for (const auto& x : vertices)
          mean += x;
        mean /= static_cast<Real>(vertices.size());
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (const auto& x : vertices)
        {
          const Eigen::Vector3d d(x(0) - mean(0), x(1) - mean(1), x(2) - mean(2));
          cov += d * d.transpose();
        }
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
        const Eigen::Vector3d a = eig.eigenvectors().col(2);
        m_axis = Math::SpatialVector<Real>{{ a(0), a(1), a(2) }};
        Real lo = std::numeric_limits<Real>::max(), hi = std::numeric_limits<Real>::lowest();
        for (const auto& x : vertices)
        {
          const Real s = x.dot(m_axis);
          lo = std::min(lo, s);
          hi = std::max(hi, s);
        }
        m_axialExtent = hi - lo;
        Alert::Info() << "[mesh] " << vertices.size() << " vertices, " << cells.size()
                      << " tetrahedra from " << pieces.size() << " piece(s); axis ("
                      << m_axis(0) << ", " << m_axis(1) << ", " << m_axis(2) << "), extent "
                      << m_axialExtent << " m" << Alert::Raise;
      }

      // ---- the velocity snapshots ----------------------------------------------
      void loadVelocity()
      {
        m_fes = std::make_unique<VectorFESType>(m_mesh, 3);
        m_u.reserve(m_times.size());
        for (size_t s = 0; s < m_times.size(); ++s)
        {
          auto& u = m_u.emplace_back(*m_fes);
          u.getData().setZero();
          for (size_t r = 0; r < m_pieces[s].size(); ++r)
          {
            const auto& [file, path] = m_pieces[s][r].velocity.front();
            const auto vals = readReal(m_cfg.caseDir / file, path);
            const auto& map = m_pieceVertex[r];
            if (vals.size() != 3 * map.size())
              throw std::runtime_error("velocity of piece " + std::to_string(r) + " at step "
                                       + std::to_string(s) + " does not match its vertices");
            for (size_t i = 0; i < map.size(); ++i)
              for (size_t c = 0; c < 3; ++c)
                u.getData()(m_fes->getGlobalIndex({ 0, map[i] }, c)) = vals[3 * i + c];
          }
        }
        Alert::Info() << "[velocity] " << m_u.size() << " snapshots loaded" << Alert::Raise;
      }

      /// @brief Velocity frozen at time t (periodic), linear between snapshots.
      struct FrozenVelocity
      {
        const VectorGridFunctionType* a;
        const VectorGridFunctionType* b;
        Real theta;
        Math::SpatialVector<Real> operator()(const Point& p) const
        {
          const auto va = a->getValue(p);
          const auto vb = b->getValue(p);
          Math::SpatialVector<Real> out(3);
          for (Index c = 0; c < 3; ++c)
            out(c) = (1.0 - theta) * va(c) + theta * vb(c);
          return out;
        }
      };

      FrozenVelocity velocityAt(Real t) const
      {
        // Periodic time within the snapshot span.
        Real tau = std::fmod(t - m_times.front(), m_span);
        if (tau < 0.0)
          tau += m_span;
        const Real tt = m_times.front() + tau;
        size_t k = 0;
        while (k + 2 < m_times.size() && m_times[k + 1] <= tt)
          ++k;
        const Real dt = m_times[k + 1] - m_times[k];
        const Real theta = dt > 0.0 ? std::clamp((tt - m_times[k]) / dt, 0.0, 1.0) : 0.0;
        return FrozenVelocity{ &m_u[k], &m_u[k + 1], theta };
      }

      // ---- the caps and the inlet -------------------------------------------------
      void findCaps()
      {
        const size_t D = m_mesh.getDimension();
        const Variational::BoundaryNormal normal(m_mesh);
        const Polytope::Traits triangle(Polytope::Type::Triangle);
        // Boundary faces aligned with the axis, grouped by their side.
        struct FaceData { Index index; Real s; Math::SpatialVector<Real> n; Real area; Math::SpatialVector<Real> c; };
        std::vector<FaceData> aligned;
        Real sMin = std::numeric_limits<Real>::max(), sMax = std::numeric_limits<Real>::lowest();
        for (auto it = m_mesh.getBoundary(); it; ++it)
        {
          const Point p(*it, triangle.getCentroid());
          const auto nv = normal.getValue(p);
          Math::SpatialVector<Real> n(3);
          for (Index c = 0; c < 3; ++c)
            n(c) = nv(c);
          if (std::abs(n.dot(m_axis)) < CapAlignment)
            continue;
          FaceData f{ it->getIndex(), 0.0, n, it->getMeasure(), p.getPhysicalCoordinates() };
          f.s = f.c.dot(m_axis);
          sMin = std::min(sMin, f.s);
          sMax = std::max(sMax, f.s);
          aligned.push_back(std::move(f));
        }
        if (aligned.empty())
          throw std::runtime_error("no boundary face is aligned with the axis: no caps found");
        const Real sMid = 0.5 * (sMin + sMax);
        Cap caps[2];
        for (const auto& f : aligned)
        {
          Cap& cap = caps[f.s < sMid ? 0 : 1];
          if (cap.faces.empty())
          {
            cap.centroid = Math::SpatialVector<Real>::Zero(3);
            cap.normal = Math::SpatialVector<Real>::Zero(3);
          }
          cap.faces.push_back(f.index);
          cap.centroid += f.area * f.c;
          cap.normal += f.area * f.n;
          cap.area += f.area;
        }
        for (Cap& cap : caps)
        {
          if (cap.faces.empty())
            throw std::runtime_error("only one cap found: the mesh has no inlet/outlet pair along the axis");
          cap.centroid /= cap.area;
          cap.normal.normalize();
          // Radius: the farthest vertex of the cap faces from its centroid.
          for (Index f : cap.faces)
            for (auto v = m_mesh.getFace(f)->getVertex(); v; ++v)
              cap.radius = std::max(cap.radius, (m_mesh.getVertexCoordinates(v->getIndex()) - cap.centroid).norm());
          // Cycle-mean flux: trapezoid in time of int_cap u . n_out, u at the
          // face centroid (P1: the mean of its vertices).
          Real acc = 0.0;
          for (size_t s = 0; s + 1 < m_times.size(); ++s)
          {
            Real q0 = 0.0, q1 = 0.0;
            for (Index f : cap.faces)
            {
              const auto face = m_mesh.getFace(f);
              const Point p(*face, triangle.getCentroid());
              const auto u0 = m_u[s].getValue(p);
              const auto u1 = m_u[s + 1].getValue(p);
              Real d0 = 0.0, d1 = 0.0;
              for (Index c = 0; c < 3; ++c)
              {
                d0 += u0(c) * cap.normal(c);
                d1 += u1(c) * cap.normal(c);
              }
              const Real area = face->getMeasure();
              q0 += d0 * area;
              q1 += d1 * area;
            }
            acc += 0.5 * (q0 + q1) * (m_times[s + 1] - m_times[s]);
          }
          cap.meanFlux = acc / m_span;
        }
        const bool firstIsInlet = caps[0].meanFlux < caps[1].meanFlux;
        m_inlet = firstIsInlet ? caps[0] : caps[1];
        m_outlet = firstIsInlet ? caps[1] : caps[0];
        if (!(m_inlet.meanFlux < 0.0))
          throw std::runtime_error("neither cap has a negative mean flux: no inflow found");
        m_meanSpeed = -m_inlet.meanFlux / m_inlet.area;
        for (Index f : m_outlet.faces)
          m_outletFaces.insert(f);
        Alert::Info() << "[inlet] centroid (" << m_inlet.centroid(0) << ", " << m_inlet.centroid(1)
                      << ", " << m_inlet.centroid(2) << ") m, radius " << m_inlet.radius
                      << " m, " << m_inlet.faces.size() << " faces, mean flow "
                      << -m_inlet.meanFlux << " m^3/s, mean speed " << m_meanSpeed
                      << " m/s; outlet " << m_outlet.faces.size() << " faces, mean flux "
                      << m_outlet.meanFlux << " m^3/s" << Alert::Raise;
        (void)D;
      }

      // ---- the lesion ---------------------------------------------------------------
      void findLesion()
      {
        std::set<Index> capFaces(m_outletFaces.begin(), m_outletFaces.end());
        capFaces.insert(m_inlet.faces.begin(), m_inlet.faces.end());
        std::set<Index> wallVertices;
        for (auto it = m_mesh.getBoundary(); it; ++it)
          if (!capFaces.count(it->getIndex()))
            for (auto v = it->getVertex(); v; ++v)
              wallVertices.insert(v->getIndex());
        std::vector<std::pair<Real, Real>> sr;   // (axial, radial) of each wall vertex
        sr.reserve(wallVertices.size());
        for (Index v : wallVertices)
          sr.push_back(axialRadial(m_mesh.getVertexCoordinates(v)));
        std::vector<Real> r(sr.size());
        for (size_t i = 0; i < sr.size(); ++i)
          r[i] = sr[i].second;
        std::nth_element(r.begin(), r.begin() + r.size() / 2, r.end());
        m_parentRadius = r[r.size() / 2];
        Real lo = std::numeric_limits<Real>::max(), hi = std::numeric_limits<Real>::lowest(), dev = 0.0;
        size_t n = 0;
        for (const auto& [sv, rv] : sr)
          if (std::abs(rv - m_parentRadius) > LesionDeviation * m_parentRadius)
          {
            lo = std::min(lo, sv);
            hi = std::max(hi, sv);
            dev += rv - m_parentRadius;
            ++n;
          }
        if (n == 0)
        {
          Alert::Warning() << "[lesion] no wall departs from the parent radius: no lesion statistics" << Alert::Raise;
          m_lesionLo = 1.0;
          m_lesionHi = 0.0;
          return;
        }
        m_lesionLo = lo;
        m_lesionHi = hi;
        m_dilation = dev > 0.0;
        const Real D = 2.0 * m_parentRadius;
        Alert::Info() << "[lesion] parent D = " << D << " m; " << (m_dilation ? "dilation" : "constriction")
                      << " over s/D in [" << lo / D << ", " << hi / D << "] from the inlet (length "
                      << (hi - lo) / D << " D)" << (m_dilation ? "; sac: r > R inside it" : "") << Alert::Raise;
      }

      /// @brief Axial position from the inlet centroid and distance to the axis.
      std::pair<Real, Real> axialRadial(const Math::SpatialVector<Real>& x) const
      {
        const Math::SpatialVector<Real> d = x - m_inlet.centroid;
        const Real sv = d.dot(m_axis);
        return { sv, (d - sv * m_axis).norm() };
      }

      // ---- the seeds ---------------------------------------------------------------
      void seed()
      {
        // An orthonormal basis of the inlet plane.
        const Math::SpatialVector<Real>& n = m_inlet.normal;
        Math::SpatialVector<Real> e1(3), e2(3);
        {
          const Eigen::Vector3d nn(n(0), n(1), n(2));
          const Eigen::Vector3d helper = std::abs(nn(0)) < 0.9 ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
          const Eigen::Vector3d a = (helper - helper.dot(nn) * nn).normalized();
          const Eigen::Vector3d b = nn.cross(a);
          for (Index c = 0; c < 3; ++c)
          {
            e1(c) = a(c);
            e2(c) = b(c);
          }
        }
        m_locator = std::make_unique<Location::AABB<MeshType>>(m_mesh);
        const Real golden = M_PI * (3.0 - std::sqrt(5.0));
        const Real rMax = SeedRadiusFraction * m_inlet.radius;
        size_t misses = 0;
        for (size_t k = 0; k < m_cfg.seeds; ++k)
        {
          // Fibonacci disc: near-uniform density over the area.
          const Real r = rMax * std::sqrt((k + 0.5) / static_cast<Real>(m_cfg.seeds));
          const Real phi = k * golden;
          Math::SpatialVector<Real> x = m_inlet.centroid + r * (std::cos(phi) * e1 + std::sin(phi) * e2)
            - (SeedInwardOffset * m_inlet.radius) * n;
          if (const auto p = m_locator->locate(x))
            m_seeds.push_back(*p);
          else
            ++misses;
        }
        if (m_seeds.empty())
          throw std::runtime_error("no seed could be located inside the mesh");
        Alert::Info() << "[seeds] " << m_seeds.size() << " seeds over the inlet (" << misses
                      << " outside the mesh), " << m_cfg.releases << " release phases" << Alert::Raise;
      }

      // ---- the integration --------------------------------------------------------
      void integrate()
      {
        const Real transit = m_axialExtent / m_meanSpeed;
        const Real tMax = std::max(m_cfg.transits * transit, m_cfg.horizon * m_period);
        // Sub-step: a fraction of the (smallest) snapshot interval.
        Real snapDt = m_span;
        for (size_t s = 0; s + 1 < m_times.size(); ++s)
          snapDt = std::min(snapDt, m_times[s + 1] - m_times[s]);
        const Real dtau = snapDt / static_cast<Real>(std::max<size_t>(1, m_cfg.substeps));
        const Real stagnation = StagnationFraction * m_axialExtent;
        const size_t stagnationSteps = static_cast<size_t>(std::ceil(m_period / dtau));
        Alert::Info() << "[integrate] mean transit L/U = " << transit << " s, t_max = " << tMax
                      << " s (" << tMax / m_period << " cycles), sub-step " << dtau << " s (" << m_cfg.substeps
                      << " per snapshot interval)" << Alert::Raise;
        // The tracer composes an operand with the flow map; only the traced
        // point is used here, the operand is a constant.
        const Variational::RealFunction zero([](const Point&) { return Real(0); });
        const Math::RungeKutta::RK4 rk4;
        size_t exited = 0, trapped = 0, wall = 0;
        m_trajectories.reserve(m_seeds.size() * m_cfg.releases);
        for (size_t m = 0; m < m_cfg.releases; ++m)
        {
          const Real t0 = m_times.front() + m * m_period / static_cast<Real>(m_cfg.releases);
          for (size_t k = 0; k < m_seeds.size(); ++k)
          {
            Trajectory traj;
            traj.seed = k;
            traj.release = m;
            traj.releaseTime = t0;
            std::optional<Point> p;
            p.emplace(m_seeds[k]);
            Real t = t0;
            Math::SpatialVector<Real> x = p->getPhysicalCoordinates();
            {
              const auto v = velocityAt(t)(*p);
              traj.samples.push_back({ x, t, v.norm() });
            }
            Math::SpatialVector<Real> xStag = x;
            size_t sinceStag = 0;
            while (t - t0 < tMax)
            {
              const Real h = std::min(dtau, tMax - (t - t0));
              const FrozenVelocity vel = velocityAt(t + 0.5 * h);
              const Variational::Flow flow(h, zero, vel, rk4);
              const auto trace = flow.trace(*p);
              const Point& q = trace.getPoint();
              const Math::SpatialVector<Real> xq = q.getPhysicalCoordinates();
              const Real consumed = h - trace.getTime();
              t += std::max(consumed, Real(0));
              traj.length += (xq - x).norm();
              const auto v = vel(q);
              traj.maxSpeed = std::max(traj.maxSpeed, v.norm());
              traj.samples.push_back({ xq, t, v.norm() });
              x = xq;
              if (trace.exited())
              {
                // The trace stopped on a boundary face: the outlet, or the
                // wall (a limit of the tracing near the no-slip wall).
                const auto& poly = q.getPolytope();
                bool outlet = false;
                if (poly.getDimension() == m_mesh.getDimension() - 1)
                  outlet = m_outletFaces.count(poly.getIndex()) > 0;
                else
                  outlet = std::abs((xq - m_outlet.centroid).dot(m_outlet.normal)) < SeedInwardOffset * m_outlet.radius * 10.0;
                traj.exit = outlet ? "outlet" : "wall";
                (outlet ? exited : wall)++;
                break;
              }
              p.emplace(q);
              if (++sinceStag >= stagnationSteps)
              {
                if ((x - xStag).norm() < stagnation)
                {
                  traj.exit = "trapped";
                  ++trapped;
                  break;
                }
                xStag = x;
                sinceStag = 0;
              }
            }
            if (traj.exit.empty())
            {
              traj.exit = "trapped";
              ++trapped;
            }
            traj.residence = t - t0;
            m_trajectories.push_back(std::move(traj));
          }
          Alert::Info() << "[integrate] release " << m + 1 << "/" << m_cfg.releases << " at t = " << t0
                        << " s done: so far outlet " << exited << ", trapped " << trapped
                        << ", wall " << wall << Alert::Raise;
        }
      }

      // ---- residence in the lesion ------------------------------------------------
      void regions()
      {
        if (!(m_lesionHi > m_lesionLo))
          return;
        const auto inRegion = [this](const Math::SpatialVector<Real>& x, size_t k) {
          const auto [sv, rv] = axialRadial(x);
          const bool slab = sv >= m_lesionLo && sv <= m_lesionHi;
          return k == 0 ? slab : (slab && m_dilation && rv > m_parentRadius);
        };
        for (auto& tr : m_trajectories)
          for (size_t k = 0; k < 2; ++k)
          {
            bool prev = inRegion(tr.samples.front().x, k);
            if (prev)
              tr.firstIn[k] = tr.lastIn[k] = tr.samples.front().t;
            for (size_t i = 1; i < tr.samples.size(); ++i)
            {
              const bool cur = inRegion(tr.samples[i].x, k);
              const Real dt = tr.samples[i].t - tr.samples[i - 1].t;
              // Trapezoid of the indicator: half a step at a crossing.
              tr.inside[k] += 0.5 * dt * (Real(prev) + Real(cur));
              if (cur)
              {
                if (tr.firstIn[k] < 0.0)
                  tr.firstIn[k] = tr.samples[i].t;
                tr.lastIn[k] = tr.samples[i].t;
              }
              prev = cur;
            }
          }
      }

      // ---- output -----------------------------------------------------------------
      void write()
      {
        const auto stem = m_cfg.outStem.empty() ? m_cfg.caseDir / "pathlines" : m_cfg.outStem;
        // CSV, one row per trajectory.
        {
          std::ofstream csv(stem.string() + ".csv");
          csv << "release,releaseTime,seed,x0,y0,z0,exit,residence,length,maxSpeed,samples,"
                 "lesionTime,lesionSpan,sacTime,sacSpan\n";
          for (const auto& tr : m_trajectories)
          {
            const auto& x0 = tr.samples.front().x;
            csv << tr.release << ',' << tr.releaseTime << ',' << tr.seed << ',' << x0(0) << ',' << x0(1)
                << ',' << x0(2) << ',' << tr.exit << ',' << tr.residence << ',' << tr.length << ','
                << tr.maxSpeed << ',' << tr.samples.size();
            for (size_t k = 0; k < 2; ++k)
              csv << ',' << tr.inside[k] << ',' << (tr.firstIn[k] < 0.0 ? 0.0 : tr.lastIn[k] - tr.firstIn[k]);
            csv << '\n';
          }
        }
        // Residence-time statistics per release phase (the exited ones), and
        // the fractions of each fate.
        {
          std::ofstream st(stem.string() + "_residence.csv");
          st << "release,releaseTime,phase,n,outlet,trapped,wall,residenceMean,residenceMedian,residenceP90,residenceMax\n";
          for (size_t m = 0; m < m_cfg.releases; ++m)
          {
            std::vector<Real> res;
            size_t n = 0, out = 0, trap = 0, wal = 0;
            Real t0 = 0.0;
            for (const auto& tr : m_trajectories)
            {
              if (tr.release != m)
                continue;
              ++n;
              t0 = tr.releaseTime;
              if (tr.exit == "outlet")
              {
                ++out;
                res.push_back(tr.residence);
              }
              else if (tr.exit == "trapped")
                ++trap;
              else
                ++wal;
            }
            std::sort(res.begin(), res.end());
            const auto q = [&](Real f) { return res.empty() ? 0.0 : res[std::min(res.size() - 1, static_cast<size_t>(f * res.size()))]; };
            const Real mean = res.empty() ? 0.0 : std::accumulate(res.begin(), res.end(), 0.0) / res.size();
            st << m << ',' << t0 << ',' << std::fmod(t0 - m_times.front(), m_period) / m_period << ',' << n << ','
               << out << ',' << trap << ',' << wal << ',' << mean << ',' << q(0.5) << ',' << q(0.9) << ','
               << (res.empty() ? 0.0 : res.back()) << '\n';
          }
        }
        writeLesion(stem);
        // XDMF: the pathlines as a segment mesh in 3D with nodal fields.
        size_t nv = 0, ne = 0;
        for (const auto& tr : m_trajectories)
        {
          nv += tr.samples.size();
          ne += tr.samples.size() - 1;
        }
        MeshType::Builder builder;
        builder.initialize(3).nodes(nv).reserve(0, nv).reserve(1, ne);
        std::vector<Real> time, speed, seedId, releaseId, residence, age, lesionT, sacT;
        time.reserve(nv); speed.reserve(nv); seedId.reserve(nv); releaseId.reserve(nv);
        residence.reserve(nv); age.reserve(nv);
        Index next = 0;
        for (const auto& tr : m_trajectories)
        {
          for (size_t i = 0; i < tr.samples.size(); ++i)
          {
            const auto& s = tr.samples[i];
            builder.vertex(s.x);
            time.push_back(s.t);
            speed.push_back(s.speed);
            seedId.push_back(static_cast<Real>(tr.seed));
            releaseId.push_back(static_cast<Real>(tr.release));
            residence.push_back(tr.residence);
            age.push_back(s.t - tr.releaseTime);
            lesionT.push_back(tr.inside[0] / m_period);
            sacT.push_back(tr.inside[1] / m_period);
            if (i > 0)
              builder.polytope(Polytope::Type::Segment, { next - 1, next });
            ++next;
          }
        }
        MeshType lines = builder.finalize();
        lines.getConnectivity().compute(1, 0);
        ScalarFESType fes(lines);
        ScalarGridFunctionType gTime(fes), gSpeed(fes), gSeed(fes), gRelease(fes), gResidence(fes), gAge(fes), gLesion(fes), gSac(fes);
        const auto fill = [&](ScalarGridFunctionType& g, const std::vector<Real>& v) {
          for (Index i = 0; i < static_cast<Index>(v.size()); ++i)
            g.getData()(fes.getGlobalIndex({ 0, i }, 0)) = v[i];
        };
        fill(gTime, time); fill(gSpeed, speed); fill(gSeed, seedId); fill(gRelease, releaseId);
        fill(gResidence, residence); fill(gAge, age); fill(gLesion, lesionT); fill(gSac, sacT);
        IO::XDMF xdmf(boost::filesystem::path(stem.string()));
        xdmf.setMesh(lines)
            .add("time", gTime)
            .add("speed", gSpeed)
            .add("seed", gSeed)
            .add("release", gRelease)
            .add("residence", gResidence)
            .add("age", gAge)
            .add("lesionTime_T", gLesion)
            .add("sacTime_T", gSac);
        xdmf.write(0.0);
        xdmf.close();
        Alert::Info() << "[write] " << m_trajectories.size() << " pathlines, " << nv << " samples: "
                      << stem.string() << ".xdmf, .csv, _residence.csv, _lesion.csv, _washout.csv" << Alert::Raise;
      }

      /// @brief Washout of each region over the trajectories that entered it
      ///        (all release phases together), and its summary statistics.
      void writeLesion(const std::filesystem::path& stem) const
      {
        if (!(m_lesionHi > m_lesionLo))
          return;
        const char* names[2] = { "lesion", "sac" };
        std::ofstream sum(stem.string() + "_lesion.csv");
        sum << "region,period,released,entered,censored,timeMedian_T,timeP90_T,timeMean_T,"
               "spanMedian_T,spanP90_T,washoutHalf_T,remaining1T,remaining2T,remaining4T\n";
        std::ofstream wash(stem.string() + "_washout.csv");
        wash << "region,tau_T,remaining\n";
        for (size_t k = 0; k < 2; ++k)
        {
          if (k == 1 && !m_dilation)
            continue;
          // span: time from the first entry to the last exit, in cycles;
          // censored: stopped while able to return (not at the outlet).
          std::vector<Real> time, span;
          size_t released = 0, censored = 0;
          for (const auto& tr : m_trajectories)
          {
            if (tr.exit == "wall")
              continue;   // a tracing limit, not a fate
            ++released;
            if (tr.firstIn[k] < 0.0)
              continue;
            time.push_back(tr.inside[k] / m_period);
            const bool cens = tr.exit != "outlet";
            censored += cens;
            const Real last = cens ? tr.samples.back().t : tr.lastIn[k];
            span.push_back((last - tr.firstIn[k]) / m_period);
          }
          if (time.empty())
            continue;
          std::sort(time.begin(), time.end());
          std::sort(span.begin(), span.end());
          const auto q = [](const std::vector<Real>& v, Real f) {
            return v[std::min(v.size() - 1, static_cast<size_t>(f * v.size()))];
          };
          const auto remaining = [&span](Real tau) {
            return static_cast<Real>(span.end() - std::upper_bound(span.begin(), span.end(), tau))
              / static_cast<Real>(span.size());
          };
          const Real mean = std::accumulate(time.begin(), time.end(), 0.0) / time.size();
          // Half-life: the first tau at which half of the entered are gone
          // for good (nan if the horizon is reached first).
          const Real half = remaining(span.back()) <= 0.5 ? q(span, 0.5) : std::numeric_limits<Real>::quiet_NaN();
          sum << names[k] << ',' << m_period << ',' << released << ',' << time.size() << ','
              << static_cast<Real>(censored) / time.size() << ',' << q(time, 0.5) << ',' << q(time, 0.9) << ','
              << mean << ',' << q(span, 0.5) << ',' << q(span, 0.9) << ',' << half << ','
              << remaining(1.0) << ',' << remaining(2.0) << ',' << remaining(4.0) << '\n';
          const Real tauMax = span.back();
          for (Real tau = 0.0; tau <= tauMax + 1.0e-12; tau += WashoutStep)
            wash << names[k] << ',' << tau << ',' << remaining(tau) << '\n';
          Alert::Info() << "[" << names[k] << "] " << time.size() << " of " << released
                        << " trajectories entered (" << 100.0 * censored / time.size()
                        << " % censored); time inside: median " << q(time, 0.5) << " T, p90 " << q(time, 0.9)
                        << " T; washout half-life " << half << " T, remaining after 1 T "
                        << remaining(1.0) << Alert::Raise;
        }
      }

      Config m_cfg;
      std::vector<Real> m_times;
      std::vector<std::vector<Piece>> m_pieces;
      std::vector<std::vector<Index>> m_pieceVertex;
      Real m_period = 0.0, m_span = 0.0;
      MeshType m_mesh;
      Math::SpatialVector<Real> m_axis;
      Real m_axialExtent = 0.0;
      std::unique_ptr<VectorFESType> m_fes;
      std::vector<VectorGridFunctionType> m_u;
      Cap m_inlet, m_outlet;
      std::set<Index> m_outletFaces;
      Real m_meanSpeed = 0.0;
      std::unique_ptr<Location::AABB<MeshType>> m_locator;
      std::vector<Point> m_seeds;
      std::vector<Trajectory> m_trajectories;
      Real m_parentRadius = 0.0, m_lesionLo = 1.0, m_lesionHi = 0.0;
      bool m_dilation = false;
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Examples::ViscoelasticFluids;
  Pathlines3D::Config cfg;
  for (int i = 1; i < argc; ++i)
  {
    const std::string key = argv[i];
    const auto value = [&]() -> std::string {
      if (i + 1 >= argc)
        throw std::runtime_error("missing value for " + key);
      return argv[++i];
    };
    if (key == "-pl_case") cfg.caseDir = value();
    else if (key == "-pl_out") cfg.outStem = value();
    else if (key == "-pl_seeds") cfg.seeds = std::stoul(value());
    else if (key == "-pl_releases") cfg.releases = std::stoul(value());
    else if (key == "-pl_substeps") cfg.substeps = std::stoul(value());
    else if (key == "-pl_transits") cfg.transits = std::stod(value());
    else if (key == "-pl_horizon") cfg.horizon = std::stod(value());
    else if (key == "-pl_period") cfg.period = std::stod(value());
    else if (key == "-pl_cycles") cfg.cycles = std::stoul(value());
    else if (key == "-h" || key == "--help")
    {
      std::cout << "Pathlines3D -pl_case <dir> [-pl_out stem] [-pl_seeds N] [-pl_releases M]\n"
                   "            [-pl_substeps k] [-pl_transits n] [-pl_horizon cycles] [-pl_period T] [-pl_cycles c]\n";
      return 0;
    }
    else
      throw std::runtime_error("unknown option " + key);
  }
  if (cfg.caseDir.empty())
  {
    std::cerr << "Pathlines3D: -pl_case <dir> is required\n";
    return 1;
  }
  Pathlines3D(cfg).run();
  return 0;
}
