/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EIKONAL_FIM_H
#define RODIN_EIKONAL_FIM_H

/**
 * @file
 * @brief Native first-order-cell fast iterative arrival times.
 */

#include <algorithm>
#include <atomic>
#include <exception>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Configure.h"
#include "Metric.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/PolytopeIterator.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Eikonal
{
  /**
   * @brief Fast iterative arrival times with native-cell semi-Lagrangian updates.
   * @tparam Solution Arrival-time grid function.
   * @tparam Coefficient Scalar speed or squared-velocity tensor function.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref FIM "FIM<GridFunction<P1<Real, Mesh<Context::Local>>, Data>, Coefficient>" | Local scalar vertex space on native first-order straight convex cells. |
   */
  template <class Solution, class Coefficient>
  class FIM;

  /**
   * @brief Native-cell fast iterative method on local real vertex spaces.
   * @tparam Data Mutable arrival-time coefficient storage.
   * @tparam Coefficient Callable scalar speed or ambient SPD squared-velocity tensor.
   *
   * ## Mathematical formulation
   *
   * A scalar coefficient solves @f$c|\nabla T|=1@f$; a matrix coefficient
   * solves @f$\sqrt{\nabla T^T M\nabla T}=1@f$. Seeds have time zero.
   * Scalar speeds must be finite and positive. Tensor coefficients must be
   * finite, symmetric positive definite and have the physical embedding size.
   * On embedded manifolds the path metric is the restriction of @f$M^{-1}@f$
   * to physical tangents (a tangential velocity tensor can be completed by
   * any decoupled positive normal block).
   *
   * ## Architecture
   *
   * A pending solve samples the coefficient once at each cell's reference
   * centroid, before starting the active-list iteration. Scalars are validated
   * and tensor metrics factored once. Samples remain fixed throughout the solve,
   * including all OpenMP workers. Reseed after changing prescribed coefficient
   * data or geometry; the next solve rebuilds every cell sample.
   *
   * At a target physical vertex @f$x_i@f$, the candidate from a native boundary
   * entity @f$F@f$ not containing that vertex is
   * @f[
   * \min_{q\in\widehat F}\left[
   *   \sum_j N_j(q)T_j+
   *   \sqrt{d(q)^T M_K^{-1}d(q)}\right],\qquad
   * d(q)=\sum_jN_j(q)x_j-x_i.
   * @f]
   * Scalar speeds use @f$|d|/c_K@f$. Physical boundary points use the existing
   * native P1 trace basis. Lower-dimensional entities provide candidates while
   * higher-dimensional traces contain infinite times. No simplex subdivision
   * or path quadrature is used.
   *
   * Geometry is restricted to native first-order straight convex cells with
   * planar faces (planar cells for embedded surfaces). Higher-order charts are
   * rejected, including higher-order representations of straight cells.
   * Validity, convexity and planar-face conformity are caller preconditions;
   * first-order multilinear geometry alone does not certify these properties.
   * Spatially varying coefficients are explicitly approximated as cellwise
   * constants; they are not integrated along a path. Face minimization remains
   * a deterministic multistart projected pattern search, with no certificate
   * of global minimization on arbitrary native quadrilateral traces.
   *
   * OpenMP is opt-in through setThreadCount(). Dependency levels preserve
   * serial co-cell read/write ordering and active-list reconstruction, giving
   * identical iterates and iteration counts. Workers only read immutable cell
   * coefficients; user coefficient callbacks run serially during preparation.
   * Worker exceptions are rethrown after the parallel region.
   * Geometry and coefficients must remain fixed during a solve; coefficient callables must be deterministic.
   * Incidences @f$0\to D@f$ and @f$D\to0@f$ must be computed explicitly.
   * Valid, regular, mutually conforming cell maps are a caller precondition.
   */
  template <class Data, class Coefficient>
  class FIM<Variational::GridFunction<
              Variational::P1<Real, Geometry::Mesh<Context::Local>>, Data>,
    Coefficient> final
  {
    public:
      /// @brief Local mesh carrying the native cells.
      using Mesh = Geometry::Mesh<Context::Local>;

      /// @brief Scalar vertex finite element space.
      using FES = Variational::P1<Real, Mesh>;

      /// @brief Borrowed mutable arrival-time field.
      using SolutionType = Variational::GridFunction<FES, Data>;

      /// @brief Owned coefficient callable or Rodin function expression.
      using CoefficientType = Coefficient;

    private:
      using CoefficientValue =
        std::remove_cvref_t<decltype(std::declval<const CoefficientType&>()(
          std::declval<const Geometry::Point&>()))>;
      using CellCoefficient =
        std::conditional_t<std::is_arithmetic_v<CoefficientValue>, Real, Metric>;

      /// @brief Absolute time decrease below which an active update is converged.
      static constexpr Real DefaultTolerance = Real(1e-10);

      /// @brief Reference-coordinate pattern-search resolution; a heuristic policy.
      static constexpr Real DefaultLocalTolerance = Real(1e-5);

      /// @brief Safety limit on active-list rounds; exhaustion is reported as failure.
      static constexpr size_t DefaultMaxIterations = 10000;

      /// @brief Safety limit on search steps per starting point.
      static constexpr size_t MaxSearchIterations = 512;

      /// @brief Initial reference search step; all nonpoint faces use unit reference edges.
      static constexpr Real InitialSearchStep = Real(0.5);

    public:
      /**
       * @brief Binds the output field and owns the coefficient.
       * @param arrival Mutable arrival-time field, overwritten by seed().
       * @param coefficient Positive scalar speed or SPD squared-velocity tensor.
       */
      FIM(SolutionType& arrival, const CoefficientType& coefficient)
        : m_arrival(arrival),
          m_coefficient(coefficient),
          m_tolerance(DefaultTolerance),
          m_localTolerance(DefaultLocalTolerance),
          m_maxIterations(DefaultMaxIterations),
          m_iterations(0),
          m_threadCount(1),
          m_seeded(false)
      {}

      /**
       * @brief Sets the absolute update tolerance in arrival-time units.
       * @param tolerance Finite positive tolerance.
       * @returns This solver.
       */
      FIM& setTolerance(Real tolerance)
      {
        if (!(tolerance > 0) || !std::isfinite(tolerance))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a finite positive arrival-time tolerance." << Alert::Raise;
        }
        m_tolerance = tolerance;
        return *this;
      }

      /**
       * @brief Sets the reference-coordinate pattern-search resolution.
       * @param tolerance Finite resolution strictly between zero and one.
       * @returns This solver.
       */
      FIM& setLocalTolerance(Real tolerance)
      {
        if (!(tolerance > 0 && tolerance < 1) || !std::isfinite(tolerance))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a local tolerance strictly between zero and one."
            << Alert::Raise;
        }
        m_localTolerance = tolerance;
        return *this;
      }

      /**
       * @brief Sets the active-list round limit.
       * @param iterations Positive maximum number of rounds.
       * @returns This solver.
       */
      FIM& setMaxIterations(size_t iterations)
      {
        if (iterations == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a positive iteration limit." << Alert::Raise;
        }
        m_maxIterations = iterations;
        return *this;
      }

      /**
       * @brief Sets the OpenMP worker count, defaulting to serial execution.
       * @param threads Positive worker count; values above one require an OpenMP build.
       * @returns This solver.
       */
      FIM& setThreadCount(size_t threads)
      {
        if (threads == 0 ||
          threads > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected a positive representable OpenMP thread count." << Alert::Raise;
        }
#ifndef RODIN_USE_OPENMP
        if (threads != 1)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Multiple FIM threads require an OpenMP-enabled build." << Alert::Raise;
        }
#endif
        m_threadCount = threads;
        return *this;
      }

      /// @returns Requested worker count; the default is one.
      size_t getThreadCount() const
      {
        return m_threadCount;
      }

      /// @returns Number of active-list rounds in the last solve.
      size_t getIterations() const
      {
        return m_iterations;
      }

      /**
       * @brief Reinitializes arrival times and the active list from zero-time seeds.
       * @tparam Seeds Iterable container of vertex indices.
       * @param seeds Nonempty valid source-vertex container; duplicates are harmless.
       * @returns This solver.
       */
      template <class Seeds>
      FIM& seed(const Seeds& seeds)
      {
        auto& arrival = m_arrival.get();
        const auto& space = arrival.getFiniteElementSpace();
        const auto& mesh = space.getMesh();
        const size_t dimension = mesh.getDimension();
        if (dimension > 0)
        {
          RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 0, dimension);
          RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, dimension, 0);
        }
        std::vector<bool> prescribed(mesh.getVertexCount(), false);
        bool any = false;
        for (Index vertex : seeds)
        {
          if (vertex >= prescribed.size())
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Seed vertex index out of range: " << vertex << Alert::Raise;
          }
          prescribed[vertex] = true;
          any = true;
        }
        if (!any)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "At least one seed vertex is required." << Alert::Raise;
        }
        m_prescribed = std::move(prescribed);
        arrival = std::numeric_limits<Real>::infinity();
        m_active.clear();
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          if (m_prescribed[vertex])
            arrival[space.getDOFs(0, vertex)[0]] = 0;
          else
            m_active.push_back(vertex);
        }
        m_iterations = 0;
        m_seeded = true;
        return *this;
      }

      /**
       * @brief Decreases active times until the local discrete updates converge.
       *
       * Each changed vertex reactivates all nonprescribed co-cell vertices.
       * Exhausting the round limit raises an exception, retaining the last iterate.
       * Unreachable connected components retain positive infinity.
       */
      void solve()
      {
        if (!m_seeded)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Call seed() before solve()." << Alert::Raise;
        }
        auto& arrival = m_arrival.get();
        const auto& space = arrival.getFiniteElementSpace();
        const auto& mesh = space.getMesh();
        const auto& connectivity = mesh.getConnectivity();
        const size_t dimension = mesh.getDimension();
        m_iterations = 0;
        if (!m_active.empty() && dimension > 0)
        {
          // Discretize the prescribed coefficient once per cell for this solve.
          m_cellCoefficients.clear();
          m_cellCoefficients.reserve(mesh.getCellCount());
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            if (!cell->getTransformation().isFirstOrder())
            {
              Alert::MemberFunctionException(*this, __func__)
                << "FIM requires native first-order cell geometry; curved/high-order "
                << "transformations are unsupported." << Alert::Raise;
            }
            const Geometry::Polytope::Traits traits(cell->getGeometry());
            const Geometry::Point center(*cell, traits.getCentroid());
            const auto value = m_coefficient(center);
            if constexpr (std::is_arithmetic_v<CoefficientValue>)
            {
              if (!(value > 0) || !std::isfinite(value))
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Speed must be finite and strictly positive." << Alert::Raise;
              }
              m_cellCoefficients.emplace_back(value);
            }
            else
              m_cellCoefficients.emplace_back(value, mesh.getSpaceDimension());
          }
        }
        while (!m_active.empty())
        {
          if (m_iterations == m_maxIterations)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "FIM did not converge within " << m_maxIterations
              << " active-list rounds." << Alert::Raise;
          }
          ++m_iterations;
          std::vector<uint8_t> changed(m_active.size(), 0);
          auto update = [&](Index position) {
            const Index vertex = m_active[position];
            const Real candidate = local(vertex);
            Real& time = arrival[space.getDOFs(0, vertex)[0]];
            if (std::isfinite(candidate) && time - candidate > m_tolerance)
            {
              time = candidate;
              changed[position] = 1;
            }
          };
          if (m_threadCount == 1 || dimension == 0)
          {
            for (Index position = 0; position < m_active.size(); ++position)
              update(position);
          }
#ifdef RODIN_USE_OPENMP
          else
          {
            // Preserve every co-cell dependency of the serial active-list order.
            // A level contains no pair of target vertices sharing a cell.
            IndexVector cellLevel(mesh.getCellCount(), 0);
            std::vector<IndexVector> levels;
            for (Index position = 0; position < m_active.size(); ++position)
            {
              const Index vertex = m_active[position];
              Index level = 0;
              for (Index cell : connectivity.getIncidence(0, dimension).at(vertex))
                level = std::max(level, cellLevel[cell]);
              if (levels.size() <= level)
                levels.resize(level + 1);
              levels[level].push_back(position);
              for (Index cell : connectivity.getIncidence(0, dimension).at(vertex))
                cellLevel[cell] = level + 1;
            }
            std::exception_ptr error;
            std::atomic<bool> failed(false);
#pragma omp parallel num_threads(m_threadCount)
            {
              auto capture = [&]() {
#pragma omp critical(RodinEikonalFIMException)
                {
                  if (!error)
                    error = std::current_exception();
                }
                failed.store(true, std::memory_order_relaxed);
              };
              for (const auto& level : levels)
              {
#pragma omp for schedule(static)
                for (Index entry = 0; entry < level.size(); ++entry)
                {
                  if (failed.load(std::memory_order_relaxed))
                    continue;
                  try
                  {
                    update(level[entry]);
                  }
                  catch (...)
                  {
                    capture();
                  }
                }
              }
            }
            if (error)
              std::rethrow_exception(error);
          }
#endif
          IndexVector next;
          std::vector<bool> queued(mesh.getVertexCount(), false);
          for (Index position = 0; position < m_active.size(); ++position)
          {
            if (!changed[position])
              continue;
            const Index vertex = m_active[position];
            for (Index cell : connectivity.getIncidence(0, dimension).at(vertex))
            {
              for (Index neighbor : connectivity.getIncidence(dimension, 0).at(cell))
              {
                if (!m_prescribed[neighbor] && !queued[neighbor])
                {
                  next.push_back(neighbor);
                  queued[neighbor] = true;
                }
              }
            }
          }
          m_active = std::move(next);
        }
      }

    private:
      /**
       * @brief Minimizes the native trace candidate at one vertex.
       * @param vertex Target global vertex index.
       * @returns Best finite candidate found, or infinity when no support is reached.
       */
      Real local(Index vertex) const
      {
        const auto& arrival = m_arrival.get();
        const auto& space = arrival.getFiniteElementSpace();
        const auto& mesh = space.getMesh();
        const auto& connectivity = mesh.getConnectivity();
        const size_t dimension = mesh.getDimension();
        if (dimension == 0)
          return std::numeric_limits<Real>::infinity();
        Real best = std::numeric_limits<Real>::infinity();
        std::vector<typename Geometry::Connectivity<Context::Local>::SubPolytope>
          entities;
        for (Index cellIndex : connectivity.getIncidence(0, dimension).at(vertex))
        {
          const auto cell = mesh.getCell(cellIndex);
          const auto& vertices = cell->getVertices();
          const Geometry::Polytope::Traits cellTraits(cell->getGeometry());
          auto cellVertex = [&](Index global) {
            return static_cast<size_t>(
              std::find(vertices.begin(), vertices.end(), global) - vertices.begin());
          };
          const Math::SpatialPoint origin =
            Geometry::Point(*cell, cellTraits.getVertex(cellVertex(vertex)))
              .getPhysicalCoordinates();
          const auto& coefficient = m_cellCoefficients[cellIndex];
          for (size_t entityDimension = 0; entityDimension < dimension; ++entityDimension)
          {
            connectivity.getSubPolytopes(entities, cellIndex, entityDimension);
            for (const auto& entity : entities)
            {
              if (std::find(entity.vertices.begin(), entity.vertices.end(), vertex) !=
                entity.vertices.end())
                continue;
              std::array<Real, 4> times;
              std::array<Math::SpatialPoint, 4> positions;
              bool reached = true;
              Real minimumTime = std::numeric_limits<Real>::infinity();
              for (size_t localVertex = 0; localVertex < entity.vertices.size();
                ++localVertex)
              {
                const Index global = entity.vertices[localVertex];
                const Real time = arrival[space.getDOFs(0, global)[0]];
                reached = reached && std::isfinite(time);
                times[localVertex] = time;
                minimumTime = std::min(minimumTime, time);
                positions[localVertex] =
                  Geometry::Point(*cell, cellTraits.getVertex(cellVertex(global)))
                    .getPhysicalCoordinates();
              }
              if (!reached || minimumTime >= best)
                continue;
              const Variational::P1Element<Real> element(entity.geometry);
              const Geometry::Polytope::Traits traits(entity.geometry);
              auto objective = [&](const Math::SpatialPoint& parameter) {
                Math::SpatialPoint entry(origin.size());
                entry.setZero();
                Real time = 0;
                for (size_t localVertex = 0; localVertex < entity.vertices.size();
                  ++localVertex)
                {
                  const Real weight = element.getBasis(localVertex)(parameter);
                  time += weight * times[localVertex];
                  entry += weight * positions[localVertex];
                }
                const Math::SpatialVector<Real> displacement = entry - origin;
                Real cost;
                if constexpr (std::is_arithmetic_v<CellCoefficient>)
                  cost = displacement.norm() / coefficient;
                else
                  cost = std::sqrt(coefficient.dot(displacement, displacement));
                if (!(cost > 0) || !std::isfinite(cost))
                {
                  Alert::MemberFunctionException(*this, __func__)
                    << "Invalid straight-cell geometry or nonfinite travel cost."
                    << Alert::Raise;
                }
                return time + cost;
              };
              if (entityDimension == 0)
              {
                best = std::min(best, objective(traits.getCentroid()));
                continue;
              }
              for (size_t start = 0; start <= traits.getVertexCount(); ++start)
              {
                Math::SpatialPoint parameter = start == traits.getVertexCount()
                  ? traits.getCentroid()
                  : traits.getVertex(start);
                Real value = objective(parameter);
                Real step = InitialSearchStep;
                size_t iteration = 0;
                while (step > m_localTolerance && iteration < MaxSearchIterations)
                {
                  ++iteration;
                  Math::SpatialPoint next = parameter;
                  Real nextValue = value;
                  for (Integer first = -1; first <= 1; ++first)
                  {
                    for (Integer second = entityDimension == 2 ? -1 : 0;
                      second <= (entityDimension == 2 ? 1 : 0); ++second)
                    {
                      if (first == 0 && second == 0)
                        continue;
                      Math::SpatialPoint trial = parameter;
                      trial(0) += step * first;
                      if (entityDimension == 2)
                        trial(1) += step * second;
                      Math::SpatialPoint projected;
                      Geometry::Polytope::Project(entity.geometry).cell(projected, trial);
                      if ((projected - parameter).squaredNorm() == 0)
                        continue;
                      const Real candidate = objective(projected);
                      if (candidate < nextValue)
                      {
                        nextValue = candidate;
                        next = projected;
                      }
                    }
                  }
                  if (nextValue < value)
                  {
                    parameter = next;
                    value = nextValue;
                  }
                  else
                    step /= 2;
                }
                if (step > m_localTolerance)
                {
                  Alert::MemberFunctionException(*this, __func__)
                    << "Native face search exhausted its iteration limit."
                    << Alert::Raise;
                }
                best = std::min(best, value);
              }
            }
          }
        }
        return best;
      }

      std::reference_wrapper<SolutionType> m_arrival;
      CoefficientType m_coefficient;
      std::vector<CellCoefficient> m_cellCoefficients;
      Real m_tolerance;
      Real m_localTolerance;
      size_t m_maxIterations;
      size_t m_iterations;
      size_t m_threadCount;
      bool m_seeded;
      std::vector<bool> m_prescribed;
      IndexVector m_active;
  };

  /**
   * @brief Deduces the native-cell solver and an owned coefficient type.
   * @tparam Solution Mutable output field type.
   * @tparam Coefficient Callable or Rodin function expression type.
   * @param arrival Mutable arrival-time field.
   * @param coefficient Scalar speed or SPD squared-velocity tensor.
   */
  template <class Solution, class Coefficient>
  FIM(Solution& arrival, const Coefficient& coefficient) -> FIM<Solution, Coefficient>;
}

#endif
