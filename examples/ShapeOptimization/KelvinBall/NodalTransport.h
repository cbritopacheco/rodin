/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_NODAL_TRANSPORT_H
#define KELVIN_BALL_NODAL_TRANSPORT_H

#include <Rodin/Location/AABB.h>
#include <Rodin/Variational.h>

namespace KelvinBall
{
  /**
   * @brief P1 nodal characteristic transport without a variational projection.
   * Failed characteristics retain the previous nodal value. Each vertex is
   * evaluated once, independently of the ordering of its incident cells.
   * The caller owns the locator and must rebuild it when the mesh changes.
   */
  class NodalTransport
  {
    public:
      /**
       * @param mesh Chamber mesh shared by both scalar fields.
       * @param locator Shared locator on the unchanged chamber mesh.
       * @param distance Previous scalar field, distinct from result.
       * @param flow Characteristic tracer including the boundary policy.
       * @param[out] result Transported P1 scalar field.
       */
      template <class Mesh, class Distance, class Flow, class Result>
      NodalTransport(const Mesh& mesh, const Rodin::Location::AABB<Mesh>& locator,
        const Distance& distance, const Flow& flow, Result& result)
      {
        for (auto vertex = mesh.getVertex(); vertex; ++vertex)
        {
          const auto point =
            locator.locate(mesh.getDimension(), vertex->getCoordinates());
          if (!point)
            throw std::runtime_error("A transport vertex is outside its chamber mesh.");
          ++m_attempted;
          const auto trace = flow.trace(*point);
          const auto value = trace.exited()
            ? distance.getValue(*point)
            : distance.getValue(trace.getPoint()) + trace.getCorrection();
          if (trace.exited())
            ++m_fallback;
          result.project({0, vertex->getIndex()},
            Rodin::Variational::RealFunction(
              [value](const Rodin::Geometry::Point&) { return value; }));
        }
      }

      /** @returns Number of traced vertices. */
      size_t getAttemptedCount() const
      {
        return m_attempted;
      }

      /** @returns Number of vertices retaining their previous distance. */
      size_t getFallbackCount() const
      {
        return m_fallback;
      }

    private:
      size_t m_attempted = 0;
      size_t m_fallback = 0;
  };
}

#endif
