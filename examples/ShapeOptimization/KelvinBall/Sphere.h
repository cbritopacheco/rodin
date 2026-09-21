/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_SPHERE_H
#define KELVIN_BALL_SPHERE_H

#include <Rodin/MMG.h>

#include "Configuration.h"

namespace KelvinBall
{
  struct ReconstructionDiagnostics
  {
      Real minimumSize = 0;
      Real maximumSize = 0;
      Real hausdorffTolerance = 0;
      size_t requiredBoundaryTriangles = 0;
      size_t cellsBefore = 0;
      size_t cellsAfter = 0;
  };

  struct SphereDiscretization
  {
      MMG::Mesh mesh;
      ReconstructionDiagnostics diagnostics;
  };

  /**
   * @brief Constructs and discretizes the spherical Kelvin chamber.
   *
   * The construction owns the background chamber, the protection of its fixed
   * boundary, and the MMG level-set reconstruction.  Keeping these stages in
   * one object makes the geometric invariants explicit and prevents the
   * construction from being fragmented into unrelated free functions.
   */
  class Sphere
  {
    public:
      explicit Sphere(const Configuration& configuration);

      SphereDiscretization discretize() const;

      /**
       * @brief Prepares the fixed background mesh used by WNGIR.
       *
       * MMG acts only on the interface-free chamber background. The returned
       * mesh contains no design interface and remains fixed during optimization.
       */
      SphereDiscretization prepareWNGIRBackground() const;

      /**
       * @brief Adapts a mesh fitted to Gamma to the prescribed size map.
       *
       * Replaces the MMG optimization pass after a level-set cut when
       * `--mmg-adapt` is given. The optimization pass keeps the edge lengths the
       * cut created, including its short edges; adaptation remeshes towards
       * a size set by the distance to Gamma instead.
       */
      void adapt(MMG::Mesh& mesh) const;

    private:
      Mesh makeUniformChamber() const;

      size_t protectFixedGeometry(MMG::Mesh& mesh) const;

      const Configuration& m_configuration;
  };
}

#endif
