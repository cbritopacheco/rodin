/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the @ref Rodin::MMG module.
 */
#ifndef RODIN_RODININTEGRATION_MMG_FORWARDDECLS_H
#define RODIN_RODININTEGRATION_MMG_FORWARDDECLS_H

namespace Rodin::MMG
{
  /**
   * @brief Base wrapper for MMG configuration and low-level mesh/solution conversion.
   */
  class MMG5;

  /**
   * @brief MMG-aware mesh type extending @ref Rodin::Geometry::Mesh.
   */
  class Mesh;

  /**
   * @brief Metric-driven remeshing operator.
   */
  class Adapt;

  /**
   * @brief Mesh-quality optimization operator.
   */
  class Optimizer;

  /**
   * @brief Level-set based implicit-domain discretization operator.
   */
  class LevelSetDiscretizer;
}

#endif
