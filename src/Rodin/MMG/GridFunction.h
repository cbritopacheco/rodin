/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GridFunction.h
 * @brief Grid function aliases for MMG-driven workflows.
 */
#ifndef RODIN_RODINEXTERNAL_MMG_GRIDFUNCTION_H
#define RODIN_RODINEXTERNAL_MMG_GRIDFUNCTION_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/P1.h"

namespace Rodin::MMG
{
  /**
   * @brief Grid function type used by the MMG module.
   * @tparam Range Value type (`Real` for scalar fields, `Math::Vector<Real>`
   * for vector fields).
   *
   * MMG workflows in Rodin operate on first-order nodal fields defined on
   * @ref Rodin::Geometry::Mesh<Context::Local> "local meshes". This alias
   * standardizes that representation for scalar and vector quantities transferred
   * to/from MMG solution structures.
   */
  template <class Range>
  using GridFunction = Variational::GridFunction<
        Variational::P1<Range, Geometry::Mesh<Context::Local>>, Math::Vector<Real>>;

  /**
   * @brief Scalar MMG grid function alias.
   *
   * Commonly used for size maps and level-set values passed to
   * @ref Adapt and @ref LevelSetDiscretizer.
   */
  using RealGridFunction = GridFunction<Real>;

  /**
   * @brief Vector-valued MMG grid function alias.
   *
   * Useful for vector metrics or vector fields in MMG-related pre/post-processing.
   */
  using VectorGridFunction = GridFunction<Math::Vector<Real>>;
}

#endif
