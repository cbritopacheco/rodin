/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Common.h
 * @brief Common MMG integration types and utility accessors.
 */
#ifndef RODIN_EXTERNAL_MMG_COMMON_H
#define RODIN_EXTERNAL_MMG_COMMON_H

#include <variant>

#include <mmg/libmmg.h>
#include <libmmgcommon_private.h>
#include <mmg2d/libmmg2d_private.h>
#include <mmg/mmg2d/libmmg2d.h>
#include <mmg3d/libmmg3d_private.h>
#include <mmg/mmg3d/libmmg3d.h>
#include <mmgs/libmmgs_private.h>
#include <mmg/mmgs/libmmgs.h>
#include <common/mmgcommon_private.h>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Types.h"

#include "Configure.h"

/*
 * mmg includes complex.h which defines the I macro. We have to undefine it to
 * avoid a whole bunch of name clashes!
 */
#ifdef I
#undef I
#endif

namespace Rodin::MMG
{
  /**
   * @brief Tag type indicating that a material reference must not be split.
   *
   * Used by @ref SplitMap for level-set multi-material workflows in
   * @ref LevelSetDiscretizer.
   */
  struct NoSplitT {};

  /**
   * @brief Singleton tag instance for convenience in split configuration.
   *
   * Example:
   * @code
   * SplitMap split;
   * split[5] = NoSplit;
   * @endcode
   */
  static constexpr NoSplitT NoSplit;

  /**
   * @brief Interior/exterior labels produced when splitting one material.
   *
   * Used as the mapped value in @ref SplitMap.
   */
  struct Split
  {
    Geometry::Attribute   interior; ///< Target attribute for the interior side.
    Geometry::Attribute   exterior; ///< Target attribute for the exterior side.
  };

  /**
   * @brief Material splitting policy map used by level-set discretization.
   *
   * Keys are input material attributes from the source mesh.
   * Values describe whether each material is split into distinct interior/exterior
   * labels (@ref Split) or kept unchanged (@ref NoSplitT).
   *
   * This map is consumed by @ref LevelSetDiscretizer::setSplit.
   */
  using SplitMap = UnorderedMap<Geometry::Attribute, std::variant<Split, NoSplitT>>;

  /**
   * @brief Gets the configured path to the `mshdist` executable.
   * @returns Null-terminated executable path.
   */
  const char* getISCDMshdistExecutable();

  /**
   * @brief Gets the configured path to the `advect` executable.
   * @returns Null-terminated executable path.
   */
  const char* getISCDAdvectExecutable();

  /**
   * @brief Gets the MMG verbosity level selected at configuration time.
   * @returns MMG verbosity value passed to MMG `imprim` settings.
   */
  int getMMGVerbosityLevel();
}

#endif
