/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_COMMON_H
#define RODIN_EXTERNAL_MMG_COMMON_H

#include <map>
#include <variant>

#include <mmg/libmmg.h>
#include <libmmgcommon.h>
#include <mmg2d/mmg2d.h>
//#include <mmg2d/libmmg2d_private.h>
#include <mmg/mmg2d/libmmg2d.h>
#include <mmg3d/mmg3d.h>
//#include <mmg3d/libmmg3d_private.h>
#include <mmg/mmg3d/libmmg3d.h>
#include <mmgs/mmgs.h>
//#include <mmgs/libmmgs_private.h>
#include <mmg/mmgs/libmmgs.h>
#include <common/mmgcommon.h>

#include "Configure.h"

/*
 * mmg includes complex.h which defines the I macro. We have to undefine it to
 * avoid a whole bunch of name clashes!
 */
#ifdef I
#undef I
#endif

namespace Rodin::External::MMG
{
  /**
   * @brief Material reference for each element.
   */
  using MaterialReference = int;

  /**
   * @brief Empty class tag to specify that a material reference should not be
   * splitted.
   */
  struct NoSplitT {};

  /**
   * @brief Instance of the empty struct tag NoSplitT.
   */
  static constexpr NoSplitT NoSplit;

  /**
   * @brief Class to specify the interior and exterior material references.
   */
  struct Split
  {
    MaterialReference     interior; /// Reference for the interior domain
     MaterialReference    exterior; /// Reference for exterior domain
  };

  /**
   * @brief Map indicating how a material reference should be split into
   * exterior and interior material references.
   */
  using SplitMap = std::map<
    MaterialReference, std::variant<Split, NoSplitT>>;

  const char* getISCDMshdistExecutable();

  const char* getISCDAdvectExecutable();

  int getMMGVerbosityLevel();
}

#endif
