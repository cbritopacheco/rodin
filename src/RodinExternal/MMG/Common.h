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
    MaterialReference interior, /// Reference for the interior domain
                      exterior; /// Reference for exterior domain
  };

  /**
   * @brief Map indicating how a material reference should be split into
   * exterior and interior material references.
   */
  using SplitMap = std::map<
    MaterialReference, std::variant<Split, NoSplitT>>;
}

#endif
