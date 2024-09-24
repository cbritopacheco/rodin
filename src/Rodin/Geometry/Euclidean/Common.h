/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_EUCLIDEAN_COMMON_H
#define RODIN_GEOMETRY_EUCLIDEAN_COMMON_H

namespace Rodin::Geometry::Euclidean
{
  template <class P0, class ... Ps>
  auto barycenter(const P0& p0, const Ps&... ps)
  {
    return (p0 + ... + ps) / (sizeof...(ps) + 1);
  }
}

#endif
