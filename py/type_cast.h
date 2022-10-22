/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <boost/filesystem/path.hpp>

namespace pybind11
{
  namespace detail
  {
    template <>
    struct type_caster<boost::filesystem::path, path_caster<boost::filesystem::path>>
    {};
  }
}
