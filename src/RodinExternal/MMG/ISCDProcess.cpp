/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <optional>
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
  boost::filesystem::path ISCDProcess::s_tmpDirPath = boost::filesystem::temp_directory_path();
  std::random_device ISCDProcess::s_randomDevice;
  std::mt19937 ISCDProcess::s_rng(ISCDProcess::s_randomDevice());
  std::uniform_int_distribution<std::mt19937::result_type> ISCDProcess::s_dist;

  ISCDProcess::ISCDProcess(const boost::filesystem::path& executable)
   : m_executable(executable)
  {}

  ISCDProcess::ISCDProcess(const ISCDProcess& other)
   : m_executable(other.m_executable),
    m_ncpu(other.m_ncpu)
  {}

  boost::filesystem::path ISCDProcess::tmpnam(
    const boost::filesystem::path& extension,
    const boost::filesystem::path& prefix)
  {
   boost::filesystem::path res;
   int counter = 0;
   do
   {
    res = s_tmpDirPath / (prefix.string() + n2hexstr(s_dist(s_rng)));
    res.replace_extension(extension);
    if (counter++ > 255)
      throw std::runtime_error(
        "Failed to obtain a temporary filename after 256 tries.");
   } while(boost::filesystem::exists(res));
   return res;
  }
}
