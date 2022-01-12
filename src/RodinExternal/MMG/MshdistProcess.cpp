/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>
#include "MshdistProcess.h"

namespace Rodin::External::MMG
{
  std::filesystem::path MshdistProcess::s_tmpDirPath = std::filesystem::temp_directory_path();
  std::random_device MshdistProcess::s_randomDevice;
  std::mt19937 MshdistProcess::s_rng(MshdistProcess::s_randomDevice());
  std::uniform_int_distribution<std::mt19937::result_type> MshdistProcess::s_dist;

  MshdistProcess::MshdistProcess()
    : m_ncpu(std::thread::hardware_concurrency())
  {}

  MshdistProcess MshdistProcess::setCPUs(unsigned int ncpu)
  {
    m_ncpu = ncpu;
    return *this;
  }

  unsigned int MshdistProcess::getCPUs() const
  {
    return m_ncpu;
  }

  int MshdistProcess::run(const std::vector<std::string>& args) const
  {
    // Accumulate args
    std::string strArgs;
    for (const auto& a : args)
      strArgs += " " + a;

    // Run command
    std::stringstream command;
    command << RODIN_MSHDIST_EXECUTABLE
            << strArgs
            << " -v " << RODIN_MMG_VERBOSITY_LEVEL
            << " -ncpu " << getCPUs();
    return std::system(command.str().c_str());
  }

  std::filesystem::path MshdistProcess::tmpnam(
      const std::filesystem::path& extension,
      const std::filesystem::path& prefix)
  {
    std::filesystem::path res;
    int counter = 0;
    do
    {
      res = s_tmpDirPath.string() + prefix.string() + n2hexstr(s_dist(s_rng));
      res.replace_extension(extension);
      if (counter++ > 255)
        throw std::runtime_error(
            "Failed to obtain a temporary filename after 256 tries.");
    } while(std::filesystem::exists(res));
    return res;
  }
}
