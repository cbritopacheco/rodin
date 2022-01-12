/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MSHDISTPROCESS_H
#define RODIN_EXTERNAL_MMG_MSHDISTPROCESS_H

#include <random>
#include <vector>
#include <sstream>
#include <filesystem>

#include "Configure.h"

namespace Rodin::External::MMG
{
  class MshdistProcess
  {
    public:
      MshdistProcess();

      MshdistProcess setCPUs(unsigned int ncpu);

      unsigned int getCPUs() const;

      int run(const std::vector<std::string>& args) const;

      static std::filesystem::path tmpnam(
          const std::filesystem::path& extension,
          const std::filesystem::path& prefix = std::filesystem::path());

    private:
      template <typename I>
      static std::string n2hexstr(I w, size_t hex_len = sizeof(I) << 1)
      {
          static const char* digits = "0123456789ABCDEF";
          std::string rc(hex_len, '0');
          for (size_t i=0, j=(hex_len - 1) * 4; i < hex_len; ++i, j -= 4)
              rc[i] = digits[(w >> j) & 0x0f];
          return rc;
      }

      static std::filesystem::path s_tmpDirPath;
      static std::random_device s_randomDevice;
      static std::mt19937 s_rng;
      static std::uniform_int_distribution<std::mt19937::result_type> s_dist;

      unsigned int m_ncpu;
  };
}

#endif
