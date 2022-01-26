/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_ISCDPROCESS_H
#define RODIN_EXTERNAL_MMG_ISCDPROCESS_H

#include <tuple>
#include <random>
#include <vector>
#include <sstream>
#include <filesystem>
#include <type_traits>

#include "Configure.h"

namespace Rodin::External::MMG
{
  class ISCDProcess
  {
    public:
      template <class ... Ts>
      class ParameterPack
      {
        template<std::size_t I = 0, typename... Tp>
        inline typename std::enable_if<I == sizeof...(Tp), void>::type
        concatTuple(const std::tuple<Tp...>&)
        {}

        template<std::size_t I = 0, typename... Tp>
        inline typename std::enable_if<I < sizeof...(Tp), void>::type
        concatTuple(const std::tuple<Tp...>& t)
        {
          if constexpr (
              std::is_convertible_v<decltype(std::get<I>(t)), std::string>)
          {
            m_str += " " + static_cast<std::string>(std::get<I>(t));
          }
          else
          {
            m_str += " " + std::to_string(std::get<I>(t));
          }
          concatTuple<I + 1, Tp...>(t);
        }

        public:
          ParameterPack(Ts... ts)
            : m_str("")
          {
            concatTuple(std::forward_as_tuple(ts...));
          }

          std::string toString() const
          {
            return m_str;
          }

        private:
          std::string m_str;
      };

      ISCDProcess(const std::filesystem::path& executable);

      template <class ... Ts>
      int run(Ts... args) const
      {
        // Accumulate args
        ParameterPack ps(std::forward<Ts>(args)...);
        std::string strArgs = ps.toString();

        // Run command
        return std::system((m_executable.string() + " " + strArgs).c_str());
      }

      int run(const std::vector<std::string>& args) const
      {
        // Run command
        std::string strArgs = "";
        for (const auto& arg : args)
          strArgs += " " + arg;
        return std::system((m_executable.string() + strArgs).c_str());
      }

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

      std::filesystem::path m_executable;
      unsigned int m_ncpu;
  };
}

#endif


