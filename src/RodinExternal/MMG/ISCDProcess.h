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
#include <iostream>
#include <boost/filesystem.hpp>
#include <type_traits>

#include <boost/process.hpp>

#include "Rodin/Alert.h"

#include "Common.h"
#include "Configure.h"

namespace Rodin::External::MMG
{
  /**
   * @internal
   * @brief Represents a process of the ISCD Toolbox.
   *
   * This class provides utilities to deal with ISCD processes, such as getting
   * temporary filenames, running commands, passing command line arguments,
   * etc.
   */
  class ISCDProcess
  {
    public:
      template <class ... Ts>
      class ParameterPack
      {
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

          std::string m_str;
      };

      ISCDProcess(const boost::filesystem::path& executable);

      ISCDProcess(const ISCDProcess& other);

      ISCDProcess(ISCDProcess&& other) = default;

      template <class ... Ts>
      int run(Ts... args)
      {
        // Accumulate args
        ParameterPack ps(std::forward<Ts>(args)...);
        std::string strArgs = ps.toString();

        // Run command
        if (m_out)
          return boost::process::system(
              m_executable.string() + " " + strArgs,
              boost::process::std_out > *m_out
            );
        else
          return boost::process::system(
              m_executable.string() + " " + strArgs,
              boost::process::std_out > boost::process::null
            );
      }

      int run(const std::vector<std::string>& args)
      {
        // Run command
        std::string strArgs = "";
        for (const auto& arg : args)
          strArgs += " " + arg;
        if (m_out)
          return boost::process::system(
              m_executable.string() + strArgs,
              boost::process::std_out > *m_out,
              boost::process::std_err > *m_out
            );
        else
          return boost::process::system(
              m_executable.string() + strArgs,
              boost::process::std_out > boost::process::null,
              boost::process::std_out > boost::process::null
            );
      }

      ISCDProcess& logOutput(bool enable = true)
      {
        if (enable)
        {
          if (!m_out)
            m_out.emplace();
        }
        else
        {
          m_out.reset();
        }
        return *this;
      }

      std::optional<boost::process::ipstream>& getOutputLog()
      {
        return m_out;
      }

      void printOutputLog(std::ostream& os = std::cout)
      {
        if (m_out)
        {
          std::string line;
          while (*m_out && std::getline(*m_out, line))
            os << line << '\n';
        }
      }

      static boost::filesystem::path tmpnam(
          const boost::filesystem::path& extension,
          const boost::filesystem::path& prefix = boost::filesystem::path());

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

      static boost::filesystem::path s_tmpDirPath;
      static std::random_device s_randomDevice;
      static std::mt19937 s_rng;
      static std::uniform_int_distribution<std::mt19937::result_type> s_dist;

      boost::filesystem::path m_executable;
      unsigned int m_ncpu;
      std::optional<boost::process::ipstream> m_out;
  };
}

#endif


