/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MEDIT_HPP
#define RODIN_IO_MEDIT_HPP

#include <boost/algorithm/string.hpp>

#include "MEDIT.h"

namespace Rodin::IO
{
  template <class Range>
  void GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::load(std::istream& is)
  {
    readVersion(is);
    readDimension(is);
    readData(is);
  }

  template <class Range>
  std::istream& GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::getline(std::istream& is, std::string& line)
  {
    m_currentLineNumber++;
    return std::getline(is, line);
  }

  template <class Range>
  std::string GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::skipEmptyLines(std::istream& is)
  {
    std::string line;
    while (getline(is, line))
    {
      if (!MEDIT::ParseEmptyLine()(line.begin(), line.end()))
        break;
    }
    return line;
  }

  template <class Range>
  void GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::readVersion(std::istream& is)
  {
    auto line = skipEmptyLines(is);
    std::optional<unsigned int> version = MEDIT::ParseMeshVersionFormatted()(line.begin(), line.end());
    if (version) // Version was on the same line
    {
      m_version = *version;
    }
    else // Version is not on the same line
    {
      auto line = skipEmptyLines(is);
      version = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
      if (version)
        m_version = *version;
      else
        throw Alert::Exception("Failed to parse version number of mesh.");
    }
  }

  template <class Range>
  void GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::readDimension(std::istream& is)
  {
    auto line = skipEmptyLines(is);
    std::optional<unsigned int> dimension = MEDIT::ParseDimension()(line.begin(), line.end());
    if (dimension) // Version was on the same line
      m_spaceDimension = *dimension;
    else // Version is not on the same line
    {
      auto line = skipEmptyLines(is);
      dimension = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
      if (dimension)
        m_spaceDimension = *dimension;
      else
        throw Alert::Exception("Failed to parse dimension of mesh.");
    }
  }

  template <class Range>
  void GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Context::Serial, Geometry::Mesh<Context::Serial>>>
  ::readData(std::istream& is)
  {
    auto& gf = this->getObject();

    auto line = skipEmptyLines(is);
    std::optional<std::string> kw = MEDIT::ParseKeyword()(line.begin(), line.end());
    if (!kw || *kw != MEDIT::Keyword::SolAtVertices)
    {
      Alert::Exception() << "Expected keyword " << MEDIT::Keyword::SolAtVertices
                         << " on line " << m_currentLineNumber 
                         << Alert::Raise;
    }

    line = skipEmptyLines(is);
    std::optional<unsigned int> size = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
    if (!size)
    {
      Alert::Exception() << "Failed to parse solution size at line "
                         << m_currentLineNumber
                         << Alert::Raise;
    }

    line = skipEmptyLines(is);
    size_t solCount, vdim;
    using boost::spirit::x3::space;
    using boost::spirit::x3::blank;
    using boost::spirit::x3::uint_;
    using boost::spirit::x3::_attr;
    using boost::spirit::x3::repeat;
    const auto get_sol_count = [&](auto& ctx) { solCount = _attr(ctx); };
    const auto get_vdim = [&](auto& ctx) { vdim = _attr(ctx); };
    const auto p = uint_[get_sol_count] >> uint_[get_vdim];
    auto it = line.begin();
    const bool r = boost::spirit::x3::phrase_parse(it, line.end(), p, space);

    assert(solCount == 1);
    if (it != line.end() || !r)
    {
      Alert::Exception() << "Failed to parse solution count and vector dimension at line "
                         << m_currentLineNumber
                         << Alert::Raise;
    }

    auto& data = gf.getData();
    assert(data.rows() >= 0);
    assert(static_cast<size_t>(data.rows()) == vdim);
    assert(data.cols() % vdim == 0);
    assert(data.cols() / vdim == size);
    assert(data.size() >= 0);
    for (size_t i = 0; i < static_cast<size_t>(data.size()); i++)
      is >> data.coeffRef(i);
  }
}

#endif
