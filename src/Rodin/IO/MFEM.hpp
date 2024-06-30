/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MFEM_HPP
#define RODIN_IO_MFEM_HPP

#include "MFEM.h"

namespace Rodin::IO
{
  template <class Range>
  void GridFunctionLoader<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context::Sequential>>>
  ::load(std::istream& is)
  {
    using boost::spirit::x3::space;
    using boost::spirit::x3::blank;
    using boost::spirit::x3::uint_;
    using boost::spirit::x3::_attr;
    using boost::spirit::x3::char_;

    MFEM::GridFunctionHeader header;
    const auto get_fec = [&](auto& ctx) { header.fec = _attr(ctx); };
    const auto get_vdim = [&](auto& ctx) { header.vdim = _attr(ctx); };
    const auto get_ordering = [&](auto& ctx) { header.ordering = static_cast<MFEM::Ordering>(_attr(ctx)); };

    std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    auto it = line.begin();
    const auto pfes = boost::spirit::x3::string("FiniteElementSpace");
    const bool rfes = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
    assert(it == line.end() && rfes);

    line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    it = line.begin();
    const auto pfec = boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
    bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
    assert(it == line.end() && rfec);

    line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    it = line.begin();
    const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
    bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
    assert(it == line.end() && rvdim);

    line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
    it = line.begin();
    const auto pordering = boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
    bool rordering = boost::spirit::x3::phrase_parse(it, line.end(), pordering, space);
    assert(it == line.end() && rordering);

    auto& gf = this->getObject();
    const auto& fes = gf.getFiniteElementSpace();
    assert(header.vdim == fes.getVectorDimension());
    auto& data = gf.getData();
    if (data.size() > 0)
    {
      line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
      data.coeffRef(0) = std::stod(line);
      assert(data.size() >= 0);
      for (size_t i = 1; i < static_cast<size_t>(data.size()); i++)
        is >> data.coeffRef(i);
      if (header.ordering == MFEM::Ordering::Nodes)
        data.transposeInPlace();
    }
    gf.setWeights();
  }
}

#endif

