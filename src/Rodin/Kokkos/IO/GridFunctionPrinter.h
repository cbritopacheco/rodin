/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_IO_GRIDFUNCTIONPRINTER_H
#define RODIN_KOKKOS_IO_GRIDFUNCTIONPRINTER_H
/**
 * @file
 * @brief File-format adapters for backend grid-function coefficients.
 */
#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/Kokkos/Variational/GridFunction.h"

namespace Rodin::IO
{
  /**
   * @brief Stages coefficient data at the existing file-format boundary.
   * @tparam Format Existing Rodin file format. @tparam FES Finite element space.
   * @tparam Scalar Coefficient type. @tparam Backend Kokkos execution backend.
   *
   * ## Architecture
   * The format implementation operates on a temporary host coefficient vector.
   * The owning field remains backend-resident, and no host mirror is retained.
   */
  template <FileFormat Format, class FES, class Scalar, class Backend>
  class GridFunctionPrinter<Format, FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
    final
    : public Printer<
        Variational::GridFunction<FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>>
  {
    public:
      /// @brief Backend field type.
      using ObjectType =
        Variational::GridFunction<FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>;

      /// @brief Temporary host coefficient field.
      using HostField = Variational::GridFunction<FES, Math::Vector<Scalar>>;

      /// @param gf Backend field, which must outlive this adapter.
      explicit GridFunctionPrinter(const ObjectType& gf)
        : m_gf(gf)
      {}

      /// @param target Format stream.
      void print(std::ostream& target) override
      {
        HostField host(m_gf.getFiniteElementSpace());
        host.setData(m_gf.getData().toHost());
        if (const auto name = m_gf.getName())
          host.setName(std::string(*name));
        GridFunctionPrinter<Format, FES, Math::Vector<Scalar>>(host).print(target);
      }

      /// @param target Format file path.
      void print(const boost::filesystem::path& target) override
      {
        HostField host(m_gf.getFiniteElementSpace());
        host.setData(m_gf.getData().toHost());
        if (const auto name = m_gf.getName())
          host.setName(std::string(*name));
        GridFunctionPrinter<Format, FES, Math::Vector<Scalar>>(host).print(target);
      }

      /// @return Backend field owned by the caller.
      const ObjectType& getObject() const override
      {
        return m_gf;
      }

    private:
      const ObjectType& m_gf;
  };
}
#endif
