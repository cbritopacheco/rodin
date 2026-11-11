/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_IO_GRIDFUNCTIONLOADER_H
#define RODIN_KOKKOS_IO_GRIDFUNCTIONLOADER_H
/**
 * @file
 * @brief File-format adapters for backend grid-function coefficients.
 */
#include "Rodin/IO/GridFunctionLoader.h"
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
  class GridFunctionLoader<Format, FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>
    final
    : public Loader<
        Variational::GridFunction<FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>>
  {
    public:
      /// @brief Backend field type.
      using ObjectType =
        Variational::GridFunction<FES, ::Rodin::Kokkos::Math::Vector<Scalar, Backend>>;

      /// @brief Temporary host coefficient field.
      using HostField = Variational::GridFunction<FES, Math::Vector<Scalar>>;

      /// @param gf Backend field, which must outlive this adapter.
      explicit GridFunctionLoader(ObjectType& gf)
        : m_gf(gf)
      {}

      /// @param target Format stream.
      void load(std::istream& target) override
      {
        HostField host(m_gf.getFiniteElementSpace());
        GridFunctionLoader<Format, FES, Math::Vector<Scalar>>(host).load(target);
        m_gf.setData(host.getData());
        if (const auto name = host.getName())
          m_gf.setName(std::string(*name));
      }

      /// @param target Format file path.
      void load(const boost::filesystem::path& target) override
      {
        HostField host(m_gf.getFiniteElementSpace());
        GridFunctionLoader<Format, FES, Math::Vector<Scalar>>(host).load(target);
        m_gf.setData(host.getData());
        if (const auto name = host.getName())
          m_gf.setName(std::string(*name));
      }

      /// @return Backend field owned by the caller.
      ObjectType& getObject() override
      {
        return m_gf;
      }

    private:
      ObjectType& m_gf;
  };
}
#endif
