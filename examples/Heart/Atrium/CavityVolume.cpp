/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "CavityVolume.h"

namespace Rodin::Examples::Heart::Atrium
{
  CavityVolume::CavityVolume(const MeshType& parent, Geometry::Attribute fluidMaterial)
    : m_parent(parent),
      m_fluidMaterial(fluidMaterial),
      m_referenceVolume(measure(parent, fluidMaterial))
  {
    if (m_referenceVolume <= 0 || !std::isfinite(m_referenceVolume))
    {
      Alert::Exception()
        << "Atrial volume requires a nonempty fluid region in the common mesh"
        << Alert::Raise;
    }
  }

  Real CavityVolume::getReferenceVolume() const noexcept
  {
    return m_referenceVolume;
  }

  Real CavityVolume::measure(const MeshType& mesh, Optional<Geometry::Attribute> material)
  {
    if (mesh.getDimension() != 3 || mesh.getSpaceDimension() != 3)
    {
      Alert::Exception() << "Cavity volume requires a three-dimensional volume mesh"
                         << Alert::Raise;
    }
    Real volume = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (mesh.getShard().isOwned(3, cell->getIndex()) &&
        (!material || cell->getAttribute() == *material))
        volume += cell->getMeasure(QuadratureOrder);
    }
    return boost::mpi::all_reduce(
      mesh.getContext().getCommunicator(), volume, std::plus<Real>());
  }
}
