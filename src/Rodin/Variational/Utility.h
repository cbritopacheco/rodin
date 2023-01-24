#ifndef RODIN_VARIATIONAL_UTILITY_H
#define RODIN_VARIATIONAL_UTILITY_H

#include <mfem.hpp>

namespace Rodin::Variational
{
  mfem::ElementTransformation* refinedToCoarse(
    mfem::Mesh &coarse_mesh, const mfem::ElementTransformation &T,
    const mfem::IntegrationPoint &ip, mfem::IntegrationPoint &coarse_ip);
}

#endif
