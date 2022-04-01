/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_DISTANCERS_H
#define RODIN_EXTERNAL_MMG_DISTANCERS_H

#include <set>
#include <optional>

#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
  class DistancerS
  {
    public:
      /**
       * @brief Creates a DistancerS object with default values.
       */
      DistancerS();

      /**
       * @brief Computes a signed distance function to a subdomain.
       *
       * This function generates the signed distance function @f$ d_\Omega @f$
       * to a domain @f$ \Omega @f$ which is contained in the given bounding
       * box @f$ D @f$, where @f$ \Omega @f$ is identified by its material
       * reference.
       * By default, the elements of @f$ \Omega @f$ are assumed to be specified
       * by the number 3.
       *
       * @param[in] box Bounding box @f$ D @f$ containing @f$ \Omega @f$.
       * @returns Signed distance function representing @f$ \Omega @f$.
       */
      ScalarSolutionS distance(MeshS& box);

      /**
       * @brief Computes a signed distance function from a given bounding box
       * and contour.
       *
       * This function generates the signed distance function @f$ d_\Omega @f$
       * to a domain @f$ \Omega @f$ supplied by means of a mesh of its boundary
       * @f$ \partial \Omega @f$, at the vertices of a computational mesh of a
       * bounding box @f$ D @f$.
       *
       * @param[in] box Bounding box @f$ D @f$ containing the contour.
       * @param[in] contour Orientable contour @f$ \partial \Omega @f$ to distance.
       * @returns Signed distance function representing @f$ \Omega @f$.
       * @note The contour mesh is allowed to contain a volume part, in which
       * case only the edge (S) or triangle (3D) information will be retained.
       */
      ScalarSolutionS distance(MeshS& box, MeshS& contour);

      /**
       * @brief Redistances the level set function.
       *
       * Given a level set function defined on the vertices of a bounding box
       * @f$ D @f$, this method will regenerate the signed distance function
       * associated to the the subdomain @f$ \Omega \subset D @f$ which the
       * level set function represents.
       *
       * @param[in,out] sol Level set function
       */
      void redistance(ScalarSolutionS& sol);

      /**
       * @brief Specifies which material reference is to be understood as
       * the interior domains.
       * @param[in] ref Interior material reference
       * @returns Reference to self (for method chaining)
       *
       * By default, the interior elements are assumed to have the material
       * reference 3.
       *
       * @see distance(MeshS&) const
       */
      DistancerS& setInteriorDomain(const MaterialReference& ref)
      {
        return setInteriorDomain(std::set<MaterialReference>{ref});
      }

      /**
       * @brief Specifies which material references are to be understood as
       * the interior domains.
       * @param[in] refs Interior material references
       * @returns Reference to self (for method chaining)
       *
       * By default, the interior elements are assumed to have the material
       * reference 3.
       *
       * @see distance(MeshS&) const
       */
      DistancerS& setInteriorDomain(const std::set<MaterialReference>& refs);

      /**
       * @brief Specifies whether to enable the scaling of the contour mesh.
       *
       * Specifies whether the contour mesh should be scaled down so that the
       * contour's bounding box is 95% of the enclosing bounding box.
       *
       * By default, it is enabled.
       *
       * @see distance(MeshS&, MeshS&)
       */
      DistancerS& enableScaling(bool b = true);

      /**
       * @brief Specifies how many CPUs to use when distancing in parallel.
       * @param[in] ncpu Number of CPUs to use
       *
       * By default, it will utilize `std::thread::hardware_concurrency()`.
       */
      DistancerS& setCPUs(unsigned int ncpu);

      /**
       * @returns Number of CPUs which will be used when performing the
       * distancing.
       */
      unsigned int getCPUs() const;

    private:
      bool m_scale;
      unsigned int m_ncpu;
      ISCDProcess m_mshdist;
      std::optional<std::set<MaterialReference>> m_interiorDomains;
      std::optional<std::set<MaterialReference>> m_activeBorders;
  };
}

#endif
