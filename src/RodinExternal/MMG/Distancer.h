/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_DISTANCER_H
#define RODIN_EXTERNAL_MMG_DISTANCER_H

#include <set>
#include <optional>

#include "Rodin/Alert.h"
#include "Rodin/Mesh.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "MMG5.h"
#include "Common.h"
#include "Configure.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"

namespace Rodin::External::MMG
{
  template <class T>
  class Distancer;

  template <class FEC>
  class Distancer<Variational::FiniteElementSpace<FEC>>
  {
    Distancer(Variational::FiniteElementSpace<FEC, Traits::Serial>& fes)
      : m_fes(fes),
        m_scale(true),
        m_ncpu(std::thread::hardware_concurrency()),
        m_mshdist(MSHDIST_EXECUTABLE)
    {}

    /**
     * @brief Specifies which material reference is to be understood as
     * the interior domains.
     * @param[in] ref Interior material reference
     * @returns Reference to self (for method chaining)
     *
     * By default, the interior elements are assumed to have the material
     * reference 3.
     *
     * @see distance(Mesh&) const
     */
    Distancer& setInteriorDomain(const MaterialReference& ref)
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
     * @see distance(Mesh&) const
     */
    Distancer& setInteriorDomain(const std::set<MaterialReference>& refs);

    /**
     * @brief Specifies whether to enable the scaling of the contour mesh.
     *
     * Specifies whether the contour mesh should be scaled down so that the
     * contour's bounding box is 95% of the enclosing bounding box.
     *
     * By default, it is enabled.
     *
     * @see distance(Mesh&, Mesh&)
     */
    Distancer& enableScaling(bool b = true);

    /**
     * @brief Specifies how many CPUs to use when distancing in parallel.
     * @param[in] ncpu Number of CPUs to use
     *
     * By default, it will utilize `std::thread::hardware_concurrency()`.
     */
    Distancer& setCPUs(unsigned int ncpu);

    /**
     * @returns Number of CPUs which will be used when performing the
     * distancing.
     */
    unsigned int getCPUs() const;

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
    Variational::GridFunction<FEC, Traits::Serial> distance(const Mesh<Traits::Serial>& box)
    {
      auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
      box.save(boxp, IO::MeshFormat::MEDIT);

      if (m_interiorDomains)
      {
        auto paramp(boxp);
        paramp.replace_extension(".mshdist");
        std::ofstream paramf(paramp.string(), std::ios::trunc);
        if (m_interiorDomains->size() > 0)
        {
          paramf << "InteriorDomains\n"
                 << m_interiorDomains->size() << "\n\n";
          for (const auto& ref : *m_interiorDomains)
            paramf << ref << "\n";
        }
      }

      int retcode = m_mshdist.run(
          boxp.string(),
          "-dom",
          "-ncpu", m_ncpu,
          "-v 0");

      Variational::GridFunction res(m_fes);
      if (retcode != 0)
        Alert::Exception("ISCD::Mshdist invocation failed.").raise();
      else
        res.load(boxp.replace_extension(".sol"));
      return res;
    }

    private:
      Variational::FiniteElementSpace<FEC, Traits::Serial>& m_fes;
      bool m_scale;
      unsigned int m_ncpu;
      ISCDProcess m_mshdist;
      std::optional<std::set<MaterialReference>> m_interiorDomains;
  };
  template <class FEC>
  Distancer(Variational::FiniteElementSpace<FEC>&)
    -> Distancer<Variational::FiniteElementSpace<FEC>>;

}

#endif
