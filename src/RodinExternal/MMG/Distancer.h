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
#include "Rodin/Geometry.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "MMG5.h"
#include "Common.h"
#include "ForwardDecls.h"
#include "ISCDProcess.h"
#include "GridFunction.h"

namespace Rodin::External::MMG
{
  /**
  * @brief Class which generates the signed distance function of a domain.
  *
  * Given a domain @f$ \Omega \subset D @f$, with boundary @f$ \partial \Omega
  * @f$ then the signed distance function @f$ \phi : D \rightarrow \mathbb{R}
  * @f$ is defined by:
  *
  * @f[
  *  \phi(x) = \begin{cases}
  *   d(x, \partial \Omega) & \text{if } x \in \Omega\\
  *   0 & \text{if } x \in \partial \Omega\\
  *   -d(x, \partial \Omega) & \text{if } x \in \Omega^c
  *  \end{cases}
  * @f]
  * where
  * @f[
  *  d(x, \partial \Omega) := \inf_{y \in \partial \Omega} d(x, y)
  * @f]
  */
  template <>
  class Distancer<ScalarGridFunction::FES>
  {
    public:
     /**
      * @brief Creates a Distancer2D object with default values.
      * @param[in] fes Finite element space for the distance function
      */
     Distancer(const ScalarGridFunction::FES& fes)
       : m_fes(fes),
         m_scale(true),
         m_distTheBoundary(false),
         m_ncpu(std::thread::hardware_concurrency()),
         m_mshdist(getISCDMshdistExecutable())
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
     Distancer& setInteriorDomain(const Geometry::Attribute& ref)
     {
       return setInteriorDomain(FlatSet<Geometry::Attribute>{ref});
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
     Distancer& setInteriorDomain(const FlatSet<Geometry::Attribute>& refs)
     {
       m_interiorDomains = refs;
       return *this;
     }

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
     Distancer& enableScaling(bool b = true)
     {
       m_scale = b;
       return *this;
     }

     /**
      * @brief Specifies how many CPUs to use when distancing in parallel.
      * @param[in] ncpu Number of CPUs to use
      *
      * By default, it will utilize `std::thread::hardware_concurrency()`.
      */
     Distancer& setCPUs(unsigned int ncpu)
     {
       m_ncpu = ncpu;
       return *this;
     }

     /**
      * @returns Number of CPUs which will be used when performing the
      * distancing.
      */
     unsigned int getCPUs() const
     {
       return m_ncpu;
     }

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
     ScalarGridFunction distance(const Geometry::Mesh<Context::Serial>& box)
     {
       if (box != m_fes.get().getMesh())
       {
         Alert::Exception()
           << "Mesh must be the same one as that of the finite element space."
           << Alert::Raise;
       }

       auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
       box.save(boxp, IO::FileFormat::MEDIT);

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

       int retcode = 1;
       if (box.isSurface())
       {
         retcode = m_mshdist.run(
             boxp.string(),
             "-dom",
             "-surf",
             "-fmm",
             "-ncpu", m_ncpu,
             "-v 0"
             );
       }
       else
       {
         if (m_distTheBoundary)
         {
           retcode = m_mshdist.run(
               boxp.string(),
               "-surf",
               "-fmm",
               "-ncpu", m_ncpu,
               "-v 0"
               );
         }
         else
         {
           retcode = m_mshdist.run(
               boxp.string(),
               "-dom",
               "-fmm",
               "-ncpu", m_ncpu
               );
         }
       }
       ScalarGridFunction res(m_fes.get());
       if (retcode != 0)
         Alert::Exception() << "ISCD::Mshdist invocation failed." << Alert::Raise;
       else
         res.load(boxp.replace_extension(".sol"), IO::FileFormat::MEDIT);
       return res;
     }

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
     ScalarGridFunction distance(
         const Geometry::Mesh<Context::Serial>& box,
         const Geometry::Mesh<Context::Serial>& contour)
     {
       auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
       box.save(boxp);

       auto contourp = m_mshdist.tmpnam(".mesh", "RodinMMG");
       contour.save(contourp, IO::FileFormat::MEDIT);

       int retcode = 1;
       if (box.isSurface())
       {
         retcode = m_mshdist.run(
             boxp.string(), contourp.string(),
             "-surf",
             m_scale ? "" : "-noscale",
             "-ncpu", m_ncpu,
             "-v 0"
             );
       }
       else
       {
         retcode = m_mshdist.run(
             boxp.string(), contourp.string(),
             m_distTheBoundary ? "-surf" : "",
             m_scale ? "" : "-noscale",
             "-ncpu", m_ncpu,
             "-v 0");
       }

       Variational::GridFunction res(m_fes.get());
       if (retcode != 0)
         Alert::Exception() << "ISCD::Mshdist invocation failed." << Alert::Raise;
       else
         res.load(boxp.replace_extension(".sol"), IO::FileFormat::MEDIT);
       return res;
     }

     Distancer& surface(bool distTheBoundary = true)
     {
       m_distTheBoundary = distTheBoundary;
       return *this;
     }

    private:
     std::reference_wrapper<const ScalarGridFunction::FES> m_fes;
     bool m_scale;
     bool m_distTheBoundary;
     unsigned int m_ncpu;
     ISCDProcess m_mshdist;
     std::optional<FlatSet<Geometry::Attribute>> m_interiorDomains;
  };

  Distancer(const ScalarGridFunction::FES&) -> Distancer<ScalarGridFunction::FES>;

  // /**
  //  * @brief Distancer specialization for redistancing a level set function.
  //  */
  // template <>
  // class Distancer<void>
  // {
  //  public:
  //    Distancer()
  //     : m_ncpu(std::thread::hardware_concurrency()),
  //      m_mshdist(getISCDMshdistExecutable()),
  //      m_distTheBoundary(false)
  //    {}

  //    /**
  //    * @brief Redistances the level set function.
  //    *
  //    * Given a level set function defined on the vertices of a bounding box
  //    * @f$ D @f$, this method will regenerate the signed distance function
  //    * associated to the the subdomain @f$ \Omega \subset D @f$ which the
  //    * level set function represents.
  //    *
  //    * @param[in,out] sol Level set function
  //    */
  //    void redistance(Variational::GridFunction<Variational::H1<Context::Serial>>& sol)
  //    {
  //     auto meshp = m_mshdist.tmpnam(".mesh", "RodinMMG");
  //     sol.getFiniteElementSpace().getMesh().save(meshp, IO::FileFormat::MEDIT);

  //     boost::filesystem::path solp(meshp);
  //     solp.replace_extension(".sol");
  //     sol.save(solp, IO::FileFormat::MEDIT);

  //     auto name = solp;
  //     name.replace_extension();

  //     int retcode = 1;
  //     if (sol.getFiniteElementSpace().getMesh().isSurface())
  //     {
  //      retcode = m_mshdist.run(name.string(),
  //         "-surf -fmm",
  //         "-ncpu", m_ncpu,
  //         "-v 0");
  //     }
  //     else
  //     {
  //      retcode = m_mshdist.run(name.string(),
  //         m_distTheBoundary ? "-surf" : ""
  //         "-ncpu", m_ncpu,
  //         "-v 0");
  //     }

  //     if (retcode != 0)
  //      Alert::Exception("ISCD::Mshdist invocation failed.").raise();
  //     sol.load(solp, IO::FileFormat::MEDIT);
  //    }

  //    Distancer& surface(bool distTheBoundary = true)
  //    {
  //     m_distTheBoundary = distTheBoundary;
  //     return *this;
  //    }

  //    /**
  //    * @brief Specifies how many CPUs to use when distancing in parallel.
  //    * @param[in] ncpu Number of CPUs to use
  //    *
  //    * By default, it will utilize `std::thread::hardware_concurrency()`.
  //    */
  //    Distancer& setCPUs(unsigned int ncpu);

  //    /**
  //    * @returns Number of CPUs which will be used when performing the
  //    * distancing.
  //    */
  //    unsigned int getCPUs() const;

  //  private:
  //    unsigned int m_ncpu;
  //    ISCDProcess m_mshdist;
  //    bool m_distTheBoundary;
  // };
  // Distancer() -> Distancer<void>;


}

#endif
