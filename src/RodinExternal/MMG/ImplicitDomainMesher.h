/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHER_H
#define RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHER_H

#include "Mesh.h"

#include "MMG5.h"
#include "Common.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Class to perform the discretization and optimization of a
   * surface implicitly defined by a level set function.
   *
   * This class performs the discretization of a level set function @f$\phi : D
   * \rightarrow \Omega @f$ defined over a computational mesh @f$ D @f$ which
   * contains a regular domain @f$ \Omega @f$.
   * @f[
   *  \left\{
   *    \begin{aligned}
   *      \phi(x) < 0 && \text{if } \ x &\in \Omega\\
   *      \phi(x) = 0 && \text{if } \ x &\in \partial \Omega\\
   *      \phi(x) > 0 && \text{if } \ x &\in D \backslash \overline{\Omega}\\
   *    \end{aligned}
   *  \right.
   * @f]
   */
  class ImplicitDomainMesher : public MMG5
  {
    public:
      /**
       * @brief Default constructor.
       */
      ImplicitDomainMesher()
        : m_ls(0.0),
          m_meshTheSurface(false)
      {}

      ImplicitDomainMesher& surface(bool meshTheSurface = true);

      /**
       * @brief Specifies the level set to discretize.
       *
       * The default value is 0.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& setLevelSet(Scalar ls);

      /**
       * @brief Specifies the removal of small parasitic components.
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& setRMC(Scalar rmc = 1e-5);

      ImplicitDomainMesher& setBaseReferences(Geometry::Attribute ref)
      {
        return setBaseReferences(FlatSet<Geometry::Attribute>{ref});
      }

      ImplicitDomainMesher& setBaseReferences(const FlatSet<Geometry::Attribute>& refs);

      /**
       * @brief Sets the material reference for the discretized boundary
       * @f$ \partial \Omega @f$
       *
       * Specifies what the label for the boundary @f$ \partial \Omega @f$ will
       * be. By default, the value is 10.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& setBoundaryReference(const Geometry::Attribute& ref);

      /**
       * @brief Specifies how to split the materials into an interior and
       * exterior domains.
       *
       * This map specifies for each input material reference the values of the
       * 2 new domains created by the level-set splitting. By default, the
       * material references of the interior and exterior, are 2 and 3
       * respectively.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& setSplit(const SplitMap& split);

      /**
       * @brief Indicates that a material reference should be split.
       * @param[in] ref Material to split
       * @param[in] s Interior and exterior labels
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& split(const Geometry::Attribute& ref, const Split& s);

      /**
       * @brief Indicates that a material reference should not be split.
       * @param[in] ref Material to split
       * @param[in] s Interior and exterior labels
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher& noSplit(const Geometry::Attribute& ref);

      /**
       * @brief Discretizes and optimizes an implicitly defined surface defined
       * by a level set function.
       * @param[in] ls Level set function
       *
       * The material reference of the level set (edge) boundary will be 10.
       */
      MMG::Mesh discretize(const MMG::ScalarGridFunction& ls);

      ImplicitDomainMesher& setAngleDetection(bool enable = true)
      {
        MMG5::setAngleDetection(enable);
        return *this;
      }

      ImplicitDomainMesher& setHMin(Scalar hmin)
      {
        MMG5::setHMin(hmin);
        return *this;
      }

      ImplicitDomainMesher& setHMax(Scalar hmax)
      {
        MMG5::setHMax(hmax);
        return *this;
      }

      ImplicitDomainMesher& setHausdorff(Scalar hausd)
      {
        MMG5::setHausdorff(hausd);
        return *this;
      }

      ImplicitDomainMesher& setGradation(Scalar hgrad)
      {
        MMG5::setGradation(hgrad);
        return *this;
      }

      const SplitMap& getSplitMap() const
      {
        return m_split;
      }

    private:
      ReturnCode discretizeMMG2D(MMG5_pMesh mesh, MMG5_pSol sol);
      ReturnCode discretizeMMG3D(MMG5_pMesh mesh, MMG5_pSol sol);
      ReturnCode discretizeMMGS(MMG5_pMesh mesh, MMG5_pSol sol);

      void generateUniqueSplit(const FlatSet<Geometry::Attribute>& attr);

      void deleteBoundaryRef(MMG5_pMesh mesh, Geometry::Attribute ref);

      Scalar m_ls;
      SplitMap m_split;
      bool m_meshTheSurface;
      std::optional<Scalar> m_rmc;
      FlatSet<Geometry::Attribute> m_lsBaseReferences;
      std::optional<Geometry::Attribute> m_isoref;

      SplitMap m_uniqueSplit;

      /**
       * @internal
       *
       * Generated attribute to original attribute map.
       */
      UnorderedMap<Geometry::Attribute, Geometry::Attribute> m_g2om;
  };
}
#endif
