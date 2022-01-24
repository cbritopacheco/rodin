/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHER2D_H
#define RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHER2D_H

#include <utility>

#include "ForwardDecls.h"
#include "Common.h"
#include "Mesh2D.h"
#include "ScalarSolution2D.h"

#include "MMG2D.h"

namespace Rodin::External::MMG
{

  /**
   * @brief Class to perform the discretization and optimization of a
   * surface implicitly defined by a level set function.
   */
  class ImplicitDomainMesher2D : public MMG2D
  {
    public:
      /**
       * @brief Discretized mesh and solution.
       */
      struct Discretization
      {
        Mesh2D mesh; /// Discretized mesh
        ScalarSolution2D solution; /// Level set function defined on the new vertices
      };

      /**
       * @brief Constructs an ImplicitDomainMesher2D with default values.
       */
      ImplicitDomainMesher2D();

      /**
       * @brief Specifies the level set to discretize.
       *
       * The default value is 0.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& setLevelSet(double ls);

      /**
       * @brief Specifies the removal of small parasitic components.
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& setRMC(double rmc = 1e-5);

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
      ImplicitDomainMesher2D& setSplit(const SplitMap& split);

      /**
       * @brief Sets the material reference for the discretized boundary
       * @f$ \partial \Omega @f$
       *
       * Specifies what the label for the boundary @f$ \partial \Omega @f$ will
       * be. By default, the value is 10.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& setBoundaryReference(const MaterialReference& ref);

      /**
       * @brief Gets the split map
       * @returns SplitMap
       */
      const SplitMap& getSplit() const;

      /**
       * @brief Indicates that a material reference should be split.
       * @param[in] ref Material to split
       * @param[in] s Interior and exterior labels
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& split(const MaterialReference& ref, const Split& s);

      /**
       * @brief Indicates that a material reference should not be split.
       * @param[in] ref Material to split
       * @param[in] s Interior and exterior labels
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& noSplit(const MaterialReference& ref);

      /**
       * @brief Discretizes and optimizes an implicitly defined surface defined
       * by a level set function.
       *
       * The material reference of the level set (edge) boundary will be 10.
       *
       * @param[in] ls Level set function
       * @returns Discretization
       */
      Discretization discretize(ScalarSolution2D& ls);

      ImplicitDomainMesher2D& setHMin(double hmin) override;
      ImplicitDomainMesher2D& setHMax(double hmax) override;
      ImplicitDomainMesher2D& setHausdorff(double hausd) override;
      ImplicitDomainMesher2D& setGradation(double hgrad) override;

    private:
      double m_ls;
      SplitMap m_split;
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd,
                            m_rmc;
      std::optional<MaterialReference> m_isoref;

  };
}

#endif
