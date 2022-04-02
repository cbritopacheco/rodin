/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHERS_H
#define RODIN_EXTERNAL_MMG_IMPLICITDOMAINMESHERS_H

#include <utility>

#include "ForwardDecls.h"
#include "Common.h"
#include "MeshS.h"
#include "ScalarSolutionS.h"

#include "MMGS.h"

namespace Rodin::External::MMG
{

  /**
   * @brief Class to perform the discretization and optimization of a
   * surface implicitly defined by a level set function.
   */
  class ImplicitDomainMesherS : public MMGS
  {
    public:
      /**
       * @brief Constructs an ImplicitDomainMesherS with default values.
       */
      ImplicitDomainMesherS();

      /**
       * @brief Specifies the level set to discretize.
       *
       * The default value is 0.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesherS& setLevelSet(double ls);

      /**
       * @brief Sets the material reference for the discretized boundary
       * @f$ \partial \Omega @f$
       *
       * Specifies what the label for the boundary @f$ \partial \Omega @f$ will
       * be. By default, the value is 10.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesherS& setBoundaryReference(const MaterialReference& ref);

      /**
       * @brief Discretizes and optimizes an implicitly defined surface defined
       * by a level set function.
       *
       * @param[in] ls Level set function
       * @returns Discretization
       */
      MeshS discretize(ScalarSolutionS& ls);

      ImplicitDomainMesherS& setHMin(double hmin) override;
      ImplicitDomainMesherS& setHMax(double hmax) override;
      ImplicitDomainMesherS& setHausdorff(double hausd) override;
      ImplicitDomainMesherS& setGradation(double hgrad) override;

    private:
      double m_ls;
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd;
      std::optional<MaterialReference> m_isoref;
  };
}

#endif

