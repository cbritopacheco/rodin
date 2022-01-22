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

#include "Mesh2D.h"
#include "ScalarSolution2D.h"

#include "MMG2D.h"

namespace Rodin::External::MMG
{
  class ImplicitDomainMesher2D : public MMG2D
  {
    public:
      struct Discretization
      {
        Mesh2D mesh;
        ScalarSolution2D solution;
      };

      ImplicitDomainMesher2D();

      /**
       * @brief Specifies the level set to discretize.
       *
       * The default value is 0.
       *
       * @returns Reference to self (for method chaining)
       */
      ImplicitDomainMesher2D& setLevelSet(double ls);

      Discretization discretize(ScalarSolution2D& ls);

      ImplicitDomainMesher2D& setHMin(double hmin) override;
      ImplicitDomainMesher2D& setHMax(double hmax) override;
      ImplicitDomainMesher2D& setHausdorff(double hausd) override;
      ImplicitDomainMesher2D& setGradation(double hgrad) override;

    private:
      double m_ls;
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd;
  };
}

#endif
