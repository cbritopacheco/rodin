/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_MESHOPTIMIZER2D_H
#define RODIN_RODININTEGRATION_MMG_MESHOPTIMIZER2D_H

#include <optional>

#include "Rodin/Alert.h"

#include "Mesh2D.h"
#include "ScalarSolution2D.h"

#include "MMG2D.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Improves the mesh quality, mantaining the mean edge lenghts of the mesh.
   */
  class MeshOptimizer2D : public MMG2D
  {
    public:
      /**
       * @brief Performs the optimization of the mesh.
       * @param[in, out] mesh Mesh to optimize
       * @returns The sizemap computed to preserve edge sizes
       * @note The mean of the edge lengths is preserved at the vertices.
       * Hence, if the edges passing through a vertex have very different
       * sizes, the resulting mesh may be very different from the initial one.
       */
      ScalarSolution2D optimize(Mesh2D& mesh);

      MeshOptimizer2D& setHMin(double hmin) override;
      MeshOptimizer2D& setHMax(double hmax) override;
      MeshOptimizer2D& setHausdorff(double hausd) override;
      MeshOptimizer2D& setGradation(double hgrad) override;

    private:
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd;
  };
}

#endif
