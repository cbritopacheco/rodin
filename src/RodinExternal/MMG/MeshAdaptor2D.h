/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHADAPTOR2D_H
#define RODIN_EXTERNAL_MMG_MESHADAPTOR2D_H

#include <memory>

#include "Mesh2D.h"
#include "ScalarSolution2D.h"

#include "MMG2D.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Adapts a Mesh2D object according to some given metric.
   *
   * The MeshAdaptor2D can adapt a Mesh2D object according to some metric and
   * parameters.
   *
   * @see @ref examples-mmg-mesh_adaptation
   *
   * @see [Mesh adaptation to a solution](https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-tutorials/mmg-remesher-mmg2d/mesh-adaptation-to-a-solution)
   */
  class MeshAdaptor2D : public MMG2D
  {
    public:
      /**
       * @brief Creates a MeshAdaptor2D object with default parameters.
       */
      MeshAdaptor2D() = default;

      /**
       * @brief Sets the metric to which the mesh will be adapted to.
       *
       * @param[in] metric Size map
       *
       * @returns Reference to self (for method chaining)
       *
       * The metric is a size map which is used to impose a desired size
       * feature when remeshing the input mesh.  This size map is a scalar
       * function defined at the mesh vertices.  At each vertex, it associates
       * the desired size of the surrounding elements of the mesh.
       */
      MeshAdaptor2D& setMetric(const ScalarSolution2D& metric);

      /**
       * @brief Adapts the mesh given the previously set parameters.
       *
       * @param[in,out] mesh Mesh to adapt
       *
       * This function performs the actual adaptation of the mesh.
       */
      void adapt(Mesh2D& mesh) const;

      MeshAdaptor2D& setHMin(double hmin) override;
      MeshAdaptor2D& setHMax(double hmax) override;
      MeshAdaptor2D& setHausdorff(double hausd) override;
      MeshAdaptor2D& setGradation(double hgrad) override;

    private:
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd;

      std::optional<std::variant<ScalarSolution2D>> m_metric;
  };
}

#endif
