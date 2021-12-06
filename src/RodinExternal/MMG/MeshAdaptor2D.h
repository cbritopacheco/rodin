/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_MESHADAPTOR2D_H
#define RODIN_EXTERNAL_MMG_MESHADAPTOR2D_H

#include <memory>

#include "MeshAdaptor.h"
#include "Mesh2D.h"
#include "Solution.h"
#include "ScalarSolution2D.h"

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
  class MeshAdaptor2D : public MeshAdaptor<Mesh2D, MeshAdaptor2D>
  {
    public:
      /**
       * @brief Creates a MeshAdaptor2D object with default parameters.
       */
      MeshAdaptor2D() = default;

      /**
       * @brief Sets the minimal edge size.
       *
       * @param[in] hmin Minimal edge size.
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hmin` option allows to truncate the edge sizes to be lower than the
       * `hmax` one.
       *
       * The default values for this parameters are computed from the mesh
       * bounding box or, if provided, from the given metric.
       *
       * - Without metric, the minimal edge size is set to 0.01 of the bounding
       * box size.
       *
       * - With metric, the minimal edge size is set to 0.1 of the
       * smallest prescribed size.
       *
       * @see setHMax(double)
       * @see [EXTERNAL: -hmin / -hmax](https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options/hmin-hmax)
       */
      MeshAdaptor2D& setHMin(double hmin);

      /**
       * @brief Sets the maximal edge size parameter.
       *
       * @param[in] hmax Maximal edge size.
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hmax` option allows to truncate the edge sizes to be greater than
       * the `hmin` parameter.
       *
       * The default values for this parameters are computed from the mesh
       * bounding box or, if provided, from the given metric.
       *
       * - Without metric, the maximal edge size is set to two times the
       *   bounding box size.
       *
       * - With metric, the maximal one is set to 10 times the maximal
       *   prescribed size.
       *
       * @see setHMin(double)
       * @see [EXTERNAL: -hmin / -hmax](https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options/hmin-hmax)
       */
      MeshAdaptor2D& setHMax(double hmax);

      /**
       * @brief Sets the Hausdorff parameter.
       *
       * @param[in] hausd Hausdorff parameter.
       *
       * @returns Reference to self (for method chaining)
       *
       * The Hausdorff parameter controls the boundary approximation.  It
       * imposes the maximal distance between the piecewise linear
       * representation of the boundary and the reconstructed ideal boundary.
       * Thus, a low Hausdorff parameter leads to the refinement of high
       * curvature areas.
       *
       * By default, the Hausdorff value is set to 0.01, which is a suitable
       * value for an object of size 1 in each direction. For smaller (resp.
       * larger) objects, you may need to decrease (resp. increase) the
       * Hausdorff parameter.
       *
       * @see [EXTERNAL: -hausd](https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options/mmg-remesher-option-hausd)
       */
      MeshAdaptor2D& setHausdorff(double hausd);

      /**
       * @brief Sets the gradation parameter
       *
       * @param[in] hgrad Gradation parameter
       *
       * @returns Reference to self (for method chaining)
       *
       * The `hgrad` option allows to set the gradation value. It controls the
       * ratio between two adjacent edges. With a gradation of @f$ h @f$, two
       * adjacent edges @f$ e_1 @f$ and @f$ e_2 @f$ must respect the following
       * constraint:
       *
       * @f[
       *  \dfrac{1}{h} \leq \dfrac{ |e_1| }{ |e_2| } \leq h
       * @f]
       *
       * By default, the gradation value is 1.3.
       *
       * @see [EXTERNAL: -hgrad](https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options/mmg-remesher-option-hgrad)
       */
      MeshAdaptor2D& setGradation(double hgrad);

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
      MeshAdaptor2D& setMetric(const ScalarSolution2D<>& metric);

      /**
       * @brief Adapts the mesh given the previously set parameters.
       *
       * @param[in,out] mesh Mesh to adapt
       *
       * This function performs the actual adaptation of the mesh.
       */
      void adapt(Mesh2D& mesh) const;

    private:
      std::optional<double> m_hmin,
                            m_hmax,
                            m_hgrad,
                            m_hausd;

      std::optional<std::variant<ScalarSolution2D<>>> m_metric;
  };
}

#endif
