/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_VECTORSOLUTION2D_H
#define RODIN_EXTERNAL_MMG_VECTORSOLUTION2D_H

#include <cstddef>
#include <iterator>
#include <functional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "Mesh2D.h"

#include "VectorSolution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Vector solution supported on a 2D mesh.
   *
   * An object of type VectorSolution2D represents a function
   * @f[
   * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
   * @f]
   * whose known values are given on vertices of some mesh @f$ \Omega @f$.
   */
  class VectorSolution2D :  public VectorSolution
  {
    public:
      /**
       * @brief Reads the solution text file.
       *
       * The file is read using MMGv2 format.
       *
       * @param[in] filename Name of file to read.
       */
      VectorSolution2D& load(const boost::filesystem::path& filename) override;

      /**
       * @internal
       */
      VectorSolution2D(MMG5_pSol sol, Mesh2D& mesh);

      /**
       * @brief Initializes the object with no data
       *
       * @param[in] mesh Reference to the underlying mesh.
       */
      VectorSolution2D(Mesh2D& mesh);

      /**
       * @brief Performs a move construction from the `other` solution object.
       *
       * @param[in] other Object to move.
       */
      VectorSolution2D(VectorSolution2D&& other);

      /**
       * @brief Performs a copy of the `other` solution object.
       *
       * @param[in] other Object to copy.
       * @note It does not perform a copy the mesh. Instead the new object
       * will have a reference to the original mesh.
       */
      VectorSolution2D(const VectorSolution2D& other);

      /**
       * @brief Move assigns the `other` solution object to this object.
       *
       * @param[in] other Object to move.
       */
      VectorSolution2D& operator=(VectorSolution2D&& other) = default;

      /**
       * @brief Copy assigns the `other` solution object to this object.
       *
       * @param[in] other Object to copy.
       */
      VectorSolution2D& operator=(const VectorSolution2D& other);

      /**
       * @brief Sets the associated mesh.
       *
       * @param[in] mesh Reference to mesh.
       *
       * @returns Reference to self (for method chaining).
       *
       * @note The method does not check to see if the mesh is compatible
       * with the current data in the solution. In general, it is up to the
       * user to ensure that the number of points are the same, keep track
       * of the modifications to the underlying mesh, etc.
       */
      VectorSolution2D& setMesh(Mesh2D& mesh);

      /**
       * @brief Gets the constant reference to the underlying mesh.
       *
       * @returns Constant reference to the underlying mesh.
       */
      const Mesh2D& getMesh() const;

      /**
       * @brief Gets the reference to the underlying mesh.
       *
       * @returns Reference to the underlying mesh.
       */
      Mesh2D& getMesh();

      void save(const boost::filesystem::path& filename) override;

    private:
      std::reference_wrapper<Mesh2D> m_mesh;
  };
}
#endif
