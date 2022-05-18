/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_VECTORSOLUTION3D_H
#define RODIN_EXTERNAL_MMG_VECTORSOLUTION3D_H

#include <cstddef>
#include <iterator>
#include <functional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "Mesh3D.h"

#include "VectorSolution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Vector solution supported on a 3D mesh.
   *
   * An object of type VectorSolution3D represents a function
   * @f[
   * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
   * @f]
   * whose known values are given on vertices of some mesh @f$ \Omega @f$.
   */
  class VectorSolution3D :  public VectorSolution
  {
    public:
      /**
       * @brief Reads the solution text file.
       *
       * The file is read using MMGv2 format.
       *
       * @param[in] filename Name of file to read.
       */
      VectorSolution3D& load(const boost::filesystem::path& filename) override;

      /**
       * @internal
       */
      VectorSolution3D(MMG5_pSol sol, Mesh3D& mesh);

      /**
       * @brief Initializes the object with no data
       *
       * @param[in] mesh Reference to the underlying mesh.
       */
      VectorSolution3D(Mesh3D& mesh);

      /**
       * @brief Performs a move construction from the `other` solution object.
       *
       * @param[in] other Object to move.
       */
      VectorSolution3D(VectorSolution3D&& other);

      /**
       * @brief Performs a copy of the `other` solution object.
       *
       * @param[in] other Object to copy.
       * @note It does not perform a copy the mesh. Instead the new object
       * will have a reference to the original mesh.
       */
      VectorSolution3D(const VectorSolution3D& other);

      /**
       * @brief Move assigns the `other` solution object to this object.
       *
       * @param[in] other Object to move.
       */
      VectorSolution3D& operator=(VectorSolution3D&& other) = default;

      /**
       * @brief Copy assigns the `other` solution object to this object.
       *
       * @param[in] other Object to copy.
       */
      VectorSolution3D& operator=(const VectorSolution3D& other);

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
      VectorSolution3D& setMesh(Mesh3D& mesh);

      /**
       * @brief Gets the constant reference to the underlying mesh.
       *
       * @returns Constant reference to the underlying mesh.
       */
      const Mesh3D& getMesh() const;

      /**
       * @brief Gets the reference to the underlying mesh.
       *
       * @returns Reference to the underlying mesh.
       */
      Mesh3D& getMesh();

      void save(const boost::filesystem::path& filename) override;

    private:
      std::reference_wrapper<Mesh3D> m_mesh;
  };
}
#endif

