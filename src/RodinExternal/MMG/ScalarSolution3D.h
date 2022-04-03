/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTION3D_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTION3D_H

#include <cstddef>
#include <iterator>
#include <functional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "Mesh3D.h"

#include "ScalarSolution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Scalar solution supported on a 3D mesh.
   *
   * An object of type ScalarSolution3D represents a function
   * @f[
   * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
   * @f]
   * whose known values are given on vertices of some mesh @f$ \Omega @f$.
   */
  class ScalarSolution3D :  public ScalarSolution
  {
    public:
      /**
       * @brief Reads the solution text file.
       *
       * The file is read using MMGv2 format.
       *
       * @param[in] filename Name of file to read.
       */
      static IncompleteScalarSolution3D load(const boost::filesystem::path& filename);

      /**
       * @internal
       */
      ScalarSolution3D(MMG5_pSol sol, Mesh3D& mesh);

      /**
       * @brief Initializes the object with no data
       *
       * @param[in] mesh Reference to the underlying mesh.
       */
      ScalarSolution3D(Mesh3D& mesh);

      /**
       * @brief Performs a move construction from the `other` solution object.
       *
       * @param[in] other Object to move.
       */
      ScalarSolution3D(ScalarSolution3D&& other);

      /**
       * @brief Performs a copy of the `other` solution object.
       *
       * @param[in] other Object to copy.
       * @note It does not perform a copy the mesh. Instead the new object
       * will have a reference to the original mesh.
       */
      ScalarSolution3D(const ScalarSolution3D& other);

      /**
       * @brief Move assigns the `other` solution object to this object.
       *
       * @param[in] other Object to move.
       */
      ScalarSolution3D& operator=(ScalarSolution3D&& other) = default;

      /**
       * @brief Copy assigns the `other` solution object to this object.
       *
       * @param[in] other Object to copy.
       */
      ScalarSolution3D& operator=(const ScalarSolution3D& other);

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
      ScalarSolution3D& setMesh(Mesh3D& mesh);

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

  /**
   * @brief A scalar solution which does not have a mesh assigned to it.
   *
   * To unlock the full functionality of the class you must call the
   * @ref setMesh(Mesh3D&) method. For example, when loading it from file:
   *
   * @code{.cpp}
   * auto sol = ScalarSolution3D::load(filename).setMesh(mesh);
   * @endcode
   */
  class IncompleteScalarSolution3D : public IncompleteSolutionBase
  {
    public:
      IncompleteScalarSolution3D();

      IncompleteScalarSolution3D(const IncompleteScalarSolution3D& other)
        : IncompleteSolutionBase(other)
      {}

      IncompleteScalarSolution3D(IncompleteScalarSolution3D&& other)
        : IncompleteSolutionBase(std::move(other))
      {}

      /**
       * @brief Sets the associated mesh and moves ownership to the new
       * object.
       *
       * @param[in] mesh Reference to mesh.
       *
       * @returns An object of type ScalarSolution3D which represents the
       * object with all its functionality.
       *
       * @note The method does not incur any significant performance penalty
       * since no data is copied.
       *
       * @warning The method does not check to see if the mesh is compatible
       * with the current data in the solution. In general, it is up to the
       * user to ensure that the number of points are the same and keep
       * track of the modifications to the underlying mesh.
       */
      ScalarSolution3D setMesh(Mesh3D& mesh)
      {
        return ScalarSolution3D(release(), mesh);
      }
  };
}
#endif
