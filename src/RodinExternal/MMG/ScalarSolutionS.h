/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTIONS_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTIONS_H

#include <cstddef>
#include <iterator>
#include <functional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "MeshS.h"

#include "ScalarSolution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Scalar solution supported on a S mesh.
   *
   * An object of type ScalarSolutionS represents a function
   * @f[
   * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
   * @f]
   * whose known values are given on vertices of some mesh @f$ \Omega @f$.
   */
  class ScalarSolutionS :  public ScalarSolution
  {
    public:
      /**
       * @brief Reads the solution text file.
       *
       * The file is read using MMGv2 format.
       *
       * @param[in] filename Name of file to read.
       */
      static IncompleteScalarSolutionS load(const boost::filesystem::path& filename);

      /**
       * @internal
       */
      ScalarSolutionS(MMG5_pSol sol, MeshS& mesh);

      /**
       * @brief Initializes the object with no data
       *
       * @param[in] mesh Reference to the underlying mesh.
       */
      ScalarSolutionS(MeshS& mesh);

      /**
       * @brief Performs a move construction from the `other` solution object.
       *
       * @param[in] other Object to move.
       */
      ScalarSolutionS(ScalarSolutionS&& other);

      /**
       * @brief Performs a copy of the `other` solution object.
       *
       * @param[in] other Object to copy.
       * @note It does not perform a copy the mesh. Instead the new object
       * will have a reference to the original mesh.
       */
      ScalarSolutionS(const ScalarSolutionS& other);

      /**
       * @brief Move assigns the `other` solution object to this object.
       *
       * @param[in] other Object to move.
       */
      ScalarSolutionS& operator=(ScalarSolutionS&& other) = default;

      /**
       * @brief Copy assigns the `other` solution object to this object.
       *
       * @param[in] other Object to copy.
       */
      ScalarSolutionS& operator=(const ScalarSolutionS& other);

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
      ScalarSolutionS& setMesh(MeshS& mesh);

      /**
       * @brief Gets the constant reference to the underlying mesh.
       *
       * @returns Constant reference to the underlying mesh.
       */
      const MeshS& getMesh() const;

      /**
       * @brief Gets the reference to the underlying mesh.
       *
       * @returns Reference to the underlying mesh.
       */
      MeshS& getMesh();

      void save(const boost::filesystem::path& filename) override;

    private:
      std::reference_wrapper<MeshS> m_mesh;
  };

  /**
   * @brief A scalar solution which does not have a mesh assigned to it.
   *
   * To unlock the full functionality of the class you must call the
   * @ref setMesh(MeshS&) method. For example, when loading it from file:
   *
   * @code{.cpp}
   * auto sol = ScalarSolutionS::load(filename).setMesh(mesh);
   * @endcode
   */
  class IncompleteScalarSolutionS : public IncompleteSolutionBase
  {
    public:
      IncompleteScalarSolutionS();

      IncompleteScalarSolutionS(const IncompleteScalarSolutionS& other)
        : IncompleteSolutionBase(other)
      {}

      IncompleteScalarSolutionS(IncompleteScalarSolutionS&& other)
        : IncompleteSolutionBase(std::move(other))
      {}

      /**
       * @brief Sets the associated mesh and moves ownership to the new
       * object.
       *
       * @param[in] mesh Reference to mesh.
       *
       * @returns An object of type ScalarSolutionS which represents the
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
      ScalarSolutionS setMesh(MeshS& mesh)
      {
        return ScalarSolutionS(release(), mesh);
      }
  };
}
#endif
