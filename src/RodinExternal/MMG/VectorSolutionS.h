/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_VECTORSOLUTIONS_H
#define RODIN_EXTERNAL_MMG_VECTORSOLUTIONS_H

#include <cstddef>
#include <iterator>
#include <functional>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "MeshS.h"

#include "VectorSolution.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Vector solution supported on a S mesh.
   *
   * An object of type VectorSolutionS represents a function
   * @f[
   * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
   * @f]
   * whose known values are given on vertices of some mesh @f$ \Omega @f$.
   */
  class VectorSolutionS :  public VectorSolution
  {
    public:
      /**
       * @brief Reads the solution text file.
       *
       * The file is read using MMGv2 format.
       *
       * @param[in] filename Name of file to read.
       */
      static IncompleteVectorSolutionS load(const boost::filesystem::path& filename);

      /**
       * @internal
       */
      VectorSolutionS(MMG5_pSol sol, MeshS& mesh);

      /**
       * @brief Initializes the object with no data
       *
       * @param[in] mesh Reference to the underlying mesh.
       */
      VectorSolutionS(MeshS& mesh);

      /**
       * @brief Performs a move construction from the `other` solution object.
       *
       * @param[in] other Object to move.
       */
      VectorSolutionS(VectorSolutionS&& other);

      /**
       * @brief Performs a copy of the `other` solution object.
       *
       * @param[in] other Object to copy.
       * @note It does not perform a copy the mesh. Instead the new object
       * will have a reference to the original mesh.
       */
      VectorSolutionS(const VectorSolutionS& other);

      /**
       * @brief Move assigns the `other` solution object to this object.
       *
       * @param[in] other Object to move.
       */
      VectorSolutionS& operator=(VectorSolutionS&& other) = default;

      /**
       * @brief Copy assigns the `other` solution object to this object.
       *
       * @param[in] other Object to copy.
       */
      VectorSolutionS& operator=(const VectorSolutionS& other);

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
      VectorSolutionS& setMesh(MeshS& mesh);

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
   * @brief A Vector solution which does not have a mesh assigned to it.
   *
   * To unlock the full functionality of the class you must call the
   * @ref setMesh(MeshS&) method. For example, when loading it from file:
   *
   * @code{.cpp}
   * auto sol = VectorSolutionS::load(filename).setMesh(mesh);
   * @endcode
   */
  class IncompleteVectorSolutionS : public IncompleteSolutionBase
  {
    public:
      IncompleteVectorSolutionS();

      IncompleteVectorSolutionS(const IncompleteVectorSolutionS& other)
        : IncompleteSolutionBase(other)
      {}

      IncompleteVectorSolutionS(IncompleteVectorSolutionS&& other)
        : IncompleteSolutionBase(std::move(other))
      {}

      /**
       * @brief Sets the associated mesh and moves ownership to the new
       * object.
       *
       * @param[in] mesh Reference to mesh.
       *
       * @returns An object of type VectorSolutionS which represents the
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
      VectorSolutionS setMesh(MeshS& mesh)
      {
        return VectorSolutionS(release(), mesh);
      }
  };
}
#endif
