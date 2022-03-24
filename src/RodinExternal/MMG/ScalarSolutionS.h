/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTIONS_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTIONS_H

#include "ForwardDecls.h"

#include "MeshS.h"
#include "ScalarSolution.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Scalar solution supported on a surface mesh.
    *
    * An object of type ScalarSolutionS represents a function
    * @f[
    * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
    * @f]
    * whose known values are given on vertices of some surfacic mesh @f$ \Omega @f$.
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
         static IncompleteScalarSolutionS load(const std::filesystem::path& filename);

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
          * @brief Frees the data.
          */
         virtual ~ScalarSolutionS();

         /**
          * @brief Move assigns the `other` solution object to this object.
          *
          * @param[in] other Object to move.
          */
         ScalarSolutionS& operator=(ScalarSolutionS&& other);

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

         void save(const std::filesystem::path& filename) override;

         MMG5_pSol& getHandle() override;

         const MMG5_pSol& getHandle() const override;

      private:
         std::reference_wrapper<MeshS> m_mesh;
         MMG5_pSol m_sol;
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
   class IncompleteScalarSolutionS
   {
      public:
         /**
          * @brief Constructs an empty scalar solution object without a mesh.
          */
         IncompleteScalarSolutionS();

         /**
          * @brief Constructs a scalar solution with `n` unitialized entries.
          * @param[in] n Number of entries that the solution has.
          */
         IncompleteScalarSolutionS(int n);

         /**
          * @brief Frees the data if it still owns the data, i.e. the
          * setMesh(MeshS&) method has not been called.
          */
         virtual ~IncompleteScalarSolutionS();

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
         ScalarSolutionS setMesh(MeshS& mesh);

         /**
          * @internal
          * @returns Reference to underlying solution handle.
          */
         MMG5_pSol& getHandle();

         /**
          * @internal
          * @returns Constant reference to underlying solution handle.
          */
         const MMG5_pSol& getHandle() const;

      private:
         MMG5_pSol m_sol;
         bool m_isOwner;
   };
}

#endif
