/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_VECTORSOLUTION2D_H
#define RODIN_EXTERNAL_MMG_VECTORSOLUTION2D_H

#include <filesystem>

#include "ForwardDecls.h"
#include "VectorSolution.h"

namespace Rodin::External::MMG
{
  class VectorSolution2D : public VectorSolution
  {
    public:
         /**
          * @brief Reads the solution text file.
          *
          * The file is read using MMGv2 format.
          *
          * @param[in] filename Name of file to read.
          */
         static IncompleteVectorSolution2D load(const std::filesystem::path& filename);

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
          * @brief Frees the data.
          */
         virtual ~VectorSolution2D();

         /**
          * @brief Move assigns the `other` solution object to this object.
          *
          * @param[in] other Object to move.
          */
         VectorSolution2D& operator=(VectorSolution2D&& other);

         /**
          * @brief Copy assigns the `other` solution object to this object.
          *
          * @param[in] other Object to copy.
          */
         VectorSolution2D& operator=(const VectorSolution2D& other);

         /**
          * @brief Write the solution to a text file.
          *
          * The file is written using MMGv2 format.
          *
          * @param[in] filename Name of file to write.
          */
         void save(const std::filesystem::path& filename) override;

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

         virtual MMG5_pSol& getHandle() override;
         virtual const MMG5_pSol& getHandle() const override;

      private:
         std::reference_wrapper<Mesh2D> m_mesh;
         MMG5_pSol m_sol;
  };

   /**
    * @brief A vector solution which does not have a mesh assigned to it.
    *
    * To unlock the full functionality of the class you must call the
    * @ref setMesh(Mesh2D&) method. For example, when loading it from file:
    *
    * @code{.cpp}
    * auto sol = VectorSolution2D::load(filename).setMesh(mesh);
    * @endcode
    */
   class IncompleteVectorSolution2D
   {
      public:
         /**
          * @brief Constructs an empty Vector solution object without a mesh.
          */
         IncompleteVectorSolution2D();

         /**
          * @brief Constructs a Vector solution with `n` unitialized entries.
          * @param[in] n Number of entries that the solution has.
          */
         IncompleteVectorSolution2D(int n);

         /**
          * @brief Frees the data if it still owns the data, i.e. the
          * setMesh(Mesh2D&) method has not been called.
          */
         virtual ~IncompleteVectorSolution2D();

         /**
          * @brief Sets the associated mesh and moves ownership to the new
          * object.
          *
          * @param[in] mesh Reference to mesh.
          *
          * @returns An object of type VectorSolution2D which represents the
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
         VectorSolution2D setMesh(Mesh2D& mesh);

         double* begin()
         {
            return getHandle()->m + getHandle()->size;
         }

         double* end()
         {
            return getHandle()->m + getHandle()->size * (getHandle()->np + 1);
         }

         MMG5_pSol& getHandle();
         const MMG5_pSol& getHandle() const;

      private:
         MMG5_pSol m_sol;
         bool m_isOwner;
   };
}

#endif
