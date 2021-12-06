/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTION2D_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTION2D_H

#include <cstddef>
#include <iterator>
#include <functional>

#include <mmg/mmg2d/libmmg2d.h>

#include "Rodin/Alert.h"

#include "ForwardDecls.h"
#include "Mesh2D.h"

#include "ScalarSolution.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Scalar solution supported on a 2D mesh.
    *
    * An object of type ScalarSolution2D represents a function
    * @f[
    * f : \Omega \subset \mathbb{R}^2 \rightarrow \mathbb{R}
    * @f]
    * whose known values are given on vertices of some mesh @f$ \Omega @f$.
    */
   template <>
   class ScalarSolution2D<true>
      :  public ScalarSolution<2, ScalarSolution2D<true>>
   {
      public:
         class Iterator
         {
            public:
               using iterator_category = std::random_access_iterator_tag;
               using difference_type   = std::ptrdiff_t;
               using value_type        = double;
               using pointer           = double*;
               using reference         = double&;

               Iterator(pointer ptr);
               reference operator*() const;
               pointer operator->();
               Iterator& operator++();
               Iterator operator++(int);
               friend bool operator==(const Iterator& a, const Iterator& b);
               friend bool operator!=(const Iterator& a, const Iterator& b);

            private:
                 pointer m_ptr;
         };

         class ConstIterator
         {
            public:
               using iterator_category = std::forward_iterator_tag;
               using difference_type   = std::ptrdiff_t;
               using value_type        = double;
               using pointer           = double*;
               using reference         = double&;
               using const_reference   = const double&;

               ConstIterator(pointer ptr);

               const_reference operator*() const;
               pointer operator->();
               ConstIterator& operator++();
               ConstIterator operator++(int);
               friend bool operator==(const ConstIterator& a, const
                     ConstIterator& b);
               friend bool operator!=(const ConstIterator& a, const
                     ConstIterator& b);

            private:
                 pointer m_ptr;
         };

         /**
          * @brief Reads the solution text file.
          *
          * The file is read using MMGv2 format.
          *
          * @param[in] filename Name of file to read.
          */
         static ScalarSolution2D<false> load(const std::string& filename);

         /**
          * @brief Write the solution to a text file.
          *
          * The file is written using MMGv2 format.
          *
          * @param[in] filename Name of file to write.
          */
         void save(const std::string& filename);

         /**
          * @brief Initializes the object with no data
          *
          * @param[in] mesh Reference to the underlying mesh.
          */
         ScalarSolution2D(Mesh2D& mesh);

         /**
          * @brief Performs a move construction from the `other` solution object.
          *
          * @param[in] other Object to move.
          */
         ScalarSolution2D(ScalarSolution2D&& other);

         /**
          * @brief Performs a deep copy of the `other` solution object.
          *
          * @param[in] other Object to copy.
          */
         ScalarSolution2D(const ScalarSolution2D& other);

         /**
          * @brief Move assigns the `other` solution object to this object.
          *
          * @param[in] other Object to move.
          */
         ScalarSolution2D& operator=(ScalarSolution2D&& other);

         /**
          * @brief Copy assigns the `other` solution object to this object.
          *
          * @param[in] other Object to copy.
          */
         ScalarSolution2D& operator=(const ScalarSolution2D& other);

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
         ScalarSolution2D& setMesh(Mesh2D& mesh);

         /**
          * @brief Gets the constant reference to the underlying mesh.
          *
          * @returns Constant reference to the underlying mesh.
          */
         const Mesh2D& getMesh() const;

         /**
          * @brief Gets the number of points on which the scalar solution is
          * given.
          *
          * @returns Number of points on which the scalar solution is given.
          *
          * @note In general one should not assume that the condition
          * @code{.cpp}
          *    sol.getMesh().count() == sol.count()
          * @endcode
          * is satisfied since the underlying mesh might have been modified.
          */
         int count() const;

         Iterator begin();
         Iterator end();

         ConstIterator begin() const;
         ConstIterator end() const;

         MMG5_pSol& getHandle();
         const MMG5_pSol& getHandle() const;

      private:
         std::reference_wrapper<Mesh2D> m_mesh;
         MMG5_pSol m_sol;
   };

   /**
    * @brief A scalar solution which does not have a mesh assigned to it.
    *
    * To unlock the full functionality of the class you must call the
    * @ref setMesh(Mesh2D&) method. For example, when loading it from file:
    *
    * @code{.cpp}
    * ScalarSolution2D sol = ScalarSolution2D::load(filename).setMesh(mesh);
    * @endcode
    */
   template <>
   class ScalarSolution2D<false>
   {
      public:
         ScalarSolution2D();

         /**
          * @brief Sets the associated mesh.
          *
          * @param[in] mesh Reference to mesh.
          *
          * @returns An object of type ScalarSolution2D which represents the
          * object with all its functionality.
          *
          * @note The method does not incur any significant performance since
          * no data is copied.
          *
          * @note The method does not check to see if the mesh is compatible
          * with the current data in the solution. In general, it is up to the
          * user to ensure that the number of points are the same, keep track
          * of the modifications to the underlying mesh, etc.
          */
         ScalarSolution2D<true> setMesh(Mesh2D& mesh) const;

         MMG5_pSol& getHandle();
         const MMG5_pSol& getHandle() const;

      private:
         MMG5_pSol m_sol;
   };
}
#endif
