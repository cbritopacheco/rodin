/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <string>
#include <cstdio>
#include <cstring>

#include <mmg2d/mmg2d.h>
#include <libmmgcommon.h>
#include <common/mmgcommon.h>

#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
   // ---- ScalarSolution2D --------------------------------------------------
   ScalarSolution2D::ScalarSolution2D(Mesh2D& mesh)
      : m_mesh(mesh)
   {
      auto calloc =
         [this]()
         {
            // Returns false on fail
            MMG5_SAFE_CALLOC(m_sol, 1, MMG5_Sol, return false);
            // Returns true on success
            return true;
         };

      if (!calloc())
         Alert::Exception("Failed to allocate memory for the mesh").raise();

      m_sol->dim  = 2; // Supported on 2D mesh
      m_sol->ver  = 2;
      m_sol->size = 1; // Scalar solution
      m_sol->type = MMG5_Scalar;
   }

   ScalarSolution2D::~ScalarSolution2D()
   {
      if (m_sol)
      {
         if (m_sol->m)
            MMG5_SAFE_FREE(m_sol->m);
         if (m_sol->namein)
            MMG5_SAFE_FREE(m_sol->namein);
         if (m_sol->nameout)
            MMG5_SAFE_FREE(m_sol->nameout);

         MMG5_SAFE_FREE(m_sol);
      }
   }

   ScalarSolution2D::ScalarSolution2D(ScalarSolution2D&& other)
      :  m_mesh(other.m_mesh),
         m_sol(other.m_sol)
   {
      other.m_sol = nullptr;
   }

   ScalarSolution2D& ScalarSolution2D::operator=(ScalarSolution2D&& other)
   {
      m_mesh = other.m_mesh;
      m_sol = other.m_sol;
      other.m_sol = nullptr;
      return *this;
   }

   IncompleteScalarSolution2D ScalarSolution2D::load(
         const std::filesystem::path& filename)
   {
      IncompleteScalarSolution2D res;
      MMG5_pSol sol = res.getHandle();

      /*
       * To load the solution file, we use basically the same methodology in
       * MMG2D_loadSol in mmg2d/inout_2d.c
       *
       * We cannot use the actual function since it requires an MMG5_pMesh
       * object which is not used in the call and subcalls of this method,
       * except for keeping track of some maximum memory constraints and
       * propagating the level of verbosity. As such, we have reimplemented the
       * loading utilizing bits of the same code.
       */
      FILE *inm;
      long posnp;
      int iswp, ier, dim, meshDim, ver, bin, np, nsols;
      int *type;

      // Read header
      meshDim = 2;

      // Verbosity is high when in Debug
      ier =  MMG5_loadSolHeader(
            filename.c_str(), meshDim,
            &inm, &ver, &bin, &iswp, &np, &dim, &nsols, &type, &posnp,
            VERBOSITY_LEVEL);

      switch (ier)
      {
         case -1:
            fclose(inm);
            MMG5_SAFE_FREE(type);
            Alert::Exception(
                  "Failed to load solution. Invalid data.").raise();
            break;
         case 0:
            fclose(inm);
            if (type)
               MMG5_SAFE_FREE(type);
            Alert::Exception(
                  "Failed to load solution. File not found: " + filename.string()).raise();
            break;
         case 1:
            // Success
            break;
         default:
            fclose(inm);
            MMG5_SAFE_FREE(type);
            Alert::Exception(
                  "Failed to load solution. Invalid error code returned.").raise();
      }

      if (nsols != 1)
      {
         fclose(inm);
         MMG5_SAFE_FREE(type);
         Alert::Exception(
               "Failed to load solution. Multiple solutions not supported.").raise();
      }

      sol->type = type[0];
      switch (type[0])
      {
         case MMG5_Scalar:
            sol->size = 1;
            break;
         case MMG5_Vector:
            sol->size = 2;
            Alert::Exception(
                  "Failed to load solution. Expected sol->size == 1, got 2.").raise();
            break;
         case MMG5_Tensor:
            sol->size = 3;
            Alert::Exception(
                  "Failed to load solution. Expected sol->size == 1, got 3.").raise();
            break;
         default:
            Alert::Exception(
                  "Failed to load solution. Unknown solution type.").raise();
      }

      assert(!sol->m);
      if (np)
      {
         sol->np  = np;
         sol->npi = np;
         sol->npmax = std::max(static_cast<int>(1.5 * sol->np), MMG2D_NPMAX);
         MMG5_SAFE_CALLOC(
               sol->m, (sol->size * (sol->npmax + 1)), double, /* No op */);
      }

      MMG5_SAFE_FREE(type);

      // Read the solutions
      rewind(inm);
      fseek(inm, posnp, SEEK_SET);

      assert(sol->ver > 1);
      assert(!bin);
      for (int k = 1; k <= sol->np; k++)
      {
        double dbuf;
        for (int i = 0; i < sol->size; i++)
        {
           auto read = [&inm, &dbuf] () { MMG_FSCANF(inm, "%lf", &dbuf); return 1; };
           if (read() < 0)
              Alert::Exception("Failed to load mesh. Error while reading.").raise();
            sol->m[sol->size * k + i] = dbuf;
        }
      }

      fclose(inm);

      return res;
   }

   void ScalarSolution2D::save(const std::filesystem::path& filename)
   {
      if (!m_sol->np || !m_sol->m)
      {
         Alert::Exception(
               "Failed to write to scalar solution to file. No data!").raise();
      }

      if (!MMG2D_saveSol(m_mesh.get().getHandle(), m_sol, filename.c_str()))
      {
         Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
      }
   }

   ScalarSolution2D::ScalarSolution2D(const ScalarSolution2D& other)
      : m_mesh(other.m_mesh)
   {
      assert(other.m_sol);

      // Allocate memory for the m_sol object
      auto calloc =
         [this]()
         {
            // Returns false on fail
            MMG5_SAFE_CALLOC(m_sol, 1, MMG5_Sol, return false);
            // Returns true on success
            return true;
         };

      if (!calloc())
         Alert::Exception("Failed to allocate memory for the mesh").raise();

      // Copy the fields
      m_sol->dim  = 2; // Supported on 2D mesh
      m_sol->ver  = 2; // Version 2
      m_sol->type = MMG5_Scalar;
      m_sol->size = 1; // One solution per entity

      m_sol->np = other.m_sol->np;
      m_sol->npmax = other.m_sol->npmax;
      m_sol->npi = other.m_sol->npi;
      m_sol->entities = other.m_sol->entities;

      if (m_sol->np)
      {
         /*
          * We should be keeping track of the memory usage in the mesh object
          * that would be associated with this solution. However, we don't do
          * that since it would require initiliazing the object with a Mesh2D.
          * For now we just assume that they're independent.
          * MMG5_ADD_MEM(
          *      mesh, (m_sol->size * (m_sol->npmax + 1)) * sizeof(double),"", ;);
          */
         MMG5_SAFE_CALLOC(m_sol->m, m_sol->npmax + 1, double, /* No op */);
         std::memcpy(m_sol->m, other.m_sol->m, m_sol->size * (m_sol->npmax + 1));
      }

      m_sol->umin = other.m_sol->umin;
      m_sol->umax = other.m_sol->umax;

      if (other.m_sol->namein)
      {
         auto nameInLength = std::strlen(other.m_sol->namein);
         MMG5_SAFE_CALLOC(m_sol->namein, nameInLength, char, /* No op */);
         std::memcpy(m_sol->namein, other.m_sol->namein, nameInLength + 1);
      }

      if (other.m_sol->nameout)
      {
         auto nameOutLength = std::strlen(other.m_sol->namein);
         MMG5_SAFE_CALLOC(m_sol->nameout, nameOutLength, char, /* No op */);
         std::memcpy(m_sol->namein, other.m_sol->namein, nameOutLength + 1);
      }
   }

   ScalarSolution2D&
   ScalarSolution2D::operator=(const ScalarSolution2D& other)
   {
      if (this != &other)
      {
         ScalarSolution2D tmp(other);

         std::swap(m_sol->dim, tmp.m_sol->dim);
         std::swap(m_sol->entities, tmp.m_sol->entities);
         std::swap(m_sol->m, tmp.m_sol->m);
         std::swap(m_sol->namein, tmp.m_sol->namein);
         std::swap(m_sol->nameout, tmp.m_sol->nameout);
         std::swap(m_sol->np, tmp.m_sol->np);
         std::swap(m_sol->npi, tmp.m_sol->npi);
         std::swap(m_sol->npmax, tmp.m_sol->npmax);
         std::swap(m_sol->size, tmp.m_sol->size);
         std::swap(m_sol->type, tmp.m_sol->type);
         std::swap(m_sol->umax, tmp.m_sol->umax);
         std::swap(m_sol->umin, tmp.m_sol->umin);
         std::swap(m_sol->ver, tmp.m_sol->ver);
      }

      return *this;
   }

   ScalarSolution2D&
   ScalarSolution2D::setMesh(Mesh2D& mesh)
   {
      if (mesh.count() != count())
         Alert::Exception("mesh.count() != count()").raise();
      m_mesh = mesh;
      return *this;
   }

   const Mesh2D& ScalarSolution2D::getMesh() const
   {
      return m_mesh;
   }

   Mesh2D& ScalarSolution2D::getMesh()
   {
      return m_mesh;
   }

   ScalarSolution2D::Iterator ScalarSolution2D::begin()
   {
      return Iterator(&m_sol->m[1]);
   }

   ScalarSolution2D::Iterator ScalarSolution2D::end()
   {
      return Iterator(&m_sol->m[count() + 1]);
   }

   ScalarSolution2D::ConstIterator ScalarSolution2D::begin() const
   {
      return ConstIterator(&m_sol->m[1]);
   }

   ScalarSolution2D::ConstIterator ScalarSolution2D::end() const
   {
      return ConstIterator(&m_sol->m[count() + 1]);
   }

   MMG5_pSol& ScalarSolution2D::getHandle()
   {
      return m_sol;
   }

   const MMG5_pSol& ScalarSolution2D::getHandle() const
   {
      return m_sol;
   }

   // ---- IncompleteScalarSolution2D -------------------------------------------
   IncompleteScalarSolution2D::IncompleteScalarSolution2D()
      : m_isOwner(true)
   {
      auto calloc =
         [this]()
         {
            // Returns false on fail
            MMG5_SAFE_CALLOC(m_sol, 1, MMG5_Sol, return false);
            // Returns true on success
            return true;
         };

      if (!calloc())
         Alert::Exception("Failed to allocate memory for the mesh").raise();

      m_sol->dim  = 2; // Supported on 2D mesh
      m_sol->ver  = 2;
      m_sol->size = 1; // Scalar solution
      m_sol->type = MMG5_Scalar;
   }

   IncompleteScalarSolution2D::IncompleteScalarSolution2D(int size)
      : IncompleteScalarSolution2D()
   {
      if (size)
      {
         m_sol->np  = size;
         m_sol->npi = size;
         m_sol->npmax = std::max(static_cast<int>(1.5 * m_sol->np), MMG2D_NPMAX);
         MMG5_SAFE_CALLOC(
               m_sol->m, (m_sol->size * (m_sol->npmax + 1)), double, /* No op */);
      }
   }

   IncompleteScalarSolution2D::~IncompleteScalarSolution2D()
   {
      if (m_isOwner)
      {
         if (m_sol)
         {
            if (m_sol->m)
               MMG5_SAFE_FREE(m_sol->m);
            if (m_sol->namein)
               MMG5_SAFE_FREE(m_sol->namein);
            if (m_sol->nameout)
               MMG5_SAFE_FREE(m_sol->nameout);

            MMG5_SAFE_FREE(m_sol);
         }
      }
   }

   ScalarSolution2D IncompleteScalarSolution2D::setMesh(Mesh2D& mesh)
   {
      ScalarSolution2D res(mesh);
      res.getHandle() = m_sol;
      m_isOwner = false;
      return res;
   }

   MMG5_pSol& IncompleteScalarSolution2D::getHandle()
   {
      return m_sol;
   }

   const MMG5_pSol& IncompleteScalarSolution2D::getHandle() const
   {
      return m_sol;
   }

   // ---- Iterator implementation -------------------------------------------
   ScalarSolution2D::Iterator::Iterator(pointer ptr)
      : m_ptr(ptr)
   {}

   ScalarSolution2D::Iterator::reference
   ScalarSolution2D::Iterator::operator*()
   const
   {
      return *m_ptr;
   }

   ScalarSolution2D::Iterator::pointer
   ScalarSolution2D::Iterator::operator->()
   {
      return m_ptr;
   }

   ScalarSolution2D::Iterator&
   ScalarSolution2D::Iterator::operator++()
   {
      m_ptr++;
      return *this;
   }

   ScalarSolution2D::Iterator
   ScalarSolution2D::Iterator::operator++(int)
   {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
   }

   bool operator==(
         const ScalarSolution2D::Iterator& a,
         const ScalarSolution2D::Iterator& b)
   {
      return a.m_ptr == b.m_ptr;
   };

   bool operator!=(
         const ScalarSolution2D::Iterator& a,
         const ScalarSolution2D::Iterator& b)
   {
      return a.m_ptr != b.m_ptr;
   };

   // ---- ConstIterator implementation --------------------------------------
   ScalarSolution2D::ConstIterator::ConstIterator(pointer ptr)
      : m_ptr(ptr)
   {}

   ScalarSolution2D::ConstIterator::const_reference
   ScalarSolution2D::ConstIterator::operator*() const
   {
      return *m_ptr;
   }

   ScalarSolution2D::ConstIterator::pointer
   ScalarSolution2D::ConstIterator::operator->()
   {
      return m_ptr;
   }

   ScalarSolution2D::ConstIterator&
   ScalarSolution2D::ConstIterator::operator++()
   {
      m_ptr++;
      return *this;
   }

   ScalarSolution2D::ConstIterator
   ScalarSolution2D::ConstIterator::operator++(int)
   {
      ConstIterator tmp = *this;
      ++(*this);
      return tmp;
   }

   bool
   operator==(
         const ScalarSolution2D::ConstIterator& a,
         const ScalarSolution2D::ConstIterator& b)
   {
      return a.m_ptr == b.m_ptr;
   };

   bool
   operator!=(
         const ScalarSolution2D::ConstIterator& a,
         const ScalarSolution2D::ConstIterator& b)
   {
      return a.m_ptr != b.m_ptr;
   };
}
