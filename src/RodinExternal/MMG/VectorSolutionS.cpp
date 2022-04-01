#include <cstring>

#include "Rodin/Alert.h"
#include "MeshS.h"

#include "VectorSolutionS.h"

namespace Rodin::External::MMG
{
  // ---- VectorSolutionS --------------------------------------------------
  VectorSolutionS::VectorSolutionS(MeshS& mesh)
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

     m_sol->dim  = mesh.getDimension(); // Supported on S mesh
     m_sol->ver  = 2;
     m_sol->size = 2; // Vector solution
     m_sol->type = MMG5_Vector;
  }

  VectorSolutionS::~VectorSolutionS()
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

  VectorSolutionS::VectorSolutionS(VectorSolutionS&& other)
     :  m_mesh(other.m_mesh),
        m_sol(other.m_sol)
  {
     other.m_sol = nullptr;
  }

  VectorSolutionS& VectorSolutionS::operator=(VectorSolutionS&& other)
  {
     m_mesh = other.m_mesh;
     m_sol = other.m_sol;
     other.m_sol = nullptr;
     return *this;
  }

  VectorSolutionS::VectorSolutionS(const VectorSolutionS& other)
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
     m_sol->dim  = other.getMesh().getDimension(); // Supported on S mesh
     m_sol->ver  = 2; // Version 2
     m_sol->type = MMG5_Vector;
     m_sol->size = 2; // Two solutions per entity

     m_sol->np = other.m_sol->np;
     m_sol->npmax = other.m_sol->npmax;
     m_sol->npi = other.m_sol->npi;
     m_sol->entities = other.m_sol->entities;

     if (m_sol->np)
     {
        /*
         * We should be keeping track of the memory usage in the mesh object
         * that would be associated with this solution. However, we don't do
         * that since it would require initiliazing the object with a MeshS.
         * For now we just assume that they're independent.
         * MMG5_ADD_MEM(
         *      mesh, (m_sol->size * (m_sol->npmax + 1)) * sizeof(double),"", ;);
         */
        MMG5_SAFE_CALLOC(m_sol->m, m_sol->size * (m_sol->npmax + 1), double, /* No op */);
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

  VectorSolutionS&
  VectorSolutionS::operator=(const VectorSolutionS& other)
  {
     if (this != &other)
     {
        VectorSolutionS tmp(other);

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

  IncompleteVectorSolutionS VectorSolutionS::load(
        const boost::filesystem::path& filename)
  {
     IncompleteVectorSolutionS res;
     MMG5_pSol sol = res.getHandle();

     /*
      * To load the solution file, we use basically the same methodology in
      * MMGS_loadSol in mmgS/inout_S.c
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
           Alert::Exception(
                 "Failed to load solution. Expected sol->size == 2, got 1.").raise();
           break;
        case MMG5_Vector:
           sol->size = 2;
           break;
        case MMG5_Tensor:
           sol->size = 3;
           Alert::Exception(
                 "Failed to load solution. Expected sol->size == 2, got 3.").raise();
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
        sol->npmax = std::max(static_cast<int>(1.5 * sol->np), MMGS_NPMAX);
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

  void VectorSolutionS::save(const boost::filesystem::path& filename)
  {
     if (!m_sol->np || !m_sol->m)
     {
        Alert::Exception(
              "Failed to write to Vector solution to file. No data!").raise();
     }

     if (!MMGS_saveSol(m_mesh.get().getHandle(), m_sol, filename.c_str()))
     {
        Alert::Exception("Failed to open file for writing: " + filename.string()).raise();
     }
  }

  VectorSolutionS& VectorSolutionS::setMesh(MeshS& mesh)
  {
     assert(mesh.count(MeshS::Entity::Vertex) == count());
     m_mesh = mesh;
     return *this;
  }

  const MeshS& VectorSolutionS::getMesh() const
  {
     return m_mesh;
  }

  MeshS& VectorSolutionS::getMesh()
  {
     return m_mesh;
  }

  MMG5_pSol& VectorSolutionS::getHandle()
  {
     return m_sol;
  }

  const MMG5_pSol& VectorSolutionS::getHandle() const
  {
     return m_sol;
  }

  // ---- IncompleteVectorSolutionS ----------------------------------------
  IncompleteVectorSolutionS::IncompleteVectorSolutionS()
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

     m_sol->dim  = 3; // Supported on S mesh
     m_sol->ver  = 2;
     m_sol->size = 2; // Vector solution
     m_sol->type = MMG5_Vector;
  }

  IncompleteVectorSolutionS::IncompleteVectorSolutionS(int n)
     : IncompleteVectorSolutionS()
  {
     if (n)
     {
        m_sol->np  = n;
        m_sol->npi = n;
        m_sol->npmax = std::max(static_cast<int>(1.5 * m_sol->np), MMGS_NPMAX);
        MMG5_SAFE_CALLOC(
              m_sol->m, (m_sol->size * (m_sol->npmax + 1)), double, /* No op */);
     }
  }

  IncompleteVectorSolutionS::~IncompleteVectorSolutionS()
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

  VectorSolutionS IncompleteVectorSolutionS::setMesh(MeshS& mesh)
  {
     VectorSolutionS res(mesh);
     res.getHandle() = m_sol;
     m_isOwner = false;
     return res;
  }

  MMG5_pSol& IncompleteVectorSolutionS::getHandle()
  {
     return m_sol;
  }

  const MMG5_pSol& IncompleteVectorSolutionS::getHandle() const
  {
     return m_sol;
  }
}
