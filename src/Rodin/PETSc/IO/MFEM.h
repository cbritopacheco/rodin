/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_MFEM_H
#define RODIN_PETSC_IO_MFEM_H

#include <petscvec.h>

#include <iomanip>
#include <limits>
#include <type_traits>
#include <vector>
#include <algorithm>

#include "Rodin/IO/MFEM.h"                // for MFEM::Ordering, Vandermonde*, Nodes*
#include "Rodin/IO/GridFunctionPrinter.h" // primary template + base templates
#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/ForwardDecls.h"

namespace Rodin::IO
{
  // --------------------------------------------------------------------------
  // P0 MFEM printer for PETSc Vec
  // --------------------------------------------------------------------------
  template <class Range, class Ctx>
  class GridFunctionPrinter<
      FileFormat::MFEM,
      Variational::P0<Range, Geometry::Mesh<Ctx>>,
      ::Vec> final
    : public GridFunctionPrinterBase<
        FileFormat::MFEM,
        Variational::P0<Range, Geometry::Mesh<Ctx>>,
        ::Vec>
  {
    public:
      using FESType    = Variational::P0<Range, Geometry::Mesh<Ctx>>;
      using DataType   = ::Vec;
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      using Parent     = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      static_assert(std::is_same_v<ScalarType, PetscScalar>,
        "PETSc MFEM printer: Range scalar must be PetscScalar.");

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void printData(std::ostream& os) override
      {
        os << std::setprecision(std::numeric_limits<PetscReal>::max_digits10);

        const auto& gf   = this->getObject();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        // MFEM P0: one scalar DOF per element per component.
        // Base prints Ordering: Nodes, but for P0 that’s effectively "block by component":
        // emit all scalar dofs for comp 0, then comp 1, ...
        const size_t vdim = fes.getVectorDimension();
        const size_t D    = mesh.getDimension();

        const size_t totalSize  = fes.getSize();
        const size_t scalarSize = totalSize / vdim;

        // We assume Rodin stores component-blocked: [u0...uN | v0...vN | ...]
        // So MFEM "Nodes" layout matches printing by component blocks.
        if constexpr (std::is_same_v<Ctx, Context::MPI>)
        {
          gf.acquire();
        }

        for (size_t c = 0; c < vdim; ++c)
        {
          const Index base = static_cast<Index>(c * scalarSize);
          for (Index i = 0; i < static_cast<Index>(scalarSize); ++i)
          {
            // For P0, scalar dof index corresponds to cell index ordering
            const Index dof = base + i;
            os << gf[dof] << '\n';
          }
        }

        if constexpr (std::is_same_v<Ctx, Context::MPI>)
        {
          gf.flush();
        }
      }
  };

  // --------------------------------------------------------------------------
  // P1 MFEM printer for PETSc Vec
  // --------------------------------------------------------------------------
  template <class Range, class Ctx>
  class GridFunctionPrinter<
      FileFormat::MFEM,
      Variational::P1<Range, Geometry::Mesh<Ctx>>,
      ::Vec> final
    : public GridFunctionPrinterBase<
        FileFormat::MFEM,
        Variational::P1<Range, Geometry::Mesh<Ctx>>,
        ::Vec>
  {
    public:
      using FESType    = Variational::P1<Range, Geometry::Mesh<Ctx>>;
      using DataType   = ::Vec;
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      using Parent     = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      static_assert(std::is_same_v<ScalarType, PetscScalar>,
        "PETSc MFEM printer: Range scalar must be PetscScalar.");

      using FESMeshContextType = typename FormLanguage::Traits<typename FESType::MeshType>::ContextType;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void printData(std::ostream& os) override
      {
        const auto& gf = this->getObject();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          const size_t sz = gf.getSize();
          for (size_t i = 0; i < sz; ++i)
            os << gf[i] << '\n';
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& fes   = gf.getFiniteElementSpace();
          const auto& shard = fes.getShard();

          const size_t localSz = shard.getSize();

          // Force read-acquire, not write-acquire:
          const auto& cgf = static_cast<const std::decay_t<decltype(gf)>&>(gf);
          cgf.acquire();

          for (size_t i = 0; i < localSz; ++i)
          {
            const Index g = fes.getGlobalIndex(i);   // IMPORTANT: shard-local dof i -> global dof
            os << cgf[g] << '\n';
          }

          cgf.flush();
        }
        else
        {
          assert(false);
        }
      }
  };

  // --------------------------------------------------------------------------
  // H1<K> MFEM printer for PETSc Vec (Context::Local and Context::MPI)
  // Matches the Math::Vector version, but reads coefficients from Vec.
  // Base prints Ordering: VectorDimension (1) for H1, so we emit components
  // per MFEM scalar node in MFEM node order.
  // --------------------------------------------------------------------------
  template <size_t K, class Range, class Ctx>
  class GridFunctionPrinter<
      FileFormat::MFEM,
      Variational::H1<K, Range, Geometry::Mesh<Ctx>>,
      ::Vec> final
    : public GridFunctionPrinterBase<
        FileFormat::MFEM,
        Variational::H1<K, Range, Geometry::Mesh<Ctx>>,
        ::Vec>
  {
    public:
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Ctx>>;
      using DataType   = ::Vec;
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      using Parent     = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      static_assert(std::is_same_v<ScalarType, PetscScalar>,
        "PETSc MFEM printer: Range scalar must be PetscScalar.");

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void printData(std::ostream& os) override
      {
        os << std::setprecision(std::numeric_limits<PetscReal>::max_digits10);

        const auto& gf   = this->getObject();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        const size_t vdim       = fes.getVectorDimension();
        const size_t D          = mesh.getDimension();
        const size_t totalSize  = fes.getSize();
        const size_t scalarSize = totalSize / vdim;

        std::vector<uint8_t> written(scalarSize, false);

        const auto to_scalar_dof = [&](Index dof) -> Index
        {
          if (vdim > 1 && static_cast<size_t>(dof) >= scalarSize)
            return dof % static_cast<Index>(scalarSize);
          return dof;
        };

        // Read a Rodin DOF (already component-blocked index) from PETSc Vec
        const auto read = [&](Index dof) -> PetscScalar
        {
          return gf[dof];
        };

        // Emit in MFEM Ordering::VectorDimension: for each MFEM scalar node, write all components
        const auto emit_scalar_dof = [&](Index rodin_scalar_dof)
        {
          const Index sdof = to_scalar_dof(rodin_scalar_dof);
          const size_t s = static_cast<size_t>(sdof);
          if (s >= scalarSize || written[s])
            return;

          for (size_t c = 0; c < vdim; ++c)
          {
            const Index dof = sdof + static_cast<Index>(c * scalarSize);
            os << read(dof) << '\n';
          }
          written[s] = true;
        };

        // Change-of-nodes matrices (Rodin Fekete -> MFEM nodes) for triangle/tet.
        // We only need them to output MFEM face/cell interior values for simplex elements,
        // like your Math::Vector implementation.
        const auto& tri_change = []() -> const Math::Matrix<PetscScalar>&
        {
          static thread_local Math::Matrix<PetscScalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<PetscScalar>();
          }
          return C;
        }();

        const auto& tet_change = []() -> const Math::Matrix<PetscScalar>&
        {
          static thread_local Math::Matrix<PetscScalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<PetscScalar>();
          }
          return C;
        }();

        if constexpr (std::is_same_v<Ctx, Context::MPI>)
        {
          // Ensure ghosted values are available if your GridFunction uses ghosts.
          gf.acquire();
        }

        // --------------------------------------------------------------------
        // 1) Vertices
        // --------------------------------------------------------------------
        const size_t nVertices = mesh.getConnectivity().getCount(0);
        std::vector<Index> vertexScalarDof(nVertices);
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const auto& vdofs = fes.getDOFs(0, v);
          vertexScalarDof[static_cast<size_t>(v)] = to_scalar_dof(vdofs(0));
        }
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
          emit_scalar_dof(vertexScalarDof[static_cast<size_t>(v)]);

        // --------------------------------------------------------------------
        // 2) Edges: oriented vmin -> vmax, emit interior edge nodes in MFEM order
        // --------------------------------------------------------------------
        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const size_t nEdges = mesh.getConnectivity().getCount(1);

          std::vector<Index> interior;
          for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
          {
            const auto& edgeVerts = conn10[e];
            const Index v0 = edgeVerts[0];
            const Index v1 = edgeVerts[1];

            const Index vmin = std::min(v0, v1);
            const Index vmax = std::max(v0, v1);

            const Index vminDof = vertexScalarDof[static_cast<size_t>(vmin)];
            const Index vmaxDof = vertexScalarDof[static_cast<size_t>(vmax)];

            const auto& edofs = fes.getDOFs(1, e);

            interior.clear();
            for (Index k = 0; k < static_cast<Index>(edofs.size()); ++k)
            {
              const Index sd = to_scalar_dof(edofs(k));
              if (sd != vminDof && sd != vmaxDof)
                interior.push_back(sd);
            }

            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index sd : interior)
              emit_scalar_dof(sd);
          }
        }

        // --------------------------------------------------------------------
        // 3) Faces (only special handling for 3D triangle faces)
        // --------------------------------------------------------------------
        if (D == 3)
        {
          const size_t faceCount = mesh.getConnectivity().getCount(2);

          constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
          const int p  = static_cast<int>(K);
          const int nV = 3;
          const int nE = 3 * (p - 1);
          const int triInteriorOffset = nV + nE;

          Math::Vector<PetscScalar> uR(TriN);
          std::vector<Math::Vector<PetscScalar>> uM(vdim, Math::Vector<PetscScalar>(TriN));

          for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
          {
            const auto faceGeom = mesh.getGeometry(2, f);
            const auto& fdofs   = fes.getDOFs(2, f);

            if (faceGeom == Geometry::Polytope::Type::Triangle)
            {
              // Convert Rodin nodal values -> MFEM nodal values, then emit MFEM face-interior nodes.
              assert(static_cast<size_t>(fdofs.size()) == TriN);

              for (size_t c = 0; c < vdim; ++c)
              {
                for (size_t k = 0; k < TriN; ++k)
                {
                  const Index sd = to_scalar_dof(fdofs(static_cast<Index>(k)));
                  uR(static_cast<Index>(k)) = read(sd + static_cast<Index>(c * scalarSize));
                }
                uM[c] = tri_change * uR;
              }

              // Face interior nodes in MFEM ordering (after vertices+edges)
              int loc = 0;
              for (int j = 1; j < p; ++j)
              {
                for (int i = 1; i + j < p; ++i)
                {
                  const int idx = triInteriorOffset + loc++;
                  for (size_t c = 0; c < vdim; ++c)
                    os << uM[c](static_cast<Index>(idx)) << '\n';
                }
              }
            }
            else
            {
              // Fallback: emit all face DOFs by scalar-dof
              for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                emit_scalar_dof(fdofs(k));
            }
          }
        }

        // --------------------------------------------------------------------
        // 4) Element interiors (special handling for triangles / tetrahedra)
        // --------------------------------------------------------------------
        const size_t nCells = mesh.getConnectivity().getCount(D);

        if (D == 2)
        {
          constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
          const int p  = static_cast<int>(K);
          const int nV = 3;
          const int nE = 3 * (p - 1);
          const int triInteriorOffset = nV + nE;

          Math::Vector<PetscScalar> uR(TriN);
          std::vector<Math::Vector<PetscScalar>> uM(vdim, Math::Vector<PetscScalar>(TriN));

          for (Index cidx = 0; cidx < static_cast<Index>(nCells); ++cidx)
          {
            const auto geom   = mesh.getGeometry(2, cidx);
            const auto& cdofs = fes.getDOFs(2, cidx);

            if (geom == Geometry::Polytope::Type::Triangle)
            {
              assert(static_cast<size_t>(cdofs.size()) == TriN);

              for (size_t comp = 0; comp < vdim; ++comp)
              {
                for (size_t k = 0; k < TriN; ++k)
                {
                  const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                  uR(static_cast<Index>(k)) = read(sd + static_cast<Index>(comp * scalarSize));
                }
                uM[comp] = tri_change * uR;
              }

              int loc = 0;
              for (int j = 1; j < p; ++j)
              {
                for (int i = 1; i + j < p; ++i)
                {
                  const int idx = triInteriorOffset + loc++;
                  for (size_t comp = 0; comp < vdim; ++comp)
                    os << uM[comp](static_cast<Index>(idx)) << '\n';
                }
              }
            }
            else
            {
              for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                emit_scalar_dof(cdofs(k));
            }
          }
        }
        else if (D == 3)
        {
          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV = 4;
          const int nE = 6 * (p - 1);
          const int nF = 2 * (p - 1) * (p - 2);
          const int tetInteriorOffset = nV + nE + nF;

          Math::Vector<PetscScalar> uR(TetN);
          std::vector<Math::Vector<PetscScalar>> uM(vdim, Math::Vector<PetscScalar>(TetN));

          for (Index cidx = 0; cidx < static_cast<Index>(nCells); ++cidx)
          {
            const auto geom   = mesh.getGeometry(3, cidx);
            const auto& cdofs = fes.getDOFs(3, cidx);

            if (geom == Geometry::Polytope::Type::Tetrahedron)
            {
              assert(static_cast<size_t>(cdofs.size()) == TetN);

              for (size_t comp = 0; comp < vdim; ++comp)
              {
                for (size_t k = 0; k < TetN; ++k)
                {
                  const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                  uR(static_cast<Index>(k)) = read(sd + static_cast<Index>(comp * scalarSize));
                }
                uM[comp] = tet_change * uR;
              }

              for (int idx = tetInteriorOffset; idx < static_cast<int>(TetN); ++idx)
              {
                for (size_t comp = 0; comp < vdim; ++comp)
                  os << uM[comp](static_cast<Index>(idx)) << '\n';
              }
            }
            else
            {
              // Wedge/hex/etc: fallback scalar-dof traversal
              for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                emit_scalar_dof(cdofs(k));
            }
          }
        }
        else
        {
          for (Index cidx = 0; cidx < static_cast<Index>(nCells); ++cidx)
          {
            const auto& cdofs = fes.getDOFs(D, cidx);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
              emit_scalar_dof(cdofs(k));
          }
        }

        if constexpr (std::is_same_v<Ctx, Context::MPI>)
        {
          gf.flush();
        }
      }
  };

  // Convenience aliases to match your existing specializations pattern
  template <size_t K, class Range>
  using MFEM_H1_PETSc_Local_Printer =
    GridFunctionPrinter<FileFormat::MFEM, Variational::H1<K, Range, Geometry::Mesh<Context::Local>>, ::Vec>;

  template <size_t K, class Range>
  using MFEM_H1_PETSc_MPI_Printer =
    GridFunctionPrinter<FileFormat::MFEM, Variational::H1<K, Range, Geometry::Mesh<Context::MPI>>, ::Vec>;
}

#endif
