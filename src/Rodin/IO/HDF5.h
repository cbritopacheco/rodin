/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HDF5_H
#define RODIN_IO_HDF5_H

#include <vector>
#include <cstdint>
#include <type_traits>
#include <boost/filesystem/path.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Geometry/Polytope.h"

#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
  #include <hdf5.h>
#endif

namespace Rodin::IO
{
  class HDF5
  {
    public:
      template <class MeshType>
      static void saveMesh(const MeshType& mesh, const boost::filesystem::path& filename)
      {
#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
        const auto file = H5Fcreate(
            filename.c_str(),
            H5F_ACC_TRUNC,
            H5P_DEFAULT,
            H5P_DEFAULT);
        check(file >= 0, "Failed to create HDF5 file.");

        const auto group = H5Gcreate2(file, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(group >= 0, "Failed to create /Mesh group.");

        const size_t sdim = mesh.getSpaceDimension();
        std::vector<double> vertices(mesh.getVertexCount() * sdim);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const size_t i = static_cast<size_t>(it->getIndex());
          for (size_t d = 0; d < sdim; ++d)
            vertices[i * sdim + d] = (*it)(d);
        }

        const hsize_t vdims[2] = {
          static_cast<hsize_t>(mesh.getVertexCount()),
          static_cast<hsize_t>(sdim)
        };
        writeDataset(file, "/Mesh/Vertices", H5T_NATIVE_DOUBLE, vdims, 2, vertices.data());

        const size_t cellCount = mesh.getCellCount();
        std::vector<unsigned long long> connectivity;
        std::vector<unsigned long long> offsets(cellCount + 1, 0);
        std::vector<int> types(cellCount, 0);
        std::vector<int> xdmfTopology;
        size_t c = 0;
        for (auto it = mesh.getCell(); !it.end(); ++it, ++c)
        {
          const auto& verticesKey = it->getVertices();
          const int xdmfType = getXDMFTopologyType(it->getGeometry());
          check(xdmfType >= 0, "Unsupported cell geometry for HDF5/XDMF export.");

          offsets[c] = static_cast<unsigned long long>(connectivity.size());
          types[c] = static_cast<int>(it->getGeometry());
          xdmfTopology.push_back(xdmfType);
          for (size_t j = 0; j < verticesKey.size(); ++j)
          {
            const auto v = static_cast<unsigned long long>(verticesKey[j]);
            connectivity.push_back(v);
            xdmfTopology.push_back(static_cast<int>(v));
          }
        }
        offsets[cellCount] = static_cast<unsigned long long>(connectivity.size());

        const hsize_t cdims[1] = { static_cast<hsize_t>(connectivity.size()) };
        writeDataset(file, "/Mesh/CellConnectivity", H5T_NATIVE_ULLONG, cdims, 1, connectivity.data());

        const hsize_t odims[1] = { static_cast<hsize_t>(offsets.size()) };
        writeDataset(file, "/Mesh/CellOffsets", H5T_NATIVE_ULLONG, odims, 1, offsets.data());

        const hsize_t tdims[1] = { static_cast<hsize_t>(types.size()) };
        writeDataset(file, "/Mesh/CellTypes", H5T_NATIVE_INT, tdims, 1, types.data());

        const hsize_t xdims[1] = { static_cast<hsize_t>(xdmfTopology.size()) };
        writeDataset(file, "/Mesh/XDMFTopology", H5T_NATIVE_INT, xdims, 1, xdmfTopology.data());

        H5Gclose(group);
        H5Fclose(file);
#else
        (void)mesh;
        (void)filename;
        Alert::Exception() << "Rodin was built without HDF5 support." << Alert::Raise;
#endif
      }

      template <class GridFunctionType>
      static void saveGridFunction(
          const GridFunctionType& gf,
          const boost::filesystem::path& filename)
      {
#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
        const auto file = H5Fcreate(
            filename.c_str(),
            H5F_ACC_TRUNC,
            H5P_DEFAULT,
            H5P_DEFAULT);
        check(file >= 0, "Failed to create HDF5 file.");

        const auto group = H5Gcreate2(file, "/GridFunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(group >= 0, "Failed to create /GridFunction group.");

        const auto& data = gf.getData();
        std::vector<double> values(static_cast<size_t>(data.size()));
        for (size_t i = 0; i < values.size(); ++i)
          values[i] = static_cast<double>(data[i]);

        const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        writeDataset(file, "/GridFunction/Values", H5T_NATIVE_DOUBLE, dims, 1, values.data());

        const unsigned long long dofCount = static_cast<unsigned long long>(gf.getSize());
        const unsigned long long vectorDim = static_cast<unsigned long long>(gf.getDimension());
        writeScalarAttribute(file, "/GridFunction/Values", "DofCount", H5T_NATIVE_ULLONG, &dofCount);
        writeScalarAttribute(file, "/GridFunction/Values", "VectorDimension", H5T_NATIVE_ULLONG, &vectorDim);

        H5Gclose(group);
        H5Fclose(file);
#else
        (void)gf;
        (void)filename;
        Alert::Exception() << "Rodin was built without HDF5 support." << Alert::Raise;
#endif
      }

    private:
      static int getXDMFTopologyType(Geometry::Polytope::Type g)
      {
        switch (g)
        {
          case Geometry::Polytope::Type::Segment:
            return 2;
          case Geometry::Polytope::Type::Triangle:
            return 4;
          case Geometry::Polytope::Type::Quadrilateral:
            return 5;
          case Geometry::Polytope::Type::Tetrahedron:
            return 6;
          case Geometry::Polytope::Type::Wedge:
            return 8;
          case Geometry::Polytope::Type::Hexahedron:
            return 9;
          default:
            return -1;
        }
      }

#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
      static void check(bool condition, const char* msg)
      {
        if (!condition)
          Alert::Exception() << msg << Alert::Raise;
      }

      static void writeDataset(
          hid_t file,
          const char* path,
          hid_t dtype,
          const hsize_t* dims,
          int rank,
          const void* data)
      {
        const auto space = H5Screate_simple(rank, dims, nullptr);
        check(space >= 0, "Failed to create HDF5 dataspace.");
        const auto dset = H5Dcreate2(file, path, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(dset >= 0, "Failed to create HDF5 dataset.");
        const auto status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        check(status >= 0, "Failed to write HDF5 dataset.");
        H5Dclose(dset);
        H5Sclose(space);
      }

      static void writeScalarAttribute(
          hid_t file,
          const char* datasetPath,
          const char* attributeName,
          hid_t dtype,
          const void* value)
      {
        const auto dset = H5Dopen2(file, datasetPath, H5P_DEFAULT);
        check(dset >= 0, "Failed to open HDF5 dataset for writing attribute.");
        const auto space = H5Screate(H5S_SCALAR);
        check(space >= 0, "Failed to create HDF5 scalar dataspace.");
        const auto attr = H5Acreate2(dset, attributeName, dtype, space, H5P_DEFAULT, H5P_DEFAULT);
        check(attr >= 0, "Failed to create HDF5 attribute.");
        const auto status = H5Awrite(attr, dtype, value);
        check(status >= 0, "Failed to write HDF5 attribute.");
        H5Aclose(attr);
        H5Sclose(space);
        H5Dclose(dset);
      }
#endif
  };
}

#endif
