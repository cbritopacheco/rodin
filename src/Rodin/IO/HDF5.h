/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HDF5_H
#define RODIN_IO_HDF5_H

#include <array>
#include <string>
#include <vector>
#include <limits>

#include <boost/filesystem/path.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/AttributeIndex.h"
#include "Rodin/Geometry/PolytopeTransformationIndex.h"

#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
  #include <hdf5.h>
#endif

namespace Rodin::IO
{
  class HDF5
  {
    public:
      struct Keyword
      {
        static constexpr const char* Mesh = "/Mesh";
        static constexpr const char* SpaceDimension = "/Mesh/SpaceDimension";
        static constexpr const char* Vertices = "/Mesh/Vertices";

        static constexpr const char* Connectivity = "/Mesh/Connectivity";
        static constexpr const char* ConnectivityMaximalDimension = "/Mesh/Connectivity/MaximalDimension";
        static constexpr const char* ConnectivityCounts = "/Mesh/Connectivity/Counts";

        static constexpr const char* Attributes = "/Mesh/Attributes";

        static constexpr const char* GridFunction = "/GridFunction";
        static constexpr const char* GridFunctionValues = "/GridFunction/Values";
      };

      template <class MeshType>
      static void saveMesh(const MeshType& mesh, const boost::filesystem::path& filename)
      {
#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
        const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        check(file >= 0, "Failed to create HDF5 file.");

        const auto meshGroup = H5Gcreate2(file, Keyword::Mesh, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(meshGroup >= 0, "Failed to create /Mesh group.");
        H5Gclose(meshGroup);

        const auto attributesGroup = H5Gcreate2(file, Keyword::Attributes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(attributesGroup >= 0, "Failed to create /Mesh/Attributes group.");
        H5Gclose(attributesGroup);

        writeScalarDataset(
            file,
            Keyword::SpaceDimension,
            H5T_NATIVE_ULLONG,
            static_cast<unsigned long long>(mesh.getSpaceDimension()));

        writeVertices(file, mesh);
        saveConnectivity(file, mesh.getConnectivity());
        saveAttributes(file, mesh);

        H5Fclose(file);
#else
        (void)mesh;
        (void)filename;
        Alert::Exception() << "Rodin was built without HDF5 support." << Alert::Raise;
#endif
      }

      template <class MeshType>
      static void loadMesh(MeshType& mesh, const boost::filesystem::path& filename)
      {
#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
        const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        check(file >= 0, "Failed to open HDF5 file.");

        const auto sdim = static_cast<size_t>(
            readScalarDataset<unsigned long long>(file, Keyword::SpaceDimension, H5T_NATIVE_ULLONG));

        auto vertices = readVertices(file);
        auto connectivity = loadConnectivity<typename MeshType::ContextType>(file);
        auto attributes = loadAttributes(file, connectivity);

        Geometry::PolytopeTransformationIndex transformations;
        transformations.initialize(connectivity.getDimension());
        for (size_t d = 0; d <= connectivity.getDimension(); ++d)
          transformations.resize(d, connectivity.getCount(d));

        mesh = typename MeshType::Builder()
          .initialize(sdim)
          .setVertices(std::move(vertices))
          .setConnectivity(std::move(connectivity))
          .setAttributeIndex(std::move(attributes))
          .setTransformationIndex(std::move(transformations))
          .finalize();

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
        const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        check(file >= 0, "Failed to create HDF5 file.");

        const auto group = H5Gcreate2(file, Keyword::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(group >= 0, "Failed to create /GridFunction group.");
        H5Gclose(group);

        const auto& data = gf.getData();
        std::vector<double> values(static_cast<size_t>(data.size()));
        for (size_t i = 0; i < values.size(); ++i)
          values[i] = static_cast<double>(data[i]);

        const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        writeDataset(file, Keyword::GridFunctionValues, H5T_NATIVE_DOUBLE, dims, 1, values.data());

        const unsigned long long dofCount = static_cast<unsigned long long>(gf.getSize());
        const unsigned long long vectorDim = static_cast<unsigned long long>(gf.getDimension());
        writeScalarAttribute(file, Keyword::GridFunctionValues, "DofCount", H5T_NATIVE_ULLONG, &dofCount);
        writeScalarAttribute(file, Keyword::GridFunctionValues, "VectorDimension", H5T_NATIVE_ULLONG, &vectorDim);

        H5Fclose(file);
#else
        (void)gf;
        (void)filename;
        Alert::Exception() << "Rodin was built without HDF5 support." << Alert::Raise;
#endif
      }

    private:
#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
      static constexpr unsigned long long NullAttributeMarker = std::numeric_limits<unsigned long long>::max();

      static void check(bool condition, const char* msg)
      {
        if (!condition)
          Alert::Exception() << msg << Alert::Raise;
      }

      static bool hasPath(hid_t file, const std::string& path)
      {
        return H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
      }

      static std::string polytopePath(size_t d, const char* suffix)
      {
        return std::string(Keyword::Connectivity) + "/Polytopes/d" + std::to_string(d) + "/" + suffix;
      }

      static std::string incidencePath(size_t d, size_t dp, const char* suffix)
      {
        return std::string(Keyword::Connectivity) + "/Incidence/d" + std::to_string(d) + "_d"
             + std::to_string(dp) + "/" + suffix;
      }

      static std::string attributePath(size_t d)
      {
        return std::string(Keyword::Attributes) + "/d" + std::to_string(d);
      }

      static void writeDataset(
          hid_t file,
          const std::string& path,
          hid_t dtype,
          const hsize_t* dims,
          int rank,
          const void* data)
      {
        const auto space = H5Screate_simple(rank, dims, nullptr);
        check(space >= 0, "Failed to create HDF5 dataspace.");
        const auto dset = H5Dcreate2(file, path.c_str(), dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(dset >= 0, "Failed to create HDF5 dataset.");
        const auto status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        check(status >= 0, "Failed to write HDF5 dataset.");
        H5Dclose(dset);
        H5Sclose(space);
      }

      template <class T>
      static void writeVectorDataset(
          hid_t file,
          const std::string& path,
          hid_t dtype,
          const std::vector<T>& values)
      {
        const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        writeDataset(file, path, dtype, dims, 1, values.data());
      }

      template <class T>
      static void writeScalarDataset(hid_t file, const std::string& path, hid_t dtype, const T& value)
      {
        const hsize_t dims[1] = { 1 };
        writeDataset(file, path, dtype, dims, 1, &value);
      }

      template <class T>
      static std::vector<T> readVectorDataset(hid_t file, const std::string& path, hid_t dtype)
      {
        const auto dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
        check(dset >= 0, "Failed to open HDF5 dataset.");
        const auto space = H5Dget_space(dset);
        check(space >= 0, "Failed to open HDF5 dataspace.");
        const auto n = static_cast<size_t>(H5Sget_simple_extent_npoints(space));
        std::vector<T> values(n);
        if (n > 0)
        {
          const auto status = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
          check(status >= 0, "Failed to read HDF5 dataset.");
        }
        H5Sclose(space);
        H5Dclose(dset);
        return values;
      }

      template <class T>
      static T readScalarDataset(hid_t file, const std::string& path, hid_t dtype)
      {
        const auto values = readVectorDataset<T>(file, path, dtype);
        check(values.size() == 1, "Invalid scalar dataset shape.");
        return values[0];
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

      template <class MeshType>
      static void writeVertices(hid_t file, const MeshType& mesh)
      {
        const size_t sdim = mesh.getSpaceDimension();
        std::vector<double> vertices(mesh.getVertexCount() * sdim);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const size_t i = static_cast<size_t>(it->getIndex());
          for (size_t d = 0; d < sdim; ++d)
            vertices[i * sdim + d] = (*it)(d);
        }

        const hsize_t dims[2] = {
          static_cast<hsize_t>(mesh.getVertexCount()),
          static_cast<hsize_t>(sdim)
        };
        writeDataset(file, Keyword::Vertices, H5T_NATIVE_DOUBLE, dims, 2, vertices.data());
      }

      static Geometry::PointCloud readVertices(hid_t file)
      {
        const auto dset = H5Dopen2(file, Keyword::Vertices, H5P_DEFAULT);
        check(dset >= 0, "Failed to open /Mesh/Vertices dataset.");
        const auto space = H5Dget_space(dset);
        check(space >= 0, "Failed to open /Mesh/Vertices dataspace.");

        hsize_t dims[2] = { 0, 0 };
        const int rank = H5Sget_simple_extent_ndims(space);
        check(rank == 2, "Invalid /Mesh/Vertices dataset rank.");
        const auto readDims = H5Sget_simple_extent_dims(space, dims, nullptr);
        check(readDims == 2, "Invalid /Mesh/Vertices dataset dimensions.");

        std::vector<double> packed(static_cast<size_t>(dims[0] * dims[1]));
        if (!packed.empty())
        {
          const auto status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, packed.data());
          check(status >= 0, "Failed to read /Mesh/Vertices dataset.");
        }

        Geometry::PointCloud vertices(static_cast<std::uint8_t>(dims[1]), static_cast<size_t>(dims[0]));
        for (size_t i = 0; i < static_cast<size_t>(dims[0]); ++i)
          for (size_t d = 0; d < static_cast<size_t>(dims[1]); ++d)
            vertices(static_cast<std::uint8_t>(d), i) = packed[i * static_cast<size_t>(dims[1]) + d];

        H5Sclose(space);
        H5Dclose(dset);
        return vertices;
      }

      template <class ConnectivityType>
      static void saveConnectivity(hid_t file, const ConnectivityType& connectivity)
      {
        const auto group = H5Gcreate2(file, Keyword::Connectivity, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(group >= 0, "Failed to create /Mesh/Connectivity group.");
        H5Gclose(group);

        const size_t D = connectivity.getDimension();
        writeScalarDataset(
            file,
            Keyword::ConnectivityMaximalDimension,
            H5T_NATIVE_ULLONG,
            static_cast<unsigned long long>(D));

        std::vector<unsigned long long> counts(D + 1, 0);
        for (size_t d = 0; d <= D; ++d)
          counts[d] = static_cast<unsigned long long>(connectivity.getCount(d));
        writeVectorDataset(file, Keyword::ConnectivityCounts, H5T_NATIVE_ULLONG, counts);

        for (size_t d = 1; d <= D; ++d)
        {
          std::vector<int> geometry;
          std::vector<unsigned long long> offsets(1, 0);
          std::vector<unsigned long long> data;
          geometry.reserve(static_cast<size_t>(connectivity.getCount(d)));
          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            const auto& key = connectivity.getPolytope(d, i);
            geometry.push_back(static_cast<int>(connectivity.getGeometry(d, i)));
            for (size_t j = 0; j < key.size(); ++j)
              data.push_back(static_cast<unsigned long long>(key[j]));
            offsets.push_back(static_cast<unsigned long long>(data.size()));
          }

          writeVectorDataset(file, polytopePath(d, "Geometry"), H5T_NATIVE_INT, geometry);
          writeVectorDataset(file, polytopePath(d, "Offsets"), H5T_NATIVE_ULLONG, offsets);
          writeVectorDataset(file, polytopePath(d, "Data"), H5T_NATIVE_ULLONG, data);
        }

        for (size_t d = 0; d <= D; ++d)
        {
          for (size_t dp = 0; dp <= D; ++dp)
          {
            const auto& inc = connectivity.getIncidence(d, dp);
            if (inc.size() != connectivity.getCount(d))
              continue;
            std::vector<unsigned long long> offsets(1, 0);
            std::vector<unsigned long long> data;
            offsets.reserve(inc.size() + 1);
            for (const auto& row : inc)
            {
              for (const auto i : row)
                data.push_back(static_cast<unsigned long long>(i));
              offsets.push_back(static_cast<unsigned long long>(data.size()));
            }
            writeVectorDataset(file, incidencePath(d, dp, "Offsets"), H5T_NATIVE_ULLONG, offsets);
            writeVectorDataset(file, incidencePath(d, dp, "Data"), H5T_NATIVE_ULLONG, data);
          }
        }
      }

      template <class ContextType>
      static Geometry::Connectivity<ContextType> loadConnectivity(hid_t file)
      {
        Geometry::Connectivity<ContextType> connectivity;

        const size_t D = static_cast<size_t>(
            readScalarDataset<unsigned long long>(
                file,
                Keyword::ConnectivityMaximalDimension,
                H5T_NATIVE_ULLONG));
        connectivity.initialize(D);

        const auto counts = readVectorDataset<unsigned long long>(file, Keyword::ConnectivityCounts, H5T_NATIVE_ULLONG);
        check(counts.size() == D + 1, "Invalid /Mesh/Connectivity/Counts size.");
        connectivity.nodes(static_cast<size_t>(counts[0]));
        for (size_t d = 1; d <= D; ++d)
          connectivity.reserve(d, static_cast<size_t>(counts[d]));

        for (size_t d = 1; d <= D; ++d)
        {
          const auto geometry = readVectorDataset<int>(file, polytopePath(d, "Geometry"), H5T_NATIVE_INT);
          const auto offsets = readVectorDataset<unsigned long long>(file, polytopePath(d, "Offsets"), H5T_NATIVE_ULLONG);
          const auto data = readVectorDataset<unsigned long long>(file, polytopePath(d, "Data"), H5T_NATIVE_ULLONG);
          check(geometry.size() == static_cast<size_t>(counts[d]), "Invalid connectivity geometry dataset size.");
          check(offsets.size() == static_cast<size_t>(counts[d]) + 1, "Invalid connectivity offsets dataset size.");

          for (size_t i = 0; i < static_cast<size_t>(counts[d]); ++i)
          {
            const auto begin = static_cast<size_t>(offsets[i]);
            const auto end = static_cast<size_t>(offsets[i + 1]);
            check(end >= begin && end <= data.size(), "Invalid connectivity polytope offsets.");
            Geometry::Polytope::Key key(end - begin);
            for (size_t j = 0; j < end - begin; ++j)
              key[j] = static_cast<Index>(data[begin + j]);
            connectivity.polytope(static_cast<Geometry::Polytope::Type>(geometry[i]), std::move(key));
          }
        }

        for (size_t d = 0; d <= D; ++d)
        {
          for (size_t dp = 0; dp <= D; ++dp)
          {
            const auto offsetsPath = incidencePath(d, dp, "Offsets");
            const auto dataPath = incidencePath(d, dp, "Data");
            if (!hasPath(file, offsetsPath) || !hasPath(file, dataPath))
              continue;

            const auto offsets = readVectorDataset<unsigned long long>(file, offsetsPath, H5T_NATIVE_ULLONG);
            const auto data = readVectorDataset<unsigned long long>(file, dataPath, H5T_NATIVE_ULLONG);
            check(offsets.size() == static_cast<size_t>(counts[d]) + 1, "Invalid incidence offsets size.");

            typename Geometry::Connectivity<ContextType>::Incidence incidence;
            incidence.resize(static_cast<size_t>(counts[d]));
            for (size_t i = 0; i < static_cast<size_t>(counts[d]); ++i)
            {
              const auto begin = static_cast<size_t>(offsets[i]);
              const auto end = static_cast<size_t>(offsets[i + 1]);
              check(end >= begin && end <= data.size(), "Invalid incidence offsets.");
              auto& row = incidence[i];
              row.reserve(end - begin);
              for (size_t j = begin; j < end; ++j)
                row.push_back(static_cast<Index>(data[j]));
            }
            connectivity.setIncidence({ d, dp }, std::move(incidence));
          }
        }

        return connectivity;
      }

      template <class MeshType>
      static void saveAttributes(hid_t file, const MeshType& mesh)
      {
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();
        for (size_t d = 0; d <= D; ++d)
        {
          std::vector<unsigned long long> attrs(connectivity.getCount(d), NullAttributeMarker);
          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            if (const auto attr = mesh.getAttribute(d, i))
              attrs[static_cast<size_t>(i)] = static_cast<unsigned long long>(*attr);
          }
          writeVectorDataset(file, attributePath(d), H5T_NATIVE_ULLONG, attrs);
        }
      }

      template <class ContextType>
      static Geometry::AttributeIndex loadAttributes(
          hid_t file,
          const Geometry::Connectivity<ContextType>& connectivity)
      {
        Geometry::AttributeIndex attrs;
        const size_t D = connectivity.getDimension();
        attrs.initialize(D);
        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = connectivity.getCount(d);
          attrs.resize(d, count);
          const auto path = attributePath(d);
          if (!hasPath(file, path))
            continue;

          const auto serialized = readVectorDataset<unsigned long long>(file, path, H5T_NATIVE_ULLONG);
          check(serialized.size() == count, "Invalid attribute dataset size.");
          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            const auto value = serialized[static_cast<size_t>(i)];
            if (value != NullAttributeMarker)
              attrs.set({ d, i }, Optional<Attribute>(static_cast<Attribute>(value)));
          }
        }
        return attrs;
      }
#endif
  };
}

#endif
