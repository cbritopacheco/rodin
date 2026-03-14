/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HDF5_H
#define RODIN_IO_HDF5_H

#include <H5Ipublic.h>
#include <array>
#include <iterator>
#include <string>
#include <vector>
#include <limits>
#include <functional>
#include <type_traits>

#include <boost/filesystem/path.hpp>

#include "ForwardDecls.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/AttributeIndex.h"
#include "Rodin/Geometry/PolytopeTransformationIndex.h"
#include "Rodin/IO/MeshLoader.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/GridFunctionLoader.h"
#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/FormLanguage/Traits.h"

#include <hdf5.h>

namespace Rodin::IO
{
  namespace HDF5
  {
    namespace Path
    {
      static constexpr const char* Mesh = "/Mesh";

      static constexpr const char* Mesh_SpaceDimension = "/Mesh/SpaceDimension";
      static constexpr const char* Mesh_Geometry = "/Mesh/Geometry";
      static constexpr const char* Mesh_Geometry_Vertices = "/Mesh/Geometry/Vertices";

      static constexpr const char* Connectivity = "/Mesh/Connectivity";
      static constexpr const char* ConnectivityD0 = "/Mesh/Connectivity/D_0";
      static constexpr const char* ConnectivityTopologicalDimension = "/Mesh/Connectivity/TopologicalDimension";
      static constexpr const char* ConnectivityCounts = "/Mesh/Connectivity/Counts";
      static constexpr const char* ConnectivityD0Types = "/Mesh/Connectivity/D_0/Types";
      static constexpr const char* ConnectivityD0Offsets = "/Mesh/Connectivity/D_0/Offsets";
      static constexpr const char* ConnectivityD0Indices = "/Mesh/Connectivity/D_0/Indices";

      static constexpr const char* XDMF = "/Mesh/XDMF";
      static constexpr const char* XDMF_Topology = "/Mesh/XDMF/Topology";

      static constexpr const char* Attributes = "/Mesh/Attributes";

      static constexpr const char* GridFunction = "/GridFunction";
      static constexpr const char* GridFunctionMeta = "/GridFunction/Meta";
      static constexpr const char* GridFunctionMetaName = "/GridFunction/Meta/Name";
      static constexpr const char* GridFunctionMetaSize = "/GridFunction/Meta/Size";
      static constexpr const char* GridFunctionMetaDimension = "/GridFunction/Meta/Dimension";
      static constexpr const char* GridFunctionValues = "/GridFunction/Values";
      static constexpr const char* GridFunctionValuesData = "/GridFunction/Values/Data";
    }

    // Sentinel value used in serialized attribute arrays to encode "no attribute".
    static constexpr unsigned long long NullAttributeMarker = std::numeric_limits<unsigned long long>::max();

    static void check(bool condition, const char* msg)
    {
      if (!condition)
        Alert::Exception() << msg << Alert::Raise;
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

    static bool hasPath(hid_t file, const std::string& path)
    {
      return H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
    }

    static std::string attributePath(size_t d)
    {
      return std::string(Path::Attributes) + "/d" + std::to_string(d);
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
  }

  template <>
  class MeshLoader<FileFormat::HDF5, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ObjectType = Geometry::Mesh<Context::Local>;
      using Parent = MeshLoaderBase<Context::Local>;

      MeshLoader(ObjectType& mesh)
        : Parent(mesh)
      {}

      void load(std::istream&) override
      {
        Alert::Exception()
          << "HDF5 mesh loading requires file-path based loading."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        this->load(file);
        H5Fclose(file);
      }

      void load(const hid_t& file)
      {
        auto& mesh = this->getObject();
        HDF5::check(file >= 0, "Failed to open HDF5 file.");

        const auto sdim = static_cast<size_t>(
            HDF5::readScalarDataset<unsigned long long>(file, HDF5::Path::Mesh_SpaceDimension, H5T_NATIVE_ULLONG));

        auto vertices = this->readVertices(file);
        auto connectivity = this->readConnectivity(file);
        auto attributes = this->readAttributes(file, connectivity);

        Geometry::PolytopeTransformationIndex transformations;
        transformations.initialize(connectivity.getDimension());
        for (size_t d = 0; d <= connectivity.getDimension(); ++d)
          transformations.resize(d, connectivity.getCount(d));

        mesh = mesh.Build()
          .initialize(sdim)
          .setVertices(std::move(vertices))
          .setConnectivity(std::move(connectivity))
          .setAttributeIndex(std::move(attributes))
          .setTransformationIndex(std::move(transformations))
          .finalize();
      }

      Geometry::PointCloud readVertices(hid_t file)
      {
        const auto dset = H5Dopen2(file, HDF5::Path::Mesh_Geometry_Vertices, H5P_DEFAULT);
        HDF5::check(dset >= 0, "Failed to open /Mesh/Vertices dataset.");

        const auto space = H5Dget_space(dset);
        HDF5::check(space >= 0, "Failed to open /Mesh/Vertices dataspace.");

        hsize_t dims[2] = { 0, 0 };
        const int rank = H5Sget_simple_extent_ndims(space);
        HDF5::check(rank == 2, "Invalid /Mesh/Vertices dataset rank.");
        const auto readDims = H5Sget_simple_extent_dims(space, dims, nullptr);
        HDF5::check(readDims == 2, "Invalid /Mesh/Vertices dataset dimensions.");

        std::vector<double> packed(static_cast<size_t>(dims[0] * dims[1]));
        if (!packed.empty())
        {
          const auto status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, packed.data());
          HDF5::check(status >= 0, "Failed to read /Mesh/Vertices dataset.");
        }

        Geometry::PointCloud vertices(static_cast<std::uint8_t>(dims[1]), static_cast<size_t>(dims[0]));
        for (size_t i = 0; i < static_cast<size_t>(dims[0]); ++i)
          for (size_t d = 0; d < static_cast<size_t>(dims[1]); ++d)
            vertices(static_cast<std::uint8_t>(d), i) = packed[i * static_cast<size_t>(dims[1]) + d];

        H5Sclose(space);
        H5Dclose(dset);
        return vertices;
      }

      Geometry::Connectivity<Context::Local> readConnectivity(hid_t file)
      {
        Geometry::Connectivity<Context::Local> connectivity;

        const size_t D = static_cast<size_t>(
            HDF5::readScalarDataset<unsigned long long>(
                file,
                HDF5::Path::ConnectivityTopologicalDimension,
                H5T_NATIVE_ULLONG));
        connectivity.initialize(D);

        auto counts = HDF5::readVectorDataset<unsigned long long>(file, HDF5::Path::ConnectivityCounts, H5T_NATIVE_ULLONG);
        HDF5::check(counts.size() == D + 1, "Invalid /Mesh/Connectivity/Counts size.");
        connectivity.nodes(static_cast<size_t>(counts[0]));
        const auto types = HDF5::readVectorDataset<int>(file, HDF5::Path::ConnectivityD0Types, H5T_NATIVE_INT);
        const auto offsets = HDF5::readVectorDataset<unsigned long long>(file, HDF5::Path::ConnectivityD0Offsets, H5T_NATIVE_ULLONG);
        const auto indices = HDF5::readVectorDataset<unsigned long long>(file, HDF5::Path::ConnectivityD0Indices, H5T_NATIVE_ULLONG);
        HDF5::check(offsets.size() == types.size() + 1, "Invalid /Mesh/Connectivity/D_0 offsets size.");
        counts[D] = static_cast<unsigned long long>(types.size());
        for (size_t d = 1; d <= D; ++d)
          connectivity.reserve(d, static_cast<size_t>(counts[d]));
        for (size_t i = 0; i < types.size(); ++i)
        {
          const auto begin = static_cast<size_t>(offsets[i]);
          const auto end = static_cast<size_t>(offsets[i + 1]);
          HDF5::check(end >= begin && end <= indices.size(), "Invalid /Mesh/Connectivity/D_0 offsets.");
          Geometry::Polytope::Key key(end - begin);
          for (size_t j = 0; j < end - begin; ++j)
            key[j] = static_cast<Index>(indices[begin + j]);
          connectivity.polytope(static_cast<Geometry::Polytope::Type>(types[i]), std::move(key));
        }
        return connectivity;
      }

      Geometry::AttributeIndex readAttributes(hid_t file, const Geometry::Connectivity<ContextType>& connectivity)
      {
        Geometry::AttributeIndex attrs;
        const size_t D = connectivity.getDimension();
        attrs.initialize(D);
        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = connectivity.getCount(d);
          attrs.resize(d, count);
          const auto path = HDF5::attributePath(d);
          if (!HDF5::hasPath(file, path))
            continue;

          const auto serialized = HDF5::readVectorDataset<unsigned long long>(file, path, H5T_NATIVE_ULLONG);
          HDF5::check(serialized.size() == count, "Invalid attribute dataset size.");
          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            const auto value = serialized[static_cast<size_t>(i)];
            if (value != HDF5::NullAttributeMarker)
              attrs.set({ d, i }, Optional<Geometry::Attribute>(static_cast<Geometry::Attribute>(value)));
          }
        }
        return attrs;
      }
  };

  template <>
  class MeshPrinter<FileFormat::HDF5, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ObjectType = Geometry::Mesh<Context::Local>;

      using Parent = MeshPrinterBase<Context::Local>;

      MeshPrinter(const ObjectType& mesh)
        : Parent(mesh)
      {}

      void print(std::ostream&) override
      {
        Alert::Exception()
          << "HDF5 mesh printing requires file-path based printing."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename) override
      {
        const auto& mesh = this->getObject();

        const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(file >= 0, "Failed to create HDF5 file.");

        const auto meshGroup = H5Gcreate2(file, HDF5::Path::Mesh, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(meshGroup >= 0, "Failed to create /Mesh group.");
        H5Gclose(meshGroup);

        const auto geometryGroup = H5Gcreate2(file, HDF5::Path::Mesh_Geometry, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(geometryGroup >= 0, "Failed to create /Mesh/Geometry group.");
        H5Gclose(geometryGroup);

        const auto connectivityGroup = H5Gcreate2(file, HDF5::Path::Connectivity, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(connectivityGroup >= 0, "Failed to create /Mesh/Connectivity group.");
        H5Gclose(connectivityGroup);

        const auto d0Group = H5Gcreate2(file, HDF5::Path::ConnectivityD0, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(d0Group >= 0, "Failed to create /Mesh/Connectivity/D_0 group.");
        H5Gclose(d0Group);

        // const auto xdmfGroup = H5Gcreate2(file, HDF5::Path::XDMF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // HDF5::check(xdmfGroup >= 0, "Failed to create /Mesh/XDMF group.");
        // H5Gclose(xdmfGroup);

        const auto attributesGroup = H5Gcreate2(file, HDF5::Path::Attributes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(attributesGroup >= 0, "Failed to create /Mesh/Attributes group.");
        H5Gclose(attributesGroup);

        HDF5::writeScalarDataset(
            file,
            HDF5::Path::Mesh_SpaceDimension,
            H5T_NATIVE_ULLONG,
            static_cast<unsigned long long>(mesh.getSpaceDimension()));

        this->writeVertices(file);
        this->writeConnectivity(file);
        this->writeAttributes(file);

        H5Fclose(file);
      }

      void writeVertices(hid_t file)
      {
        const auto& mesh = this->getObject();
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
        HDF5::writeDataset(file, HDF5::Path::Mesh_Geometry_Vertices, H5T_NATIVE_DOUBLE, dims, 2, vertices.data());
      }

      void writeConnectivity(hid_t file)
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();
        HDF5::writeScalarDataset(
            file,
            HDF5::Path::ConnectivityTopologicalDimension,
            H5T_NATIVE_ULLONG,
            static_cast<unsigned long long>(D));

        std::vector<unsigned long long> counts(D + 1, 0);
        counts[0] = static_cast<unsigned long long>(connectivity.getCount(0));
        counts[D] = static_cast<unsigned long long>(connectivity.getCount(D));
        HDF5::writeVectorDataset(file, HDF5::Path::ConnectivityCounts, H5T_NATIVE_ULLONG, counts);

        const auto& d0 = connectivity.getIncidence(D, 0);
        HDF5::check(d0.size() == connectivity.getCount(D), "Missing canonical D_0 connectivity.");
        std::vector<int> types;
        std::vector<unsigned long long> offsets(1, 0);
        std::vector<unsigned long long> indices;
        types.reserve(connectivity.getCount(D));
        for (Index i = 0; i < static_cast<Index>(connectivity.getCount(D)); ++i)
        {
          types.push_back(static_cast<int>(connectivity.getGeometry(D, i)));
          for (const auto v : d0[static_cast<size_t>(i)])
            indices.push_back(static_cast<unsigned long long>(v));
          offsets.push_back(static_cast<unsigned long long>(indices.size()));
        }

        HDF5::writeVectorDataset(file, HDF5::Path::ConnectivityD0Types, H5T_NATIVE_INT, types);
        HDF5::writeVectorDataset(file, HDF5::Path::ConnectivityD0Offsets, H5T_NATIVE_ULLONG, offsets);
        HDF5::writeVectorDataset(file, HDF5::Path::ConnectivityD0Indices, H5T_NATIVE_ULLONG, indices);
      }

      void writeAttributes(hid_t file)
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();
        for (size_t d = 0; d <= D; ++d)
        {
          if (d != 0 && d != D)
            continue;
          std::vector<unsigned long long> attrs(connectivity.getCount(d), HDF5::NullAttributeMarker);
          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            if (const auto attr = mesh.getAttribute(d, i))
              attrs[static_cast<size_t>(i)] = static_cast<unsigned long long>(*attr);
          }
          HDF5::writeVectorDataset(file, HDF5::attributePath(d), H5T_NATIVE_ULLONG, attrs);
        }
      }
  };

  template <class FES, class Scalar>
  class GridFunctionLoader<FileFormat::HDF5, FES, Math::Vector<Scalar>>
    : public GridFunctionLoaderBase<FES, Math::Vector<Scalar>>
  {
    public:
      using DataType = Math::Vector<Scalar>;

      using ObjectType = Variational::GridFunction<FES, DataType>;

      using Parent = GridFunctionLoaderBase<FES, DataType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction loading is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("load(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        auto& gf = this->getObject();
        const auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        HDF5::check(file >= 0, "Failed to open HDF5 file.");

        const auto values = HDF5::readVectorDataset<double>(file, HDF5::Path::GridFunctionValuesData, H5T_NATIVE_DOUBLE);
        const auto dofCount = static_cast<size_t>(
            HDF5::readScalarDataset<unsigned long long>(file, HDF5::Path::GridFunctionMetaSize, H5T_NATIVE_ULLONG));
        const auto vectorDim = static_cast<size_t>(
            HDF5::readScalarDataset<unsigned long long>(file, HDF5::Path::GridFunctionMetaDimension, H5T_NATIVE_ULLONG));
        check(values.size() == gf.getData().size(), "Invalid GridFunction data size.");
        check(dofCount == gf.getSize(), "GridFunction size mismatch.");
        check(vectorDim == gf.getDimension(), "GridFunction dimension mismatch.");

        auto& data = gf.getData();
        for (size_t i = 0; i < values.size(); ++i)
          data[i] = static_cast<typename std::remove_reference_t<decltype(data[0])>>(values[i]);

        H5Fclose(file);
      }
  };

  template <class FES, class Scalar>
  class GridFunctionPrinter<FileFormat::HDF5, FES, Math::Vector<Scalar>> final
    : public GridFunctionPrinterBase<FileFormat::HDF5, FES, Math::Vector<Scalar>>
  {
    public:
      using DataType = Math::Vector<Scalar>;

      using ObjectType = Variational::GridFunction<FES, DataType>;

      using Parent = GridFunctionPrinterBase<FileFormat::HDF5, FES, DataType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream&)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction printing is file-path based."
          << "Please use the "
          << Alert::Identifier::Function("print(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename)
      {
        const auto& gf = this->getObject();

        const auto file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(file >= 0, "Failed to create HDF5 file.");

        const auto group = H5Gcreate2(file, HDF5::Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(group >= 0, "Failed to create /GridFunction group.");
        H5Gclose(group);

        const auto metaGroup = H5Gcreate2(file, HDF5::Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(metaGroup >= 0, "Failed to create /GridFunction/Meta group.");
        H5Gclose(metaGroup);

        const auto valuesGroup = H5Gcreate2(file, HDF5::Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5::check(valuesGroup >= 0, "Failed to create /GridFunction/Values group.");
        H5Gclose(valuesGroup);

        const auto& data = gf.getData();
        std::vector<double> values(static_cast<size_t>(data.size()));
        for (size_t i = 0; i < values.size(); ++i)
          values[i] = static_cast<double>(data[i]);

        const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
        HDF5::writeDataset(file, HDF5::Path::GridFunctionValuesData, H5T_NATIVE_DOUBLE, dims, 1, values.data());

        const unsigned long long dofCount = static_cast<unsigned long long>(gf.getSize());
        const unsigned long long vectorDim = static_cast<unsigned long long>(gf.getDimension());
        HDF5::writeScalarDataset(file, HDF5::Path::GridFunctionMetaSize, H5T_NATIVE_ULLONG, dofCount);
        HDF5::writeScalarDataset(file, HDF5::Path::GridFunctionMetaDimension, H5T_NATIVE_ULLONG, vectorDim);

        H5Fclose(file);
      }
  };
}

#endif
