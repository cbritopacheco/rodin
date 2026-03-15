/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_HDF5_H
#define RODIN_IO_HDF5_H

#include <hdf5.h>

#include <array>
#include <limits>
#include <string>
#include <vector>
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

namespace Rodin::IO
{
  namespace HDF5
  {
    using U64 = unsigned long long;
    using I32 = int;
    using F64 = double;
    using U8  = unsigned char;

    static constexpr U64 NullAttributeMarker = std::numeric_limits<U64>::max();

    namespace Path
    {
      static constexpr const char* Mesh = "/Mesh";

      static constexpr const char* MeshMeta = "/Mesh/Meta";
      static constexpr const char* MeshMetaSpaceDimension = "/Mesh/Meta/SpaceDimension";

      static constexpr const char* MeshGeometry = "/Mesh/Geometry";
      static constexpr const char* MeshGeometryVertices = "/Mesh/Geometry/Vertices";

      static constexpr const char* MeshConnectivity = "/Mesh/Connectivity";
      static constexpr const char* MeshConnectivityMeta = "/Mesh/Connectivity/Meta";
      static constexpr const char* MeshConnectivityMetaMaximalDimension = "/Mesh/Connectivity/Meta/MaximalDimension";
      static constexpr const char* MeshConnectivityMetaDimension = "/Mesh/Connectivity/Meta/Dimension";

      static constexpr const char* MeshConnectivityCounts = "/Mesh/Connectivity/Counts";
      static constexpr const char* MeshConnectivityCountsByDimension = "/Mesh/Connectivity/Counts/ByDimension";
      static constexpr const char* MeshConnectivityCountsByGeometry = "/Mesh/Connectivity/Counts/ByGeometry";

      static constexpr const char* MeshConnectivityEntities = "/Mesh/Connectivity/Entities";
      static constexpr const char* MeshConnectivityState = "/Mesh/Connectivity/State";
      static constexpr const char* MeshConnectivityStatePresent = "/Mesh/Connectivity/State/Present";
      static constexpr const char* MeshConnectivityStateDirty = "/Mesh/Connectivity/State/Dirty";
      static constexpr const char* MeshConnectivityIncidence = "/Mesh/Connectivity/Incidence";

      static constexpr const char* MeshAttributes = "/Mesh/Attributes";
      static constexpr const char* MeshTransformations = "/Mesh/Transformations";

      static constexpr const char* MeshXDMF = "/Mesh/XDMF";
      static constexpr const char* MeshXDMFTopology = "/Mesh/XDMF/Topology";
      static constexpr const char* MeshXDMFTopologySize = "/Mesh/XDMF/TopologySize";

      static constexpr const char* GridFunction = "/GridFunction";
      static constexpr const char* GridFunctionMeta = "/GridFunction/Meta";
      static constexpr const char* GridFunctionMetaName = "/GridFunction/Meta/Name";
      static constexpr const char* GridFunctionMetaSize = "/GridFunction/Meta/Size";
      static constexpr const char* GridFunctionMetaDimension = "/GridFunction/Meta/Dimension";
      static constexpr const char* GridFunctionValues = "/GridFunction/Values";
      static constexpr const char* GridFunctionValuesData = "/GridFunction/Values/Data";
    }

    class Handle
    {
      public:
        Handle()
          : m_id(-1),
            m_close(nullptr)
        {}

        Handle(hid_t id, herr_t (*closeFn)(hid_t))
          : m_id(id),
            m_close(closeFn)
        {}

        Handle(const Handle&) = delete;
        Handle& operator=(const Handle&) = delete;

        Handle(Handle&& other) noexcept
          : m_id(other.m_id),
            m_close(other.m_close)
        {
          other.m_id = -1;
          other.m_close = nullptr;
        }

        Handle& operator=(Handle&& other) noexcept
        {
          if (this != &other)
          {
            reset();
            m_id = other.m_id;
            m_close = other.m_close;
            other.m_id = -1;
            other.m_close = nullptr;
          }
          return *this;
        }

        ~Handle()
        {
          reset();
        }

        void reset()
        {
          if (m_id >= 0 && m_close)
            m_close(m_id);
          m_id = -1;
          m_close = nullptr;
        }

        hid_t get() const
        {
          return m_id;
        }

        explicit operator bool() const
        {
          return m_id >= 0;
        }

      private:
        hid_t m_id;
        herr_t (*m_close)(hid_t);
    };

    inline
    Handle File(hid_t id)
    {
      return Handle(id, H5Fclose);
    }

    inline
    Handle Group(hid_t id)
    {
      return Handle(id, H5Gclose);
    }

    inline
    Handle DataSet(hid_t id)
    {
      return Handle(id, H5Dclose);
    }

    inline
    Handle Space(hid_t id)
    {
      return Handle(id, H5Sclose);
    }

    template <class T>
    hid_t getNativeType();

    template <>
    inline
    hid_t getNativeType<U64>()
    {
      return H5T_NATIVE_ULLONG;
    }

    template <>
    inline
    hid_t getNativeType<I32>()
    {
      return H5T_NATIVE_INT;
    }

    template <>
    inline
    hid_t getNativeType<F64>()
    {
      return H5T_NATIVE_DOUBLE;
    }

    template <>
    inline
    hid_t getNativeType<U8>()
    {
      return H5T_NATIVE_UCHAR;
    }

    inline
    bool exists(hid_t loc, const std::string& path)
    {
      return H5Lexists(loc, path.c_str(), H5P_DEFAULT) > 0;
    }

    inline
    std::string entityGroupPath(size_t d)
    {
      return std::string(Path::MeshConnectivityEntities) + "/" + std::to_string(d);
    }

    inline
    std::string entityTypesPath(size_t d)
    {
      return entityGroupPath(d) + "/Types";
    }

    inline
    std::string entityOffsetsPath(size_t d)
    {
      return entityGroupPath(d) + "/Offsets";
    }

    inline
    std::string entityIndicesPath(size_t d)
    {
      return entityGroupPath(d) + "/Indices";
    }

    inline
    std::string incidenceGroupPath(size_t d, size_t dp)
    {
      return std::string(Path::MeshConnectivityIncidence) + "/" + std::to_string(d) + "_" + std::to_string(dp);
    }

    inline
    std::string incidenceOffsetsPath(size_t d, size_t dp)
    {
      return incidenceGroupPath(d, dp) + "/Offsets";
    }

    inline
    std::string incidenceIndicesPath(size_t d, size_t dp)
    {
      return incidenceGroupPath(d, dp) + "/Indices";
    }

    inline
    std::string attributePath(size_t d)
    {
      return std::string(Path::MeshAttributes) + "/" + std::to_string(d);
    }

    inline
    std::string transformationGroupPath(size_t d)
    {
      return std::string(Path::MeshTransformations) + "/" + std::to_string(d);
    }

    inline
    std::string transformationKindPath(size_t d)
    {
      return transformationGroupPath(d) + "/Kind";
    }

    inline
    size_t getGeometryCountArraySize()
    {
      using PT = Geometry::Polytope::Type;
      return static_cast<size_t>(PT::Hexahedron) + 1;
    }

    inline
    U64 getXDMFMixedTopologyId(Geometry::Polytope::Type t)
    {
      using PT = Geometry::Polytope::Type;
      switch (t)
      {
        case PT::Point:         return 1;
        case PT::Segment:       return 2;
        case PT::Triangle:      return 4;
        case PT::Quadrilateral: return 5;
        case PT::Tetrahedron:   return 6;
        case PT::Wedge:         return 8;
        case PT::Hexahedron:    return 9;
      }

      Alert::Exception()
        << "Unsupported polytope type for XDMF mixed topology."
        << Alert::Raise;
    }

    template <class T>
    std::vector<T> readVectorDataset(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      const auto count = static_cast<size_t>(H5Sget_simple_extent_npoints(space.get()));
      std::vector<T> values(count);
      if (count > 0)
      {
        const auto status = H5Dread(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to read HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
      return values;
    }

    template <class T>
    T readScalarDataset(hid_t file, const std::string& path)
    {
      const auto values = readVectorDataset<T>(file, path);
      if (values.size() != 1)
      {
        Alert::Exception()
          << "Expected scalar HDF5 dataset at path: " << path
          << Alert::Raise;
      }
      return values[0];
    }

    template <class T>
    void writeVectorDataset(hid_t file, const std::string& path, const std::vector<T>& values)
    {
      const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
      const auto space = Space(H5Screate_simple(1, dims, nullptr));
      if (!space)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataspace for dataset: " << path
          << Alert::Raise;
      }

      const auto dset = DataSet(H5Dcreate2(
          file,
          path.c_str(),
          getNativeType<T>(),
          space.get(),
          H5P_DEFAULT,
          H5P_DEFAULT,
          H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataset: " << path
          << Alert::Raise;
      }

      if (!values.empty())
      {
        const auto status = H5Dwrite(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to write HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
    }

    template <class T>
    void writeScalarDataset(hid_t file, const std::string& path, const T& value)
    {
      writeVectorDataset(file, path, std::vector<T>{ value });
    }

    template <class T>
    void writeMatrixDataset(
        hid_t file,
        const std::string& path,
        const std::vector<T>& values,
        hsize_t rows,
        hsize_t cols)
    {
      if (values.size() != static_cast<size_t>(rows * cols))
      {
        Alert::Exception()
          << "Invalid HDF5 matrix payload size for dataset: " << path
          << Alert::Raise;
      }

      const hsize_t dims[2] = { rows, cols };
      const auto space = Space(H5Screate_simple(2, dims, nullptr));
      if (!space)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataspace for dataset: " << path
          << Alert::Raise;
      }

      const auto dset = DataSet(H5Dcreate2(
          file,
          path.c_str(),
          getNativeType<T>(),
          space.get(),
          H5P_DEFAULT,
          H5P_DEFAULT,
          H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataset: " << path
          << Alert::Raise;
      }

      if (!values.empty())
      {
        const auto status = H5Dwrite(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to write HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
    }

    inline
    std::pair<hsize_t, hsize_t> readMatrixShape(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      hsize_t dims[2] = { 0, 0 };
      const int rank = H5Sget_simple_extent_ndims(space.get());
      if (rank != 2)
      {
        Alert::Exception()
          << "Expected rank-2 HDF5 dataset at path: " << path
          << Alert::Raise;
      }

      const auto status = H5Sget_simple_extent_dims(space.get(), dims, nullptr);
      if (status != 2)
      {
        Alert::Exception()
          << "Failed to read matrix shape for HDF5 dataset: " << path
          << Alert::Raise;
      }

      return { dims[0], dims[1] };
    }

    inline
    std::vector<hsize_t> readDatasetShape(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      const int rank = H5Sget_simple_extent_ndims(space.get());
      if (rank < 0)
      {
        Alert::Exception()
          << "Failed to read rank of HDF5 dataset: " << path
          << Alert::Raise;
      }

      std::vector<hsize_t> dims(static_cast<size_t>(rank), 0);
      if (rank > 0)
      {
        const auto status = H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr);
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to read shape of HDF5 dataset: " << path
            << Alert::Raise;
        }
      }

      return dims;
    }
  }

  template <>
  class MeshLoader<FileFormat::HDF5, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshLoaderBase<ContextType>;

      explicit
      MeshLoader(ObjectType& mesh)
        : Parent(mesh)
      {}

      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 mesh loading requires file-path based loading."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        const auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open HDF5 mesh file: " << filename
            << Alert::Raise;
        }

        auto& mesh = this->getObject();

        const size_t sdim = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::MeshMetaSpaceDimension));

        const auto vertices = this->readVertices(file.get());
        auto connectivity = this->readConnectivity(file.get());
        auto attributes = this->readAttributes(file.get(), connectivity);

        Geometry::PolytopeTransformationIndex transformations;
        transformations.initialize(connectivity.getMaximalDimension());
        for (size_t d = 0; d <= connectivity.getMaximalDimension(); ++d)
          transformations.resize(d, connectivity.getCount(d));

        mesh = ObjectType::Build()
          .initialize(sdim)
          .setVertices(vertices)
          .setConnectivity(std::move(connectivity))
          .setAttributeIndex(std::move(attributes))
          .setTransformationIndex(std::move(transformations))
          .finalize();
      }

    private:
      Geometry::PointCloud readVertices(hid_t file) const
      {
        const auto [nv, sdim] = HDF5::readMatrixShape(file, HDF5::Path::MeshGeometryVertices);
        const auto packed = HDF5::readVectorDataset<HDF5::F64>(file, HDF5::Path::MeshGeometryVertices);
        if (packed.size() != static_cast<size_t>(nv * sdim))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid vertex payload size in HDF5 mesh file."
            << Alert::Raise;
        }

        Geometry::PointCloud vertices(static_cast<std::uint8_t>(sdim), static_cast<size_t>(nv));
        for (size_t i = 0; i < static_cast<size_t>(nv); ++i)
        {
          for (size_t d = 0; d < static_cast<size_t>(sdim); ++d)
            vertices(static_cast<std::uint8_t>(d), i) = packed[i * static_cast<size_t>(sdim) + d];
        }
        return vertices;
      }

      Geometry::Connectivity<ContextType> readConnectivity(hid_t file) const
      {
        Geometry::Connectivity<ContextType> connectivity;

        const size_t Dmax = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file, HDF5::Path::MeshConnectivityMetaMaximalDimension));
        const size_t D = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file, HDF5::Path::MeshConnectivityMetaDimension));

        const auto byDimension = HDF5::readVectorDataset<HDF5::U64>(
            file,
            HDF5::Path::MeshConnectivityCountsByDimension);
        if (byDimension.size() != Dmax + 1)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/Counts/ByDimension dataset size."
            << Alert::Raise;
        }

        connectivity.initialize(Dmax);
        connectivity.nodes(static_cast<size_t>(byDimension[0]));

        for (size_t d = 1; d <= Dmax; ++d)
        {
          const auto types = HDF5::readVectorDataset<HDF5::I32>(file, HDF5::entityTypesPath(d));
          const auto offsets = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::entityOffsetsPath(d));
          const auto indices = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::entityIndicesPath(d));

          if (types.size() != static_cast<size_t>(byDimension[d]))
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid entity type count for dimension " << d << "."
              << Alert::Raise;
          }

          if (offsets.size() != types.size() + 1)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid CSR offsets for entity dimension " << d << "."
              << Alert::Raise;
          }

          connectivity.reserve(d, static_cast<size_t>(byDimension[d]));
          for (size_t i = 0; i < types.size(); ++i)
          {
            const size_t begin = static_cast<size_t>(offsets[i]);
            const size_t end = static_cast<size_t>(offsets[i + 1]);
            if (begin > end || end > indices.size())
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid CSR entity offsets for dimension " << d << "."
                << Alert::Raise;
            }

            Geometry::Polytope::Key key(end - begin);
            for (size_t k = begin; k < end; ++k)
              key[k - begin] = static_cast<Index>(indices[k]);

            connectivity.polytope(
                static_cast<Geometry::Polytope::Type>(types[i]),
                std::move(key));
          }
        }

        const auto present = HDF5::readVectorDataset<HDF5::U8>(
            file,
            HDF5::Path::MeshConnectivityStatePresent);
        const auto dirty = HDF5::readVectorDataset<HDF5::U8>(
            file,
            HDF5::Path::MeshConnectivityStateDirty);

        if (present.size() != (Dmax + 1) * (Dmax + 1))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/State/Present dataset size."
            << Alert::Raise;
        }

        if (dirty.size() != (Dmax + 1) * (Dmax + 1))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/State/Dirty dataset size."
            << Alert::Raise;
        }

        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            connectivity.setDirty(d, dp, dirty[flat] != 0);
            if (!present[flat])
              continue;

            const auto offsets = HDF5::readVectorDataset<HDF5::U64>(
                file,
                HDF5::incidenceOffsetsPath(d, dp));
            const auto indices = HDF5::readVectorDataset<HDF5::U64>(
                file,
                HDF5::incidenceIndicesPath(d, dp));

            if (offsets.size() != connectivity.getCount(d) + 1)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid CSR offsets for incidence "
                << d << " -> " << dp << "."
                << Alert::Raise;
            }

            Geometry::Incidence inc;
            inc.resize(connectivity.getCount(d));
            for (size_t i = 0; i < connectivity.getCount(d); ++i)
            {
              const size_t begin = static_cast<size_t>(offsets[i]);
              const size_t end = static_cast<size_t>(offsets[i + 1]);
              if (begin > end || end > indices.size())
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Invalid CSR row bounds for incidence "
                  << d << " -> " << dp << "."
                  << Alert::Raise;
              }

              auto& row = inc[i];
              row.reserve(end - begin);
              for (size_t k = begin; k < end; ++k)
                row.push_back(static_cast<Index>(indices[k]));
            }

            connectivity.setIncidence({ d, dp }, std::move(inc));
          }
        }

        return connectivity;
      }

      Geometry::AttributeIndex readAttributes(
          hid_t file,
          const Geometry::Connectivity<ContextType>& connectivity) const
      {
        Geometry::AttributeIndex attrs;
        const size_t D = connectivity.getDimension();
        attrs.initialize(D);

        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = connectivity.getCount(d);
          attrs.resize(d, count);

          const auto path = HDF5::attributePath(d);
          if (!HDF5::exists(file, path))
            continue;

          const auto values = HDF5::readVectorDataset<HDF5::U64>(file, path);
          if (values.size() != count)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid attribute dataset size for dimension " << d << "."
              << Alert::Raise;
          }

          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (values[static_cast<size_t>(i)] != HDF5::NullAttributeMarker)
            {
              attrs.set(
                  { d, i },
                  Optional<Geometry::Attribute>(
                    static_cast<Geometry::Attribute>(values[static_cast<size_t>(i)])));
            }
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
      using ContextType = Context::Local;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshPrinterBase<ContextType>;

      explicit
      MeshPrinter(const ObjectType& mesh)
        : Parent(mesh)
      {}

      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 mesh printing requires file-path based printing."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename) override
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();

        const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to create HDF5 mesh file: " << filename
            << Alert::Raise;
        }

        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::Mesh, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Meta group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshGeometry, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Geometry group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivity, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivityMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity/Meta group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivityCounts, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity/Counts group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivityEntities, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity/Entities group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivityState, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity/State group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshConnectivityIncidence, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Connectivity/Incidence group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshAttributes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Attributes group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshTransformations, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/Transformations group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::MeshXDMF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /Mesh/XDMF group."
              << Alert::Raise;
          }
        }

        HDF5::writeScalarDataset(
            file.get(),
            HDF5::Path::MeshMetaSpaceDimension,
            static_cast<HDF5::U64>(mesh.getSpaceDimension()));

        this->writeVertices(file.get());
        this->writeConnectivity(file.get());
        this->writeAttributes(file.get());
        this->writeTransformations(file.get());
        this->writeXDMFTopology(file.get());
      }

    private:
      void writeVertices(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const size_t nv = mesh.getVertexCount();
        const size_t sdim = mesh.getSpaceDimension();

        std::vector<HDF5::F64> packed;
        packed.resize(nv * sdim);
        for (Index i = 0; i < static_cast<Index>(nv); ++i)
        {
          const auto x = mesh.getVertexCoordinates(i);
          for (size_t d = 0; d < sdim; ++d)
            packed[static_cast<size_t>(i) * sdim + d] = static_cast<HDF5::F64>(x(d));
        }

        HDF5::writeMatrixDataset(
            file,
            HDF5::Path::MeshGeometryVertices,
            packed,
            static_cast<hsize_t>(nv),
            static_cast<hsize_t>(sdim));
      }

      void writeConnectivity(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t Dmax = connectivity.getMaximalDimension();
        const size_t D = connectivity.getDimension();

        HDF5::writeScalarDataset(
            file,
            HDF5::Path::MeshConnectivityMetaMaximalDimension,
            static_cast<HDF5::U64>(Dmax));
        HDF5::writeScalarDataset(
            file,
            HDF5::Path::MeshConnectivityMetaDimension,
            static_cast<HDF5::U64>(D));

        std::vector<HDF5::U64> byDimension(Dmax + 1, 0);
        for (size_t d = 0; d <= Dmax; ++d)
          byDimension[d] = static_cast<HDF5::U64>(connectivity.getCount(d));
        HDF5::writeVectorDataset(file, HDF5::Path::MeshConnectivityCountsByDimension, byDimension);

        std::vector<HDF5::U64> byGeometry(HDF5::getGeometryCountArraySize(), 0);
        for (size_t i = 0; i < byGeometry.size(); ++i)
          byGeometry[i] = 0;
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Point)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Point));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Segment)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Segment));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Triangle)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Triangle));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Quadrilateral)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Quadrilateral));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Tetrahedron)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Tetrahedron));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Wedge)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Wedge));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Hexahedron)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Hexahedron));
        HDF5::writeVectorDataset(file, HDF5::Path::MeshConnectivityCountsByGeometry, byGeometry);

        for (size_t d = 1; d <= Dmax; ++d)
        {
          const auto group = HDF5::Group(H5Gcreate2(
              file,
              HDF5::entityGroupPath(d).c_str(),
              H5P_DEFAULT,
              H5P_DEFAULT,
              H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create entity group for dimension " << d << "."
              << Alert::Raise;
          }

          std::vector<HDF5::I32> types;
          std::vector<HDF5::U64> offsets;
          std::vector<HDF5::U64> indices;

          types.reserve(connectivity.getCount(d));
          offsets.reserve(connectivity.getCount(d) + 1);
          offsets.push_back(0);

          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            types.push_back(static_cast<HDF5::I32>(connectivity.getGeometry(d, i)));
            const auto& key = connectivity.getPolytope(d, i);
            for (size_t k = 0; k < key.size(); ++k)
              indices.push_back(static_cast<HDF5::U64>(key[k]));
            offsets.push_back(static_cast<HDF5::U64>(indices.size()));
          }

          HDF5::writeVectorDataset(file, HDF5::entityTypesPath(d), types);
          HDF5::writeVectorDataset(file, HDF5::entityOffsetsPath(d), offsets);
          HDF5::writeVectorDataset(file, HDF5::entityIndicesPath(d), indices);
        }

        std::vector<HDF5::U8> present((Dmax + 1) * (Dmax + 1), 0);
        std::vector<HDF5::U8> dirty((Dmax + 1) * (Dmax + 1), 0);
        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            const auto& inc = connectivity.getIncidence(d, dp);
            present[flat] = static_cast<HDF5::U8>(inc.size() == connectivity.getCount(d) ? 1 : 0);
            dirty[flat] = static_cast<HDF5::U8>(connectivity.isDirty(d, dp) ? 1 : 0);
          }
        }

        HDF5::writeMatrixDataset(
            file,
            HDF5::Path::MeshConnectivityStatePresent,
            present,
            static_cast<hsize_t>(Dmax + 1),
            static_cast<hsize_t>(Dmax + 1));
        HDF5::writeMatrixDataset(
            file,
            HDF5::Path::MeshConnectivityStateDirty,
            dirty,
            static_cast<hsize_t>(Dmax + 1),
            static_cast<hsize_t>(Dmax + 1));

        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            if (!present[flat])
              continue;

            const auto group = HDF5::Group(H5Gcreate2(
                file,
                HDF5::incidenceGroupPath(d, dp).c_str(),
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT));
            if (!group)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to create incidence group for "
                << d << " -> " << dp << "."
                << Alert::Raise;
            }

            const auto& inc = connectivity.getIncidence(d, dp);
            std::vector<HDF5::U64> offsets;
            std::vector<HDF5::U64> indices;
            offsets.reserve(connectivity.getCount(d) + 1);
            offsets.push_back(0);

            for (size_t i = 0; i < connectivity.getCount(d); ++i)
            {
              const auto& row = inc[i];
              for (const auto j : row)
                indices.push_back(static_cast<HDF5::U64>(j));
              offsets.push_back(static_cast<HDF5::U64>(indices.size()));
            }

            HDF5::writeVectorDataset(file, HDF5::incidenceOffsetsPath(d, dp), offsets);
            HDF5::writeVectorDataset(file, HDF5::incidenceIndicesPath(d, dp), indices);
          }
        }
      }

      void writeAttributes(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();

        for (size_t d = 0; d <= D; ++d)
        {
          std::vector<HDF5::U64> attrs(connectivity.getCount(d), HDF5::NullAttributeMarker);
          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            if (const auto attr = mesh.getAttribute(d, i))
              attrs[static_cast<size_t>(i)] = static_cast<HDF5::U64>(*attr);
          }
          HDF5::writeVectorDataset(file, HDF5::attributePath(d), attrs);
        }
      }

      void writeTransformations(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();
        const size_t Dmax = connectivity.getMaximalDimension();

        for (size_t d = 0; d <= Dmax; ++d)
        {
          const auto group = HDF5::Group(H5Gcreate2(
              file,
              HDF5::transformationGroupPath(d).c_str(),
              H5P_DEFAULT,
              H5P_DEFAULT,
              H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create transformation group for dimension " << d << "."
              << Alert::Raise;
          }

          std::vector<HDF5::I32> kind(connectivity.getCount(d), 0);
          HDF5::writeVectorDataset(file, HDF5::transformationKindPath(d), kind);
        }
      }

      void writeXDMFTopology(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();

        std::vector<HDF5::U64> topology;
        for (Index i = 0; i < static_cast<Index>(connectivity.getCount(D)); ++i)
        {
          const auto geometry = connectivity.getGeometry(D, i);
          const auto& key = connectivity.getPolytope(D, i);

          topology.push_back(HDF5::getXDMFMixedTopologyId(geometry));

          if (geometry == Geometry::Polytope::Type::Segment)
            topology.push_back(static_cast<HDF5::U64>(key.size()));

          for (size_t k = 0; k < key.size(); ++k)
            topology.push_back(static_cast<HDF5::U64>(key[k]));
        }

        HDF5::writeVectorDataset(file, HDF5::Path::MeshXDMFTopology, topology);
        HDF5::writeScalarDataset(
            file,
            HDF5::Path::MeshXDMFTopologySize,
            static_cast<HDF5::U64>(topology.size()));
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

      explicit
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction loading is file-path based."
          << Alert::Raise;
      }

      void load(const boost::filesystem::path& filename) override
      {
        auto& gf = this->getObject();
        const auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open HDF5 GridFunction file: " << filename
            << Alert::Raise;
        }

        const auto values = HDF5::readVectorDataset<HDF5::F64>(file.get(), HDF5::Path::GridFunctionValuesData);
        const auto dofCount = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::GridFunctionMetaSize));
        const auto vectorDim = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::GridFunctionMetaDimension));

        if (values.size() != static_cast<size_t>(gf.getData().size()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid GridFunction data size."
            << Alert::Raise;
        }
        if (dofCount != gf.getSize())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GridFunction size mismatch."
            << Alert::Raise;
        }
        if (vectorDim != gf.getDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GridFunction dimension mismatch."
            << Alert::Raise;
        }

        auto& data = gf.getData();
        for (size_t i = 0; i < values.size(); ++i)
          data[i] = static_cast<typename std::remove_reference_t<decltype(data[0])>>(values[i]);
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

      explicit
      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction printing is file-path based."
          << Alert::Raise;
      }

      void print(const boost::filesystem::path& filename) override
      {
        const auto& gf = this->getObject();

        const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to create HDF5 GridFunction file: " << filename
            << Alert::Raise;
        }

        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction group."
              << Alert::Raise;
          }
        }
        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction/Meta group."
              << Alert::Raise;
          }
        }
        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction/Values group."
              << Alert::Raise;
          }
        }

        const auto& data = gf.getData();
        std::vector<HDF5::F64> values(static_cast<size_t>(data.size()));
        for (size_t i = 0; i < values.size(); ++i)
          values[i] = static_cast<HDF5::F64>(data[i]);

        HDF5::writeVectorDataset(file.get(), HDF5::Path::GridFunctionValuesData, values);
        HDF5::writeScalarDataset(file.get(), HDF5::Path::GridFunctionMetaSize, static_cast<HDF5::U64>(gf.getSize()));
        HDF5::writeScalarDataset(file.get(), HDF5::Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(gf.getDimension()));
      }
  };
}

#endif
