/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POINTCLOUD_H
#define RODIN_GEOMETRY_POINTCLOUD_H

#include <cstdint>
#include <cassert>
#include <ostream>
#include <vector>

#include <Eigen/Core>

#include "Rodin/Math/SpatialVector.h"

namespace Rodin::Geometry
{
  /**
   * @brief Point cloud with a point-centric API and Eigen "matrix" views.
   *
   * Ownership: std::vector<Eigen::Vector3<Real>> (always 3 scalars per point).
   * Active dimension: m_rows in {0,1,2,3} tells which prefix is meaningful.
   *
   * Views:
   * - getMatrix(): (rows x N) view for computations in the active dimension.
   * - getPackedMatrix(): (3 x N) full view of the underlying packed storage.
   *
   * Storage layout (AoS):
   *   [x0 y0 z0][x1 y1 z1]...
   *
   * We map it as a column-major 3xN matrix with outer stride 3 and inner stride 1.
   */
  class PointCloud
  {
    public:
      /// @brief The coordinate scalar type.
      using Scalar = Real;
      /// @brief Maximum storable spatial dimension (RODIN_MAXIMAL_SPACE_DIMENSION, fixed at 3).
      static constexpr std::uint8_t MaxRows = RODIN_MAXIMAL_SPACE_DIMENSION;
      static_assert(MaxRows == 3, "RODIN_MAXIMAL_SPACE_DIMENSION must be equal to 3.");

      /// @brief Underlying packed storage (three coordinates per point).
      using Data = std::vector<std::array<Scalar, 3>>;

      // Map types with explicit stride (AoS -> column-major 3xN)
      /// @brief Eigen stride mapping the packed storage to a column-major view.
      using StrideType = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

      /// @brief Mutable Eigen map of the packed storage as a 3-by-N matrix.
      using MapType3xN = Eigen::Map<
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
        0, StrideType>;

      /// @brief Const Eigen map of the packed storage as a 3-by-N matrix.
      using ConstMapType3xN = Eigen::Map<
        const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
        0, StrideType>;

      /// @brief Constructs an empty point cloud (dimension 0, no points).
      PointCloud() noexcept
        : m_dimension(0),
          m_pts()
      {}

      /// @brief Constructs @p n zero-initialized points in @p rows dimensions.
      explicit
      PointCloud(std::uint8_t rows, size_t n)
        : m_dimension(rows),
          m_pts(static_cast<size_t>(n))
      {
        assert(rows <= MaxRows);
      }

      /// @brief Copy constructor.
      PointCloud(const PointCloud&) = default;

      /// @brief Move constructor.
      PointCloud(PointCloud&&) = default;

      /// @brief Copy assignment operator.
      PointCloud& operator=(const PointCloud&) = default;

      /// @brief Move assignment operator.
      PointCloud& operator=(PointCloud&&) = default;

      /// @brief Returns the active spatial dimension (0..3).
      [[nodiscard]] constexpr
      std::uint8_t getDimension() const noexcept
      {
        return m_dimension;
      }

      /// @brief Returns the active spatial dimension (Eigen-style row count).
      [[nodiscard]] constexpr
      std::uint8_t rows() const noexcept
      {
        return m_dimension;
      }

      /// @brief Returns the number of points (Eigen-style column count).
      [[nodiscard]] inline
      size_t cols() const noexcept
      {
        return m_pts.size();
      }

      /// @brief Returns the number of points.
      [[nodiscard]] inline
      size_t getCount() const noexcept
      {
        return m_pts.size();
      }

      /// @brief Sets the active spatial dimension (must not exceed MaxRows).
      void setDimension(std::uint8_t r) noexcept
      {
        assert(r <= MaxRows);
        m_dimension = r;
      }

      /// @brief Resizes to @p n points in @p r dimensions.
      void resize(std::uint8_t r, size_t n)
      {
        assert(r <= MaxRows);
        m_dimension = r;
        m_pts.resize(static_cast<size_t>(n));
      }

      /// @brief Reserves storage for @p n points.
      void reserve(size_t n)
      {
        m_pts.reserve(static_cast<size_t>(n));
      }

      /// @brief Removes all points.
      void clear() noexcept
      {
        m_pts.clear();
      }

      /// @brief Appends a 1D point (requires active dimension 1).
      void push_back(const std::array<Scalar, 1>& p)
      {
        assert(m_dimension == 1);
        m_pts.push_back({ p[0], 0, 0 });
      }

      /// @brief Appends a 2D point (requires active dimension 2).
      void push_back(const std::array<Real, 2>& p)
      {
        assert(m_dimension == 2);
        m_pts.push_back({ p[0], p[1], 0 });
      }

      /// @brief Appends a 3D point (requires active dimension 3).
      void push_back(const std::array<Scalar, 3>& p)
      {
        assert(m_dimension == 3);
        m_pts.push_back(p);
      }

      /// @brief Appends a spatial point, zero-padding unused coordinates.
      void push_back(const Math::SpatialPoint& p)
      {
        assert(p.size() <= MaxRows);
        std::array<Scalar, 3> pt = { 0, 0, 0 };
        switch (p.size())
        {
          case 3: pt[2] = p(2); [[fallthrough]];
          case 2: pt[1] = p(1); [[fallthrough]];
          case 1: pt[0] = p(0); break;
          case 0: break;
          default: assert(false);
        }
        m_pts.push_back(pt);
      }

      /// @brief Returns a reference to coordinate @p i of point @p j.
      [[nodiscard]] inline
      Scalar& operator()(std::uint8_t i, size_t j) noexcept
      {
        assert(i < m_dimension);
        assert(j < getCount());
        return m_pts[j][static_cast<size_t>(i)];
      }

      /// @brief Returns a const reference to coordinate @p i of point @p j.
      [[nodiscard]] inline
      const Scalar& operator()(std::uint8_t i, size_t j) const noexcept
      {
        assert(i < m_dimension);
        assert(j < getCount());
        return m_pts[j][i];
      }

      /// @brief Returns a reference to the packed 3-coordinate storage of point @p j.
      [[nodiscard]] inline
      std::array<Scalar, 3>& point3(size_t j) noexcept
      {
        assert(j < getCount());
        return m_pts[static_cast<size_t>(j)];
      }

      /// @brief Returns a const reference to the packed 3-coordinate storage of point @p j.
      [[nodiscard]] inline
      const std::array<Scalar, 3>& point3(size_t j) const noexcept
      {
        assert(j < getCount());
        return m_pts[static_cast<size_t>(j)];
      }

      /// @brief Returns point @p j as a SpatialPoint of the active dimension.
      [[nodiscard]] inline
      Math::SpatialPoint operator[](size_t j) const noexcept
      {
        assert(j < getCount());
        Math::SpatialPoint p(m_dimension);
        const auto& v = this->point3(j);
        switch (m_dimension)
        {
          case 3: p(2) = v[2]; [[fallthrough]];
          case 2: p(1) = v[1]; [[fallthrough]];
          case 1: p(0) = v[0]; break;
          case 0: break;
          default: assert(false);
        }
        return p;
      }

      /// @brief Returns point @p j as a SpatialPoint (Eigen-style column accessor).
      [[nodiscard]] inline
      auto col(size_t j) const noexcept
      {
        return this->operator[](j);
      }

      /// @brief Sets the active coordinates of all points to zero.
      void setZero() noexcept
      {
        // Fallthrough: touch only the active prefix (keeps inactive components unchanged).
        if (m_dimension == 0 || getCount() == 0)
          return;

        for (auto& p : m_pts)
        {
          switch (m_dimension)
          {
            case 3: p[2] = Scalar(0); [[fallthrough]];
            case 2: p[1] = Scalar(0); [[fallthrough]];
            case 1: p[0] = Scalar(0); break;
            case 0: break;
            default: assert(false);
          }
        }
      }

      // --- Eigen views ---

      /**
       * @brief Returns a (rows x N) matrix view of the active coordinates.
       */
      [[nodiscard]] inline
      auto getMatrix() noexcept
      {
        return getPackedMatrix().topRows(static_cast<Eigen::Index>(m_dimension));
      }

      /// @brief Returns a (rows x N) matrix view of the active coordinates.
      [[nodiscard]] inline
      auto getMatrix() const noexcept
      {
        return getPackedMatrix().topRows(static_cast<Eigen::Index>(m_dimension));
      }

      /**
       * @brief Returns a (3 x N) matrix view of the underlying packed storage.
       *
       * Naming rationale: this is the raw, packed, always-3D storage view.
       */
      [[nodiscard]] inline
      MapType3xN getPackedMatrix() noexcept
      {
        const StrideType stride(/*outer*/ MaxRows, /*inner*/ 1);
        return MapType3xN(rawData(), static_cast<Eigen::Index>(MaxRows), getCount(), stride);
      }

      /// @brief Returns a (3 x N) matrix view of the underlying packed storage.
      [[nodiscard]] inline
      ConstMapType3xN getPackedMatrix() const noexcept
      {
        const StrideType stride(/*outer*/ MaxRows, /*inner*/ 1);
        return ConstMapType3xN(rawData(), static_cast<Eigen::Index>(MaxRows), getCount(), stride);
      }

      // --- norms / dot using views ---

      /// @brief Returns the Frobenius inner product of the active coordinates with another point cloud.
      [[nodiscard]] inline
      Scalar dot(const PointCloud& other) const noexcept
      {
        assert(m_dimension == other.m_dimension);
        assert(getCount() == other.getCount());
        if (m_dimension == 0 || getCount() == 0)
          return Scalar(0);
        return (this->getMatrix().cwiseProduct(other.getMatrix())).sum();
      }

      /// @brief Returns the Frobenius inner product of the active coordinates with an Eigen matrix.
      template <class EigenDerived>
      [[nodiscard]] inline
      Scalar dot(const Eigen::MatrixBase<EigenDerived>& other) const noexcept
      {
        assert(static_cast<std::uint8_t>(other.rows()) == m_dimension);
        assert(static_cast<Eigen::Index>(other.cols()) == getCount());
        if (m_dimension == 0 || getCount() == 0)
          return Scalar(0);
        return (this->getMatrix().cwiseProduct(other)).sum();
      }

      /// @brief Returns the squared Frobenius norm of the active coordinates.
      [[nodiscard]] inline
      Scalar squaredNorm() const noexcept
      {
        if (m_dimension == 0 || getCount() == 0)
          return Scalar(0);
        return this->getMatrix().squaredNorm();
      }

      // raw container access (serialization, mesh APIs, etc.)
      /// @brief Returns a reference to the underlying packed point storage.
      [[nodiscard]] inline
      std::vector<std::array<Scalar, 3>>& getPoints() noexcept
      {
        return m_pts;
      }

      /// @brief Returns a const reference to the underlying packed point storage.
      [[nodiscard]] inline
      const std::vector<std::array<Scalar, 3>>& getPoints() const noexcept
      {
        return m_pts;
      }

      /// @brief Sets the active coordinates of all points to @p v.
      void setConstant(const Scalar& v) noexcept
      {
        // Set active prefix; keep inactive components unchanged.
        if (m_dimension == 0 || getCount() == 0)
          return;
        this->getMatrix().setConstant(v);
      }

      /// @brief Scales the active coordinates of all points by @p s in place.
      PointCloud& operator*=(const Scalar& s) noexcept
      {
        if (m_dimension == 0 || getCount() == 0)
          return *this;

        for (auto& p : m_pts)
        {
          switch (m_dimension)
          {
            case 3: p[2] *= s; [[fallthrough]];
            case 2: p[1] *= s; [[fallthrough]];
            case 1: p[0] *= s; break;
            case 0: break;
            default: assert(false);
          }
        }

        return *this;
      }

      /// @brief Serializes the point cloud (for boost::serialization).
      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_dimension;
        ar & m_pts;
      }

    private:
      std::uint8_t m_dimension;
      Data m_pts;

      [[nodiscard]] inline
      Scalar* rawData() noexcept
      {
        return m_pts.empty() ? nullptr : m_pts.front().data();
      }

      [[nodiscard]] inline
      const Scalar* rawData() const noexcept
      {
        return m_pts.empty() ? nullptr : m_pts.front().data();
      }
  };

  // -------- free operators --------

  /// @brief Scalar-times-point-cloud (scales the active coordinates).
  [[nodiscard]] inline
  PointCloud operator*(const Real& s, const PointCloud& A)
  {
    PointCloud C(A.rows(), A.getCount());
    if (A.rows() == 0 || A.getCount() == 0)
      return C;
    C.getPackedMatrix() = A.getPackedMatrix(); // copy all 3 rows
    C.getMatrix() *= s;                        // scale active rows
    return C;
  }

  /// @brief Point-cloud-times-scalar (scales the active coordinates).
  [[nodiscard]] inline
  PointCloud operator*(const PointCloud& A, const Real& s)
  {
    return s * A;
  }

  /// @brief Componentwise sum of two point clouds' active coordinates.
  [[nodiscard]] inline
  PointCloud operator+(const PointCloud& A, const PointCloud& B)
  {
    assert(A.rows() == B.rows());
    assert(A.getCount() == B.getCount());
    PointCloud C(A.rows(), A.getCount());
    if (A.rows() == 0 || A.getCount() == 0)
      return C;
    C.getPackedMatrix() = A.getPackedMatrix();
    C.getMatrix() += B.getMatrix();
    return C;
  }

  /// @brief Componentwise difference of two point clouds' active coordinates.
  [[nodiscard]] inline
  PointCloud operator-(const PointCloud& A, const PointCloud& B)
  {
    assert(A.rows() == B.rows());
    assert(A.getCount() == B.getCount());
    PointCloud C(A.rows(), A.getCount());
    if (A.rows() == 0 || A.getCount() == 0)
      return C;
    C.getPackedMatrix() = A.getPackedMatrix();
    C.getMatrix() -= B.getMatrix();
    return C;
  }

  /// @brief Product of the active-coordinate matrix and an Eigen matrix.
  template <class EigenDerived>
  [[nodiscard]] inline
  auto operator*(const PointCloud& A, const Eigen::MatrixBase<EigenDerived>& B)
  {
    assert(static_cast<Eigen::Index>(A.getCount()) == B.rows());
    return A.getMatrix() * B;
  }

  /// @brief Product of an Eigen matrix and the active-coordinate matrix.
  template <class EigenDerived>
  [[nodiscard]] inline
  auto operator*(const Eigen::MatrixBase<EigenDerived>& A, const PointCloud& B)
  {
    assert(A.cols() == static_cast<Eigen::Index>(B.rows()));
    return A * B.getMatrix();
  }

  /// @brief Streams the active-coordinate matrix to an output stream.
  inline
  std::ostream& operator<<(std::ostream& os, const PointCloud& P)
  {
    os << P.getMatrix();
    return os;
  }
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for a Geometry::PointCloud: exposes the scalar type.
  template <>
  struct Traits<Rodin::Geometry::PointCloud>
  {
      /// @brief Scalar value type.
      using ScalarType = Rodin::Real;
  };
}

#endif
