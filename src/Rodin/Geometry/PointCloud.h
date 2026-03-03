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
      using Scalar = Real;
      static constexpr std::uint8_t MaxRows = RODIN_MAXIMAL_SPACE_DIMENSION;
      static_assert(MaxRows == 3, "RODIN_MAXIMAL_SPACE_DIMENSION must be equal to 3.");

      using Data = std::vector<std::array<Scalar, 3>>;

      // Map types with explicit stride (AoS -> column-major 3xN)
      using StrideType = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

      using MapType3xN = Eigen::Map<
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
        0, StrideType>;

      using ConstMapType3xN = Eigen::Map<
        const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
        0, StrideType>;

      PointCloud() noexcept
        : m_dimension(0),
          m_pts()
      {}

      explicit
      PointCloud(std::uint8_t rows, size_t n)
        : m_dimension(rows),
          m_pts(static_cast<size_t>(n))
      {
        assert(rows <= MaxRows);
      }

      PointCloud(const PointCloud&) = default;

      PointCloud(PointCloud&&) = default;

      PointCloud& operator=(const PointCloud&) = default;

      PointCloud& operator=(PointCloud&&) = default;

      [[nodiscard]] constexpr
      std::uint8_t getDimension() const noexcept
      {
        return m_dimension;
      }

      [[nodiscard]] constexpr
      std::uint8_t rows() const noexcept
      {
        return m_dimension;
      }

      [[nodiscard]] inline
      size_t cols() const noexcept
      {
        return m_pts.size();
      }

      [[nodiscard]] inline
      size_t getCount() const noexcept
      {
        return m_pts.size();
      }

      void setDimension(std::uint8_t r) noexcept
      {
        assert(r <= MaxRows);
        m_dimension = r;
      }

      void resize(std::uint8_t r, size_t n)
      {
        assert(r <= MaxRows);
        m_dimension = r;
        m_pts.resize(static_cast<size_t>(n));
      }

      void reserve(size_t n)
      {
        m_pts.reserve(static_cast<size_t>(n));
      }

      void clear() noexcept
      {
        m_pts.clear();
      }

      void push_back(const std::array<Scalar, 1>& p)
      {
        assert(m_dimension == 1);
        m_pts.push_back({ p[0], 0, 0 });
      }

      void push_back(const std::array<Real, 2>& p)
      {
        assert(m_dimension == 2);
        m_pts.push_back({ p[0], p[1], 0 });
      }

      void push_back(const std::array<Scalar, 3>& p)
      {
        assert(m_dimension == 3);
        m_pts.push_back(p);
      }

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

      [[nodiscard]] inline
      Scalar& operator()(std::uint8_t i, size_t j) noexcept
      {
        assert(i < m_dimension);
        assert(j < getCount());
        return m_pts[j][static_cast<size_t>(i)];
      }

      [[nodiscard]] inline
      const Scalar& operator()(std::uint8_t i, size_t j) const noexcept
      {
        assert(i < m_dimension);
        assert(j < getCount());
        return m_pts[j][i];
      }

      [[nodiscard]] inline
      std::array<Scalar, 3>& point3(size_t j) noexcept
      {
        assert(j < getCount());
        return m_pts[static_cast<size_t>(j)];
      }

      [[nodiscard]] inline
      const std::array<Scalar, 3>& point3(size_t j) const noexcept
      {
        assert(j < getCount());
        return m_pts[static_cast<size_t>(j)];
      }

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

      [[nodiscard]] inline
      auto col(size_t j) const noexcept
      {
        return this->operator[](j);
      }

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

      [[nodiscard]] inline
      ConstMapType3xN getPackedMatrix() const noexcept
      {
        const StrideType stride(/*outer*/ MaxRows, /*inner*/ 1);
        return ConstMapType3xN(rawData(), static_cast<Eigen::Index>(MaxRows), getCount(), stride);
      }

      // --- norms / dot using views ---

      [[nodiscard]] inline
      Scalar dot(const PointCloud& other) const noexcept
      {
        assert(m_dimension == other.m_dimension);
        assert(getCount() == other.getCount());
        if (m_dimension == 0 || getCount() == 0)
          return Scalar(0);
        return (this->getMatrix().cwiseProduct(other.getMatrix())).sum();
      }

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

      [[nodiscard]] inline
      Scalar squaredNorm() const noexcept
      {
        if (m_dimension == 0 || getCount() == 0)
          return Scalar(0);
        return this->getMatrix().squaredNorm();
      }

      // raw container access (serialization, mesh APIs, etc.)
      [[nodiscard]] inline
      std::vector<std::array<Scalar, 3>>& getPoints() noexcept
      {
        return m_pts;
      }

      [[nodiscard]] inline
      const std::vector<std::array<Scalar, 3>>& getPoints() const noexcept
      {
        return m_pts;
      }

      void setConstant(const Scalar& v) noexcept
      {
        // Set active prefix; keep inactive components unchanged.
        if (m_dimension == 0 || getCount() == 0)
          return;
        this->getMatrix().setConstant(v);
      }

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

  [[nodiscard]] inline
  PointCloud operator*(const PointCloud& A, const Real& s)
  {
    return s * A;
  }

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

  template <class EigenDerived>
  [[nodiscard]] inline
  auto operator*(const PointCloud& A, const Eigen::MatrixBase<EigenDerived>& B)
  {
    assert(static_cast<Eigen::Index>(A.getCount()) == B.rows());
    return A.getMatrix() * B;
  }

  template <class EigenDerived>
  [[nodiscard]] inline
  auto operator*(const Eigen::MatrixBase<EigenDerived>& A, const PointCloud& B)
  {
    assert(A.cols() == static_cast<Eigen::Index>(B.rows()));
    return A * B.getMatrix();
  }

  inline
  std::ostream& operator<<(std::ostream& os, const PointCloud& P)
  {
    os << P.getMatrix();
    return os;
  }
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<Rodin::Geometry::PointCloud>
  {
    using ScalarType = Rodin::Real;
  };
}

#endif
