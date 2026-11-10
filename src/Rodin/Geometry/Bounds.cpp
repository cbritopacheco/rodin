/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <Eigen/LU>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Bounds.h"
#include "PolytopeTransformation.h"

namespace Rodin::Geometry
{
  /**
   * @brief Evaluates a binomial coefficient.
   * @param n Polynomial degree.
   * @param i Basis index.
   * @returns Binomial coefficient for the supplied degree and index.
   */

  Real Bounds::binomial(size_t n, size_t i)
  {
    if (i > n)
      return 0;
    if (i > n - i)
      i = n - i;
    Real out = 1;
    for (size_t p = 1; p <= i; ++p)
    {
      out *= static_cast<Real>(n + 1 - p);
      out /= static_cast<Real>(p);
    }
    return out;
  }

  /**
   * @brief Univariate Bernstein polynomial of degree @f$ n @f$ on [0, 1].
   * @param n Polynomial degree.
   * @param i Basis index.
   * @param x Evaluation coordinate or physical query point.
   * @returns Value of the indexed univariate Bernstein polynomial at the supplied coordinate.
   */
  Real Bounds::bernstein(size_t n, size_t i, Real x)
  {
    if (i > n)
      return 0;
    Real out = binomial(n, i);
    for (size_t p = 0; p < i; ++p)
      out *= x;
    for (size_t p = 0; p < n - i; ++p)
      out *= Real(1) - x;
    return out;
  }

  /**
   * @brief Multinomial coefficient @f$ k! / \prod_i \alpha_i! @f$.
   * @param k Polynomial degree.
   * @param alpha Barycentric multi-index.
   * @param count Number of active barycentric coordinates or control points.
   * @returns Multinomial coefficient of the active multi-index entries.
   */
  Real Bounds::multinomial(size_t k, const std::array<size_t, 4>& alpha, size_t count)
  {
    Real out = 1;
    for (size_t p = 2; p <= k; ++p)
      out *= static_cast<Real>(p);
    for (size_t i = 0; i < count; ++i)
    {
      for (size_t p = 2; p <= alpha[i]; ++p)
        out /= static_cast<Real>(p);
    }
    return out;
  }

  /**
   * @brief Bernstein polynomial of degree @f$ k @f$ on the reference
   * simplex of dimension @p rdim, indexed by a barycentric multi-index.
   * @param k Polynomial degree.
   * @param alpha Barycentric multi-index.
   * @param rdim Reference dimension.
   * @param r Reference evaluation point.
   * @returns Value of the indexed simplex Bernstein polynomial at the reference point.
   */
  Real Bounds::simplexBernstein(size_t k, const std::array<size_t, 4>& alpha, size_t rdim,
    const Math::SpatialPoint& r)
  {
    Real lambda = 1;
    for (size_t i = 0; i < rdim; ++i)
      lambda -= r[static_cast<Eigen::Index>(i)];
    Real out = multinomial(k, alpha, rdim + 1);
    for (size_t p = 0; p < alpha[0]; ++p)
      out *= lambda;
    for (size_t i = 0; i < rdim; ++i)
    {
      for (size_t p = 0; p < alpha[i + 1]; ++p)
        out *= r[static_cast<Eigen::Index>(i)];
    }
    return out;
  }

  /**
   * @brief Evaluates one control-point basis function of degree @p k.
   * @param g Reference-cell geometry.
   * @param k Polynomial degree.
   * @param mode Multi-index identifying the control-basis function.
   * @param r Reference evaluation point.
   * @returns Value of the selected control-point basis function at the reference point.
   */
  Real Bounds::controlBasisValue(Geometry::Polytope::Type g, size_t k,
    const std::array<size_t, 4>& mode, const Math::SpatialPoint& r)
  {
    using G = Geometry::Polytope::Type;
    switch (g)
    {
      case G::Segment:
        return bernstein(k, mode[0], r.x());
      case G::Triangle:
        return simplexBernstein(k, mode, 2, r);
      case G::Tetrahedron:
        return simplexBernstein(k, mode, 3, r);
      case G::Quadrilateral:
        return bernstein(k, mode[0], r.x()) * bernstein(k, mode[1], r.y());
      case G::Hexahedron:
        return bernstein(k, mode[0], r.x()) * bernstein(k, mode[1], r.y()) *
          bernstein(k, mode[2], r.z());
      case G::Wedge:
        return simplexBernstein(k, mode, 2, r) * bernstein(k, mode[3], r.z());
      case G::Pyramid:
      {
        // Collapsed coordinates a = x / (1 - z), b = y / (1 - z), as in
        // the pyramid element itself. The apex layer carries the whole
        // partition of unity at z = 1.
        const size_t m = mode[2];
        const size_t n = k - m;
        const Real z = r.z();
        const Real q = Real(1) - z;
        const Real bz = bernstein(k, m, z);
        if (n == 0)
          return bz;
        if (!(q > std::numeric_limits<Real>::epsilon()))
          return 0;
        return bernstein(n, mode[0], r.x() / q) * bernstein(n, mode[1], r.y() / q) * bz;
      }
      case G::Point:
        return 1;
    }
    assert(false);
    return 0;
  }

  /**
   * @brief Reference lattice and control-point basis indices of degree
   * @p k on @p g.
   *
   * The lattice is the equispaced (collapsed, on pyramids) node set of
   * the corresponding finite element, so it is unisolvent for the space
   * the geometry transformation lives in.
   * @param g Reference-cell geometry.
   * @param k Polynomial degree.
   * @param modes Output multi-indices for the control lattice.
   * @param samples Output reference points for sampling the control basis.
   */
  void Bounds::makeControlLattice(Geometry::Polytope::Type g, size_t k,
    std::vector<std::array<size_t, 4>>& modes, std::vector<Math::SpatialPoint>& samples)
  {
    using G = Geometry::Polytope::Type;
    const Real kr = static_cast<Real>(k);
    const size_t rdim = Geometry::Polytope::Traits(g).getDimension();
    auto emplace = [&](const std::array<size_t, 4>& mode,
                     std::initializer_list<Real> coords) {
      modes.push_back(mode);
      Math::SpatialPoint r;
      r.resize(static_cast<Eigen::Index>(rdim));
      Eigen::Index i = 0;
      for (Real c : coords)
        r[i++] = c;
      assert(i == static_cast<Eigen::Index>(rdim));
      samples.push_back(std::move(r));
    };

    switch (g)
    {
      case G::Segment:
      {
        for (size_t i = 0; i <= k; ++i)
          emplace({i, 0, 0, 0}, {static_cast<Real>(i) / kr});
        return;
      }
      case G::Triangle:
      {
        for (size_t a2 = 0; a2 <= k; ++a2)
        {
          for (size_t a1 = 0; a1 + a2 <= k; ++a1)
          {
            emplace({k - a1 - a2, a1, a2, 0},
              {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr});
          }
        }
        return;
      }
      case G::Tetrahedron:
      {
        for (size_t a3 = 0; a3 <= k; ++a3)
        {
          for (size_t a2 = 0; a2 + a3 <= k; ++a2)
          {
            for (size_t a1 = 0; a1 + a2 + a3 <= k; ++a1)
            {
              emplace({k - a1 - a2 - a3, a1, a2, a3},
                {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                  static_cast<Real>(a3) / kr});
            }
          }
        }
        return;
      }
      case G::Quadrilateral:
      {
        for (size_t j = 0; j <= k; ++j)
        {
          for (size_t i = 0; i <= k; ++i)
          {
            emplace({i, j, 0, 0}, {static_cast<Real>(i) / kr, static_cast<Real>(j) / kr});
          }
        }
        return;
      }
      case G::Hexahedron:
      {
        for (size_t l = 0; l <= k; ++l)
        {
          for (size_t j = 0; j <= k; ++j)
          {
            for (size_t i = 0; i <= k; ++i)
            {
              emplace({i, j, l, 0},
                {static_cast<Real>(i) / kr, static_cast<Real>(j) / kr,
                  static_cast<Real>(l) / kr});
            }
          }
        }
        return;
      }
      case G::Wedge:
      {
        for (size_t m = 0; m <= k; ++m)
        {
          for (size_t a2 = 0; a2 <= k; ++a2)
          {
            for (size_t a1 = 0; a1 + a2 <= k; ++a1)
            {
              emplace({k - a1 - a2, a1, a2, m},
                {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                  static_cast<Real>(m) / kr});
            }
          }
        }
        return;
      }
      case G::Pyramid:
      {
        for (size_t m = 0; m <= k; ++m)
        {
          const size_t n = k - m;
          const Real z = static_cast<Real>(m) / kr;
          const Real q = Real(1) - z;
          if (n == 0)
          {
            emplace({0, 0, m, 0}, {Real(0), Real(0), z});
            continue;
          }
          const Real nr = static_cast<Real>(n);
          for (size_t j = 0; j <= n; ++j)
          {
            for (size_t i = 0; i <= n; ++i)
            {
              emplace({i, j, m, 0},
                {q * static_cast<Real>(i) / nr, q * static_cast<Real>(j) / nr, z});
            }
          }
        }
        return;
      }
      case G::Point:
        return;
    }
    assert(false);
  }

  /**
   * @brief Reference samples and the map taking sampled coordinates to
   * control points of a geometry transformation.
   *
   * The transformation of a polytope of per-factor degree @f$ k @f$ lies
   * in the span of a basis @f$ \psi_n @f$ that is non-negative and sums
   * to one on the reference domain, using Bernstein polynomials on simplices,
   * their tensor products on quadrilaterals, hexahedra and wedges, and
   * the collapsed-coordinate product basis on pyramids. Writing
   * @f[
   *   T(r) = \sum_n C_n \, \psi_n(r),
   * @f]
   * every mapped point is a convex combination of the control points
   * @f$ C_n @f$, so their componentwise extrema bound the entire curved
   * image in exact arithmetic. Floating-point conversion requires a
   * roundoff allowance; ill-conditioned conversions cannot certify a box.
   *
   * The control points are recovered from values sampled on a unisolvent
   * lattice. The relation @f$ F_m = T(r_m) = \sum_n C_n \psi_n(r_m) @f$ reads
   * @f$ F = C V^\top @f$ with @f$ V_{mn} = \psi_n(r_m) @f$, so
   * @f$ C = F (V^{-1})^\top @f$. Only the conversion depends on the
   * geometry and the degree, so it is built once and shared by every
   * polytope of that kind.
   */
  /**
   * @brief Returns the cached control-point basis of degree @p k on @p g.
   * @param g Reference-cell geometry.
   * @param k Polynomial degree.
   * @returns Reference to the cached control-point basis for the requested geometry and degree.
   */
  const Bounds::ControlBasis& Bounds::getControlBasis(
    Geometry::Polytope::Type g, size_t k)
  {
    using Key = std::pair<int, size_t>;
    static std::map<Key, ControlBasis> s_cache;
    static std::mutex s_mutex;

    const Key key{static_cast<int>(g), k};
    std::lock_guard<std::mutex> lock(s_mutex);
    const auto it = s_cache.find(key);
    if (it != s_cache.end())
      return it->second;

    ControlBasis basis;
    std::vector<std::array<size_t, 4>> modes;
    makeControlLattice(g, k, modes, basis.samples);
    const size_t n = modes.size();
    using G = Geometry::Polytope::Type;
    basis.tensorFactors = k >= MinSeparableTensorDegree
      ? (g == G::Quadrilateral ? 2 : (g == G::Hexahedron ? 3 : 0))
      : 0;
    if (basis.tensorFactors > 0)
    {
      // Tensor conversion is a succession of small one-dimensional solves,
      // not an inverse of the full (k+1)^d by (k+1)^d lattice matrix.
      Math::Matrix<Real> vandermonde(k + 1, k + 1);
      for (size_t m = 0; m <= k; ++m)
      {
        for (size_t c = 0; c <= k; ++c)
          vandermonde(m, c) = bernstein(k, c, static_cast<Real>(m) / k);
      }
      const Eigen::FullPivLU<Math::Matrix<Real>> lu(vandermonde);
      if (lu.isInvertible())
      {
        basis.conversion = lu.inverse().transpose();
        basis.roundoffAmplification = ControlRoundoffFactor *
          static_cast<Real>((k + 1) * basis.tensorFactors) *
          std::numeric_limits<Real>::epsilon() /
          std::pow(lu.rcond(), static_cast<Real>(basis.tensorFactors));
        basis.valid = basis.conversion.allFinite() &&
          std::isfinite(basis.roundoffAmplification) &&
          basis.roundoffAmplification < Real(1);
      }
      return s_cache.emplace(key, std::move(basis)).first->second;
    }
    if (n > 0)
    {
      Math::Matrix<Real> vandermonde(n, n);
      for (size_t m = 0; m < n; ++m)
      {
        for (size_t c = 0; c < n; ++c)
        {
          vandermonde(m, c) = controlBasisValue(g, k, modes[c], basis.samples[m]);
        }
      }
      const Eigen::FullPivLU<Math::Matrix<Real>> lu(vandermonde);
      if (lu.isInvertible())
      {
        basis.conversion = lu.inverse().transpose();
        basis.roundoffAmplification = ControlRoundoffFactor * static_cast<Real>(n) *
          std::numeric_limits<Real>::epsilon() / lu.rcond();
        basis.valid = basis.conversion.allFinite() &&
          std::isfinite(basis.roundoffAmplification) &&
          basis.roundoffAmplification < Real(1);
      }
    }
    return s_cache.emplace(key, std::move(basis)).first->second;
  }

  Bounds::Bounds(
    const PolytopeTransformation& transformation, Polytope::Type geometry, size_t degree)
    : m_count(0),
      m_lower{},
      m_upper{},
      m_error(0),
      m_scale(0),
      m_dimension(transformation.getPhysicalDimension())
  {
    m_lower.fill(-std::numeric_limits<Real>::infinity());
    m_upper.fill(std::numeric_limits<Real>::infinity());
    const Polytope::Traits traits(geometry);
    if (degree <= 1)
    {
      m_count = traits.getVertexCount();
      std::array<Real, Math::SpatialVector<Real>::MaxSize> lower;
      std::array<Real, Math::SpatialVector<Real>::MaxSize> upper;
      lower.fill(std::numeric_limits<Real>::infinity());
      upper.fill(-std::numeric_limits<Real>::infinity());
      for (size_t j = 0; j < m_count; ++j)
      {
        transformation.transform(m_vertices[j], traits.getVertex(j));
        for (size_t i = 0; i < m_dimension; ++i)
        {
          const Real value = m_vertices[j][i];
          if (!std::isfinite(value))
          {
            m_count = 0;
            return;
          }
          m_scale = std::max(m_scale, std::abs(value));
          lower[i] = std::min(lower[i], value);
          upper[i] = std::max(upper[i], value);
        }
      }
      const Real error =
        ControlRoundoffFactor * std::numeric_limits<Real>::epsilon() * m_scale;
      for (size_t i = 0; i < m_dimension; ++i)
      {
        if (std::isfinite(lower[i]) && std::isfinite(upper[i]))
        {
          m_lower[i] =
            std::nextafter(lower[i] - error, -std::numeric_limits<Real>::infinity());
          m_upper[i] =
            std::nextafter(upper[i] + error, std::numeric_limits<Real>::infinity());
        }
      }
      return;
    }
    else
    {
      const auto& basis = getControlBasis(geometry, degree);
      if (!basis.valid)
        return;
      const size_t n = basis.samples.size();
      Math::Matrix<Real> sampled(m_dimension, n);
      Math::SpatialPoint point;
      for (size_t j = 0; j < n; ++j)
      {
        transformation.transform(point, basis.samples[j]);
        for (size_t i = 0; i < m_dimension; ++i)
          sampled(i, j) = point[i];
      }
      const Math::Vector<Real> origin = sampled.col(0);
      sampled.colwise() -= origin;
      if (basis.tensorFactors == 0)
        m_control = sampled * basis.conversion;
      else
      {
        m_control = sampled;
        Math::Matrix<Real> next(m_dimension, n);
        const size_t width = degree + 1;
        size_t stride = 1;
        for (size_t axis = 0; axis < basis.tensorFactors; ++axis)
        {
          for (size_t block = 0; block < n; block += width * stride)
          {
            for (size_t offset = 0; offset < stride; ++offset)
            {
              for (size_t j = 0; j < width; ++j)
              {
                const size_t dst = block + offset + j * stride;
                next.col(dst).setZero();
                for (size_t c = 0; c < width; ++c)
                {
                  next.col(dst) +=
                    m_control.col(block + offset + c * stride) * basis.conversion(c, j);
                }
              }
            }
          }
          m_control.swap(next);
          stride *= width;
        }
      }
      m_error = basis.roundoffAmplification * sampled.cwiseAbs().maxCoeff();
      m_control.colwise() += origin;
      if (!m_control.allFinite() || !std::isfinite(m_error))
      {
        m_control.resize(0, 0);
        return;
      }
      m_scale = m_control.cwiseAbs().maxCoeff();
    }
    for (size_t j = 0; j < m_count; ++j)
    {
      for (size_t i = 0; i < m_dimension; ++i)
      {
        if (!std::isfinite(m_vertices[j][i]))
        {
          m_count = 0;
          return;
        }
        m_scale = std::max(m_scale, std::abs(m_vertices[j][i]));
      }
    }
    const Real error =
      m_error + ControlRoundoffFactor * std::numeric_limits<Real>::epsilon() * m_scale;
    for (size_t i = 0; i < m_dimension; ++i)
    {
      Real lower = std::numeric_limits<Real>::infinity();
      Real upper = -std::numeric_limits<Real>::infinity();
      for (size_t j = 0; j < m_count; ++j)
      {
        lower = std::min(lower, m_vertices[j][i]);
        upper = std::max(upper, m_vertices[j][i]);
      }
      for (Eigen::Index j = 0; j < m_control.cols(); ++j)
      {
        lower = std::min(lower, m_control(i, j));
        upper = std::max(upper, m_control(i, j));
      }
      if (std::isfinite(lower) && std::isfinite(upper))
      {
        m_lower[i] =
          std::nextafter(lower - error, -std::numeric_limits<Real>::infinity());
        m_upper[i] = std::nextafter(upper + error, std::numeric_limits<Real>::infinity());
      }
    }
  }

  Real Bounds::getDirectionalSupport(
    const Math::SpatialVector<Real>& direction, Real length) const
  {
    Real upper = -std::numeric_limits<Real>::infinity();
    for (size_t j = 0; j < m_count; ++j)
      upper = std::max(upper, m_vertices[j].dot(direction));
    for (Eigen::Index j = 0; j < m_control.cols(); ++j)
    {
      Real product = 0;
      for (size_t i = 0; i < m_dimension; ++i)
        product += m_control(i, j) * direction[i];
      upper = std::max(upper, product);
    }
    if (!std::isfinite(upper))
      return std::numeric_limits<Real>::infinity();
    const Real error = length *
      (m_error + ControlRoundoffFactor * std::numeric_limits<Real>::epsilon() * m_scale);
    return std::nextafter(upper + error, std::numeric_limits<Real>::infinity());
  }
}
