/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_H1_H
#define RODIN_KOKKOS_VARIATIONAL_H1_H

/**
 * @file
 * @brief H1 Lagrange spaces on a resident Kokkos mesh.
 */

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Rodin/Variational/H1/H1.h"
#include "Rodin/Variational/H1/Div.h"
#include "Rodin/Variational/H1/Grad.h"
#include "Rodin/Variational/H1/Jacobian.h"
#include "Rodin/Variational/H1/ShapeFunction.h"

#include "Rodin/Kokkos/Geometry/Mesh.h"

#include "DOFs.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Quadratic Lagrange basis of the simplices.
   *
   * Node @f$ q @f$ of a simplex follows the node order of
   * Rodin::Variational::H1Element of degree two. It is either the vertex
   * @f$ i @f$ or the midpoint of the edge @f$ (i, j) @f$, and its basis
   * function is a product of barycentric coordinates
   * @f[
   *   \phi_q = \lambda_i (2 \lambda_i - 1) \quad \text{or} \quad
   *   \phi_q = 4 \lambda_i \lambda_j .
   * @f]
   * The barycentric coordinates of a simplex are its P1 basis, so the
   * values and physical gradients are combinations of the P1 basis cached
   * by IntegrationPointView.
   */
  class QuadraticBasis
  {
    public:
      /// @brief Maximum number of nodes, reached by the tetrahedron.
      static constexpr Index Capacity = 10;

      /**
       * @brief Gets the basis of every simplex.
       *
       * The table is built once on the host from the reference nodes of the
       * host elements.
       *
       * @return Quadratic basis of the point, segment, triangle and
       * tetrahedron.
       */
      static const QuadraticBasis& get()
      {
        static const QuadraticBasis s_basis = [] {
          using Type = ::Rodin::Geometry::Polytope::Type;
          constexpr std::array<Type, 4> simplices = {
            Type::Point, Type::Segment, Type::Triangle, Type::Tetrahedron};
          QuadraticBasis basis;
          for (size_t d = 0; d < simplices.size(); ++d)
          {
            const auto& nodes =
              ::Rodin::Variational::H1Element<2, Real>::getNodes(simplices[d]);
            assert(nodes.size() == getCount(d));
            for (size_t q = 0; q < nodes.size(); ++q)
            {
              std::array<Real, 4> lambda = {1, 0, 0, 0};
              for (size_t k = 0; k < d; ++k)
              {
                lambda[k + 1] = nodes[q](k);
                lambda[0] -= nodes[q](k);
              }
              std::array<std::uint8_t, 2> pair = {0, 0};
              size_t found = 0;
              for (size_t a = 0; a <= d; ++a)
              {
                if (lambda[a] > 0.25)
                  pair[found++] = static_cast<std::uint8_t>(a);
              }
              assert(found == 1 || found == 2);
              basis.m_first[d][q] = pair[0];
              basis.m_second[d][q] = found == 1 ? pair[0] : pair[1];
            }
          }
          return basis;
        }();
        return s_basis;
      }

      /**
       * @brief Gets the number of nodes of a simplex.
       * @param d Dimension of the simplex.
       * @return @f$ (d + 1)(d + 2) / 2 @f$.
       */
      KOKKOS_INLINE_FUNCTION
      static constexpr Index getCount(Index d)
      {
        return (d + 1) * (d + 2) / 2;
      }

      /**
       * @brief Evaluates a basis function.
       * @param ip Integration point on a simplex.
       * @param q Local node index.
       * @return Basis value @f$ \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
      {
        const Index d = ip.getPoint().getPolytope().getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        const Index i = m_first[d][q];
        const Index j = m_second[d][q];
        const Real li = ip.getBasis(i);
        if (i == j)
          return li * (2 * li - 1);
        return 4 * li * ip.getBasis(j);
      }

      /**
       * @brief Evaluates the physical gradient of a basis function.
       * @param ip Integration point on a simplex.
       * @param q Local node index.
       * @return Gradient @f$ \nabla \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(const Point& ip, Index q) const
      {
        const Index d = ip.getPoint().getPolytope().getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        const Index i = m_first[d][q];
        const Index j = m_second[d][q];
        const Real li = ip.getBasis(i);
        const auto& gi = ip.getGradient(i);
        Geometry::Vector result(gi.size());
        if (i == j)
        {
          for (Index c = 0; c < gi.size(); ++c)
            result(c) = (4 * li - 1) * gi(c);
          return result;
        }
        const Real lj = ip.getBasis(j);
        const auto& gj = ip.getGradient(j);
        for (Index c = 0; c < gi.size(); ++c)
          result(c) = 4 * (lj * gi(c) + li * gj(c));
        return result;
      }

      /**
       * @brief Gets the reference coordinates of a node.
       * @param polytope Simplex of the mesh.
       * @param q Local node index.
       * @return Vertex @f$ i @f$ or midpoint of the edge @f$ (i, j) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
      {
        const Index d = polytope.getDimension();
        KOKKOS_ASSERT(q < getCount(d));
        Geometry::Vector result = getReferenceVertex(polytope, m_first[d][q]);
        const Geometry::Vector other = getReferenceVertex(polytope, m_second[d][q]);
        for (Index k = 0; k < result.size(); ++k)
          result(k) = (result(k) + other(k)) / 2;
        return result;
      }

    private:
      ::Kokkos::Array<::Kokkos::Array<std::uint8_t, Capacity>, 4> m_first;
      ::Kokkos::Array<::Kokkos::Array<std::uint8_t, Capacity>, 4> m_second;
  };

  /**
   * @brief Lagrange basis of degree K on the device.
   *
   * The basis function of node @f$ q @f$ of a reference geometry is
   * expanded in monomials of the reference coordinates,
   * @f[
   *   \phi_q(r) = \sum_m C_{qm} r^{\alpha_m} ,
   * @f]
   * over @f$ P_K @f$ on the simplices, @f$ Q_K @f$ on the quadrilateral and
   * the hexahedron, and @f$ P_K \otimes P_K @f$ on the wedge, so that it
   * spans the space of Rodin::Variational::H1Element. The coefficients
   * @f$ C = V^{-T} @f$ invert the Vandermonde matrix @f$ V_{im} @f$ of the
   * monomials at the host nodes, computed once on the host in extended
   * precision.
   *
   * The nodes of the segment, the quadrilateral and the hexahedron are the
   * tensor products of the nodes @f$ x_0, \dots, x_K @f$ of the segment, so
   * their basis functions factor into the one dimensional Lagrange
   * polynomials
   * @f[
   *   \phi_q(r) = \prod_k \ell_{i_k}(r_k) , \quad
   *   \ell_i(x) = w_i \prod_{m \neq i} (x - x_m) ,
   * @f]
   * where @f$ (i_k) @f$ are the indices of the coordinates of node
   * @f$ q @f$ and @f$ w_i @f$ the barycentric weights. They are evaluated in
   * this form, which takes a fixed number of operations whatever the node.
   *
   * The physical gradient follows from the reference one through
   * @f$ \nabla r_k = \sum_a \hat{x}_{a,k} \nabla \lambda_a @f$, where
   * @f$ \lambda_a @f$ is the P1 basis cached by IntegrationPointView and
   * @f$ \hat{x}_a @f$ its reference vertex, the P1 basis reproducing the
   * reference coordinates. Pyramids are not supported.
   *
   * @tparam K Polynomial degree.
   * @tparam MemorySpace Kokkos memory space of the tables.
   */
  template <size_t K, class MemorySpace>
  class LagrangeBasis
  {
    public:
      /// @brief Number of reference geometries.
      static constexpr Index GeometryCount = 8;

      /**
       * @brief Gets the number of nodes of a reference geometry.
       * @param g Reference geometry.
       * @return Number of nodes, zero for the unsupported pyramid.
       */
      KOKKOS_INLINE_FUNCTION
      static constexpr Index getCount(::Rodin::Geometry::Polytope::Type g)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        constexpr Index n = K + 1;
        switch (g)
        {
          case Type::Point:
            return 1;
          case Type::Segment:
            return n;
          case Type::Triangle:
            return n * (n + 1) / 2;
          case Type::Quadrilateral:
            return n * n;
          case Type::Tetrahedron:
            return n * (n + 1) * (n + 2) / 6;
          case Type::Hexahedron:
            return n * n * n;
          case Type::Wedge:
            return n * n * (n + 1) / 2;
          case Type::Pyramid:
            return 0;
        }
        return 0;
      }

      /// @brief Empty tables.
      LagrangeBasis() = default;

      /**
       * @brief Builds the tables of every supported geometry.
       * @param execution Execution space ordering the transfers.
       */
      template <class ExecutionSpace>
      explicit LagrangeBasis(const ExecutionSpace& execution)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        Index nodes = 0;
        Index coefficients = 0;
        for (Index g = 0; g < GeometryCount; ++g)
        {
          const Index n = getCount(static_cast<Type>(g));
          m_count[g] = n;
          m_nodeOffsets[g] = nodes;
          m_coefficientOffsets[g] = coefficients;
          nodes += n;
          coefficients += n * n;
        }
        m_exponents = ::Kokkos::View<std::uint8_t*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin H1 monomials")), 3 * nodes);
        m_nodes = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin H1 nodes")), 3 * nodes);
        m_coefficients = ::Kokkos::View<Real*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin H1 coefficients")),
          coefficients);
        m_lines = ::Kokkos::View<std::uint8_t*, MemorySpace>(
          ::Kokkos::view_alloc(execution, std::string("Rodin H1 line indices")),
          3 * nodes);
        auto exponents = ::Kokkos::create_mirror_view(m_exponents);
        auto points = ::Kokkos::create_mirror_view(m_nodes);
        auto values = ::Kokkos::create_mirror_view(m_coefficients);
        auto lines = ::Kokkos::create_mirror_view(m_lines);
        const auto& segment =
          ::Rodin::Variational::H1Element<K, Real>::getNodes(Type::Segment);
        assert(segment.size() == K + 1);
        for (Index i = 0; i <= K; ++i)
          m_line[i] = segment[i](0);
        for (Index i = 0; i <= K; ++i)
        {
          long double product = 1;
          for (Index m = 0; m <= K; ++m)
          {
            if (m != i)
              product *= static_cast<long double>(segment[i](0)) - segment[m](0);
          }
          m_weights[i] = static_cast<Real>(1 / product);
        }
        for (Index g = 0; g < GeometryCount; ++g)
        {
          const Type type = static_cast<Type>(g);
          const Index n = m_count[g];
          if (n == 0)
            continue;
          const auto alpha = getExponents(type);
          const auto& reference =
            ::Rodin::Variational::H1Element<K, Real>::getNodes(type);
          assert(alpha.size() == n && reference.size() == n);
          // Transpose of the Vandermonde matrix, row-major.
          std::vector<long double> vandermonde(n * n);
          for (Index i = 0; i < n; ++i)
          {
            const Index rdim = reference[i].size();
            for (Index k = 0; k < 3; ++k)
            {
              points(3 * (m_nodeOffsets[g] + i) + k) =
                k < rdim ? reference[i](k) : Real(0);
            }
            for (Index m = 0; m < n; ++m)
            {
              long double monomial = 1;
              for (Index k = 0; k < rdim; ++k)
              {
                for (Index p = 0; p < alpha[m][k]; ++p)
                  monomial *= static_cast<long double>(reference[i](k));
              }
              vandermonde[m * n + i] = monomial;
            }
          }
          const auto inverse = invert(std::move(vandermonde), n);
          for (Index m = 0; m < n; ++m)
          {
            for (Index k = 0; k < 3; ++k)
              exponents(3 * (m_nodeOffsets[g] + m) + k) = alpha[m][k];
          }
          // A tensor geometry factors when each node coordinate is a segment node.
          m_tensor[g] = type == Type::Segment || type == Type::Quadrilateral ||
            type == Type::Hexahedron;
          for (Index i = 0; i < n && m_tensor[g]; ++i)
          {
            const Index rdim = reference[i].size();
            for (Index k = 0; k < 3; ++k)
            {
              Index match = 0;
              if (k < rdim)
              {
                match = K + 1;
                for (Index m = 0; m <= K; ++m)
                {
                  if (std::abs(reference[i](k) - m_line[m]) < 1e-12)
                    match = m;
                }
                if (match > K)
                  m_tensor[g] = false;
              }
              lines(3 * (m_nodeOffsets[g] + i) + k) = static_cast<std::uint8_t>(match);
            }
          }
          for (Index q = 0; q < n; ++q)
          {
            for (Index m = 0; m < n; ++m)
            {
              values(m_coefficientOffsets[g] + q * n + m) =
                static_cast<Real>(inverse[q * n + m]);
            }
          }
        }
        ::Kokkos::deep_copy(execution, m_exponents, exponents);
        ::Kokkos::deep_copy(execution, m_nodes, points);
        ::Kokkos::deep_copy(execution, m_coefficients, values);
        ::Kokkos::deep_copy(execution, m_lines, lines);
        execution.fence("Rodin H1 basis transfer");
      }

      /**
       * @brief Evaluates a basis function.
       * @param ip Integration point.
       * @param q Local node index.
       * @return Basis value @f$ \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
      {
        const Geometry::PointView& p = ip.getPoint();
        const Index g = static_cast<Index>(p.getPolytope().getGeometry());
        const Geometry::Vector& r = p.getReferenceCoordinates();
        const Index n = m_count[g];
        KOKKOS_ASSERT(q < n);
        if (m_tensor[g])
        {
          const std::uint8_t* line = m_lines.data() + 3 * (m_nodeOffsets[g] + q);
          Real result = 1;
          for (Index k = 0; k < r.size(); ++k)
            result *= getLine(line[k], r(k)).value;
          return result;
        }
        const Powers powers = getPowers(r);
        const Real* c = m_coefficients.data() + m_coefficientOffsets[g] + q * n;
        const std::uint8_t* alpha = m_exponents.data() + 3 * m_nodeOffsets[g];
        Real result = 0;
        for (Index m = 0; m < n; ++m)
        {
          Real monomial = 1;
          for (Index k = 0; k < r.size(); ++k)
            monomial *= powers[k][alpha[3 * m + k]];
          result += c[m] * monomial;
        }
        return result;
      }

      /**
       * @brief Evaluates the physical gradient of a basis function.
       * @param ip Integration point.
       * @param q Local node index.
       * @return Gradient @f$ \nabla \phi_q @f$.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(const Point& ip, Index q) const
      {
        const Geometry::PointView& p = ip.getPoint();
        const auto type = p.getPolytope().getGeometry();
        const Index g = static_cast<Index>(type);
        const Geometry::Vector& r = p.getReferenceCoordinates();
        const Index n = m_count[g];
        KOKKOS_ASSERT(q < n);
        const Index rdim = r.size();
        ::Kokkos::Array<Real, 3> reference = {0, 0, 0};
        if (m_tensor[g])
        {
          const std::uint8_t* line = m_lines.data() + 3 * (m_nodeOffsets[g] + q);
          ::Kokkos::Array<Line, 3> factors;
          for (Index k = 0; k < rdim; ++k)
            factors[k] = getLine(line[k], r(k));
          for (Index k = 0; k < rdim; ++k)
          {
            Real derivative = factors[k].derivative;
            for (Index l = 0; l < rdim; ++l)
              derivative *= l == k ? Real(1) : factors[l].value;
            reference[k] = derivative;
          }
          return map(ip, type, reference);
        }
        const Powers powers = getPowers(r);
        const Real* c = m_coefficients.data() + m_coefficientOffsets[g] + q * n;
        const std::uint8_t* alpha = m_exponents.data() + 3 * m_nodeOffsets[g];
        for (Index m = 0; m < n; ++m)
        {
          for (Index k = 0; k < rdim; ++k)
          {
            const Index e = alpha[3 * m + k];
            if (e == 0)
              continue;
            Real derivative = c[m] * e * powers[k][e - 1];
            for (Index l = 0; l < rdim; ++l)
            {
              if (l != k)
                derivative *= powers[l][alpha[3 * m + l]];
            }
            reference[k] += derivative;
          }
        }
        return map(ip, type, reference);
      }

      /**
       * @brief Gets the reference coordinates of a node.
       * @param polytope Polytope of the mesh.
       * @param q Local node index.
       * @return Node of Rodin::Variational::H1Element.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
      {
        const Index g = static_cast<Index>(polytope.getGeometry());
        KOKKOS_ASSERT(q < m_count[g]);
        const Index rdim = polytope.getDimension();
        Geometry::Vector result(rdim);
        for (Index k = 0; k < rdim; ++k)
          result(k) = m_nodes(3 * (m_nodeOffsets[g] + q) + k);
        return result;
      }

    private:
      using Powers = ::Kokkos::Array<::Kokkos::Array<Real, K + 1>, 3>;

      /**
       * @brief Inverts a matrix by Gauss-Jordan elimination with partial
       * pivoting.
       * @param a Row-major square matrix.
       * @param n Order of the matrix.
       * @return Row-major inverse.
       *
       * Eigen marks its functions host and device under nvcc, which then warns
       * that its long double code calls host functions and loses precision on
       * the device, although the basis is built on the host only.
       */
      static std::vector<long double> invert(std::vector<long double> a, Index n)
      {
        std::vector<long double> inverse(n * n, 0);
        for (Index i = 0; i < n; ++i)
          inverse[i * n + i] = 1;
        for (Index j = 0; j < n; ++j)
        {
          Index pivot = j;
          for (Index i = j + 1; i < n; ++i)
          {
            if (std::abs(a[i * n + j]) > std::abs(a[pivot * n + j]))
              pivot = i;
          }
          assert(a[pivot * n + j] != 0);
          if (pivot != j)
          {
            for (Index k = 0; k < n; ++k)
            {
              std::swap(a[j * n + k], a[pivot * n + k]);
              std::swap(inverse[j * n + k], inverse[pivot * n + k]);
            }
          }
          const long double scale = 1 / a[j * n + j];
          for (Index k = 0; k < n; ++k)
          {
            a[j * n + k] *= scale;
            inverse[j * n + k] *= scale;
          }
          for (Index i = 0; i < n; ++i)
          {
            const long double factor = a[i * n + j];
            if (i == j || factor == 0)
              continue;
            for (Index k = 0; k < n; ++k)
            {
              a[i * n + k] -= factor * a[j * n + k];
              inverse[i * n + k] -= factor * inverse[j * n + k];
            }
          }
        }
        return inverse;
      }

      /**
       * @brief Maps a reference gradient to the physical polytope.
       * @param ip Integration point.
       * @param type Reference geometry of the polytope.
       * @param reference Reference gradient.
       * @return Physical gradient.
       */
      template <class Point>
      KOKKOS_INLINE_FUNCTION static Geometry::Vector map(const Point& ip,
        ::Rodin::Geometry::Polytope::Type type, const ::Kokkos::Array<Real, 3>& reference)
      {
        const Index rdim = ip.getPoint().getReferenceCoordinates().size();
        // The P1 basis reproduces the reference coordinates.
        const Index sdim = ip.getGradient(0).size();
        Geometry::Vector result(sdim);
        for (Index a = 0; a < ip.getCount(); ++a)
        {
          Real weight = 0;
          for (Index k = 0; k < rdim; ++k)
            weight += reference[k] * Geometry::getReferenceVertexCoordinate(type, a, k);
          const Geometry::Vector& gradient = ip.getGradient(a);
          for (Index s = 0; s < sdim; ++s)
            result(s) += weight * gradient(s);
        }
        return result;
      }

      /// @brief Value and derivative of a one dimensional Lagrange polynomial.
      struct Line
      {
        Real value;
        Real derivative;
      };

      /**
       * @brief Evaluates a one dimensional Lagrange polynomial.
       * @param i Index of its segment node.
       * @param x Coordinate.
       * @return @f$ \ell_i(x) @f$ and @f$ \ell_i'(x) @f$.
       *
       * The product skips node @f$ i @f$ by selection, without a branch, and
       * the weight is selected the same way, since an index differing
       * between the threads of a warp would serialize its read.
       */
      KOKKOS_INLINE_FUNCTION
      Line getLine(Index i, Real x) const
      {
        Real value = 1;
        Real derivative = 0;
        Real weight = 0;
        for (Index m = 0; m <= K; ++m)
        {
          const Real e = x - m_line[m];
          const bool skip = m == i;
          derivative = skip ? derivative : derivative * e + value;
          value = skip ? value : value * e;
          weight = skip ? m_weights[m] : weight;
        }
        return {weight * value, weight * derivative};
      }

      /**
       * @brief Gets the powers @f$ r_k^p @f$, @f$ p \leq K @f$.
       * @param r Reference coordinates.
       * @return Table whose entry @f$ (k, p) @f$ is @f$ r_k^p @f$.
       */
      KOKKOS_INLINE_FUNCTION
      static Powers getPowers(const Geometry::Vector& r)
      {
        Powers powers;
        for (Index k = 0; k < r.size(); ++k)
        {
          powers[k][0] = 1;
          for (Index p = 1; p <= K; ++p)
            powers[k][p] = powers[k][p - 1] * r(k);
        }
        return powers;
      }

      /**
       * @brief Gets the exponents of the monomials of a reference geometry.
       * @param g Reference geometry.
       * @return Exponents @f$ (i, j, k) @f$ of the monomials of degree at most
       * K in each variable, of total degree at most K on a simplex and of
       * degree at most K in the first two variables together on a wedge.
       */
      static std::vector<std::array<std::uint8_t, 3>> getExponents(
        ::Rodin::Geometry::Polytope::Type g)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        const bool simplex = g == Type::Triangle || g == Type::Tetrahedron;
        const size_t dim = ::Rodin::Geometry::Polytope::Traits(g).getDimension();
        const size_t ki = dim > 0 ? K : 0;
        const size_t kj = dim > 1 ? K : 0;
        const size_t kk = dim > 2 ? K : 0;
        std::vector<std::array<std::uint8_t, 3>> result;
        for (size_t k = 0; k <= kk; ++k)
        {
          for (size_t j = 0; j <= kj; ++j)
          {
            for (size_t i = 0; i <= ki; ++i)
            {
              if ((simplex && i + j + k > K) || (g == Type::Wedge && i + j > K))
                continue;
              result.push_back({static_cast<std::uint8_t>(i),
                static_cast<std::uint8_t>(j), static_cast<std::uint8_t>(k)});
            }
          }
        }
        return result;
      }

      ::Kokkos::Array<Index, GeometryCount> m_count = {};
      ::Kokkos::Array<Index, GeometryCount> m_nodeOffsets = {};
      ::Kokkos::Array<Index, GeometryCount> m_coefficientOffsets = {};
      ::Kokkos::View<std::uint8_t*, MemorySpace> m_exponents;
      ::Kokkos::View<Real*, MemorySpace> m_nodes;
      ::Kokkos::View<Real*, MemorySpace> m_coefficients;
      ::Kokkos::Array<bool, GeometryCount> m_tensor = {};
      ::Kokkos::Array<Real, K + 1> m_line = {};
      ::Kokkos::Array<Real, K + 1> m_weights = {};
      ::Kokkos::View<std::uint8_t*, MemorySpace> m_lines;
  };

  /**
   * @brief Degrees of freedom of an H1 space of degree K on the device.
   *
   * Holds the local to global numbering of every polytope of the mesh,
   * copied from the host space, and the basis. Component @f$ r @f$ of the
   * scalar degree of freedom @f$ q @f$ has the global index @f$ q + r n @f$,
   * where @f$ n @f$ is the stride, the size of the scalar space.
   *
   * The quadratic simplices keep the closed form of QuadraticBasis, the
   * other polytopes take LagrangeBasis.
   *
   * @tparam K Polynomial degree.
   * @tparam MemorySpace Kokkos memory space of the tables.
   */
  template <size_t K, class MemorySpace>
  class LagrangeDOFs
  {
    public:
      /// @brief Empty tables.
      LagrangeDOFs() = default;

      /**
       * @brief Copies the numbering of a scalar host space to the device.
       * @param fes Scalar host space of degree K.
       * @param execution Execution space ordering the transfers.
       */
      template <class FES, class ExecutionSpace>
      LagrangeDOFs(const FES& fes, const ExecutionSpace& execution)
        : m_stride(fes.getSize()),
          m_basis(execution)
      {
        if constexpr (K == 2)
          m_quadratic = QuadraticBasis::get();
        const auto& mesh = fes.getMesh();
        assert(mesh.getDimension() < m_dofs.size());
        for (size_t d = 0; d <= mesh.getDimension(); ++d)
        {
          // Polytopes of one node count take a fixed stride, others offsets.
          const Index count = mesh.getPolytopeCount(d);
          Index uniform = count > 0 ? fes.getDOFs(d, 0).size() : 0;
          Index total = 0;
          for (Index i = 0; i < count; ++i)
          {
            const Index n = fes.getDOFs(d, i).size();
            total += n;
            if (n != uniform)
              uniform = 0;
          }
          m_counts[d] = uniform;
          m_dofs[d] = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, std::string("Rodin H1 degrees of freedom")),
            total);
          auto host = ::Kokkos::create_mirror_view(m_dofs[d]);
          Index position = 0;
          for (Index i = 0; i < count; ++i)
          {
            const auto& dofs = fes.getDOFs(d, i);
            for (Index q = 0; q < static_cast<Index>(dofs.size()); ++q)
              host(position++) = dofs(q);
          }
          ::Kokkos::deep_copy(execution, m_dofs[d], host);
          if (uniform == 0)
          {
            m_offsets[d] = ::Kokkos::View<Index*, MemorySpace>(
              ::Kokkos::view_alloc(execution, std::string("Rodin H1 offsets")),
              count + 1);
            auto offsets = ::Kokkos::create_mirror_view(m_offsets[d]);
            offsets(0) = 0;
            for (Index i = 0; i < count; ++i)
              offsets(i + 1) = offsets(i) + fes.getDOFs(d, i).size();
            ::Kokkos::deep_copy(execution, m_offsets[d], offsets);
          }
        }
        execution.fence("Rodin H1 degrees of freedom transfer");
      }

      /**
       * @brief Gets the number of local degrees of freedom of a polytope.
       * @param vertices Number of vertices of a polytope.
       * @return Maximum node count over the geometries of at most
       * @p vertices vertices.
       */
      KOKKOS_INLINE_FUNCTION
      static constexpr Index getCapacity(Index vertices)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        using Basis = LagrangeBasis<K, MemorySpace>;
        Index result = 1;
        if (vertices >= 2)
          result = Basis::getCount(Type::Segment);
        if (vertices >= 3)
          result = Basis::getCount(Type::Triangle);
        if (vertices >= 4)
        {
          result =
            Basis::getCount(Type::Quadrilateral) > Basis::getCount(Type::Tetrahedron)
            ? Basis::getCount(Type::Quadrilateral)
            : Basis::getCount(Type::Tetrahedron);
        }
        if (vertices >= 6 && Basis::getCount(Type::Wedge) > result)
          result = Basis::getCount(Type::Wedge);
        if (vertices >= 8 && Basis::getCount(Type::Hexahedron) > result)
          result = Basis::getCount(Type::Hexahedron);
        return result;
      }

      /// @brief Gets the offset between components. @return Scalar space size.
      KOKKOS_INLINE_FUNCTION
      Index getStride() const
      {
        return m_stride;
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param polytope Polytope of the mesh.
       * @return Scalar global indices in local order.
       */
      KOKKOS_INLINE_FUNCTION
      Geometry::IndexView bind(const Geometry::PolytopeView& polytope) const
      {
        const Index d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const Index n = m_counts[d];
        if (n > 0)
          return {m_dofs[d].data() + i * n, n};
        const Index begin = m_offsets[d](i);
        return {m_dofs[d].data() + begin, m_offsets[d](i + 1) - begin};
      }

      /// @param ip Integration point. @param q Local index. @return Basis value.
      template <class Point>
      KOKKOS_INLINE_FUNCTION Real getBasis(const Point& ip, Index q) const
      {
        if constexpr (K == 2)
        {
          if (isSimplex(ip.getPoint().getPolytope()))
            return m_quadratic.getBasis(ip, q);
        }
        return m_basis.getBasis(ip, q);
      }

      /// @param ip Integration point. @param q Local index. @return Physical gradient.
      template <class Point>
      KOKKOS_INLINE_FUNCTION Geometry::Vector getGradient(const Point& ip, Index q) const
      {
        if constexpr (K == 2)
        {
          if (isSimplex(ip.getPoint().getPolytope()))
            return m_quadratic.getGradient(ip, q);
        }
        return m_basis.getGradient(ip, q);
      }

      /// @param polytope Polytope. @param q Local index. @return Reference node.
      KOKKOS_INLINE_FUNCTION
      Geometry::Vector getNode(const Geometry::PolytopeView& polytope, Index q) const
      {
        if constexpr (K == 2)
        {
          if (isSimplex(polytope))
            return m_quadratic.getNode(polytope, q);
        }
        return m_basis.getNode(polytope, q);
      }

    private:
      /**
       * @brief Tells whether a polytope is a simplex.
       * @param polytope Device polytope.
       * @return Whether @p polytope is a point, a segment, a triangle or a
       * tetrahedron.
       */
      KOKKOS_INLINE_FUNCTION
      static bool isSimplex(const Geometry::PolytopeView& polytope)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        const Type g = polytope.getGeometry();
        return g == Type::Point || g == Type::Segment || g == Type::Triangle ||
          g == Type::Tetrahedron;
      }

      Index m_stride = 0;
      QuadraticBasis m_quadratic;
      LagrangeBasis<K, MemorySpace> m_basis;
      ::Kokkos::Array<Index, 4> m_counts = {};
      ::Kokkos::Array<::Kokkos::View<Index*, MemorySpace>, 4> m_dofs;
      ::Kokkos::Array<::Kokkos::View<Index*, MemorySpace>, 4> m_offsets;
  };
}

namespace Rodin::Kokkos::Variational::Internal
{
  /**
   * @brief Host copy of a resident mesh with the H1 spaces of degree K built
   * on it.
   *
   * The host spaces number the degrees of freedom and answer the form
   * language queries; the scalar numbering is copied to the device.
   *
   * @tparam K Polynomial degree.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Backend>
  class LagrangeMirror
  {
    public:
      /// @brief Resident mesh type.
      using MeshType = ::Rodin::Geometry::Mesh<
        ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>;

      /// @brief Host mesh type.
      using HostMeshType = ::Rodin::Geometry::Mesh<::Rodin::Context::Local>;

      /// @brief Memory space of the device tables.
      using MemorySpace = typename MeshType::MemorySpace;

      /// @brief Scalar host space.
      using ScalarSpace = ::Rodin::Variational::H1<K, Real, HostMeshType>;

      /// @brief Vector host space.
      using VectorSpace =
        ::Rodin::Variational::H1<K, ::Rodin::Math::SpatialVector<Real>, HostMeshType>;

      /**
       * @brief Builds the host spaces and copies the numbering to the device.
       * @param mesh Resident mesh without pyramids. The incidences
       * @f$ d \to d - 1 @f$ it lacks are computed on the host copy.
       * @param vdim Vector dimension, zero for a scalar space only.
       */
      LagrangeMirror(const MeshType& mesh, size_t vdim)
        : m_mesh(supported(mesh)),
          m_scalar(std::integral_constant<size_t, K>{}, m_mesh),
          m_dofs(m_scalar, mesh.getContext().getExecutionSpace())
      {
        if (vdim > 0)
          m_vector.emplace(std::integral_constant<size_t, K>{}, m_mesh, vdim);
      }

      /**
       * @brief Not copyable, the host spaces referring to the owned host mesh.
       * @param other Mirror not copied.
       */
      LagrangeMirror(const LagrangeMirror& other) = delete;

      /**
       * @brief Not copy assignable, the host spaces referring to the owned host mesh.
       * @param other Mirror not copied.
       * @return Reference to this mirror.
       */
      LagrangeMirror& operator=(const LagrangeMirror& other) = delete;

      /// @brief Gets the scalar host space. @return Scalar space.
      const ScalarSpace& getScalar() const
      {
        return m_scalar;
      }

      /// @brief Gets the vector host space. @return Vector space.
      const VectorSpace& getVector() const
      {
        assert(m_vector);
        return *m_vector;
      }

      /// @brief Gets the device numbering. @return Scalar degrees of freedom.
      const LagrangeDOFs<K, MemorySpace>& getDOFs() const
      {
        return m_dofs;
      }

    private:
      /**
       * @brief Copies a mesh to the host after checking its geometries.
       * @param mesh Resident mesh.
       * @return Host copy of @p mesh, with the incidences @f$ d \to d - 1 @f$.
       */
      static HostMeshType supported(const MeshType& mesh)
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        if (mesh.getPolytopeCount(Type::Pyramid) > 0)
        {
          Alert::Exception() << "Kokkos H1 spaces do not support pyramids."
                             << Alert::Raise;
        }
        // The incidences missing on the device only number the host
        // closure, the kernels never bind the polytopes they create.
        HostMeshType host = mesh.toHost();
        for (size_t d = host.getDimension(); d > 1; --d)
        {
          if (host.getConnectivity().getIncidence(d, d - 1).size() == 0)
            host.getConnectivity().compute(d, d - 1);
        }
        return host;
      }

      HostMeshType m_mesh;
      ScalarSpace m_scalar;
      std::optional<VectorSpace> m_vector;
      LagrangeDOFs<K, MemorySpace> m_dofs;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup H1Specializations
   * @brief Real valued H1 Lagrange space on a resident Kokkos mesh.
   *
   * The space numbers its degrees of freedom as the host space of the same
   * degree built on Geometry::Mesh::toHost(), so the two share their global
   * numbering. The host space is kept for the form language queries, and
   * its numbering is copied to the device, where the basis is evaluated by
   * Rodin::Kokkos::Variational::LagrangeDOFs.
   *
   * Every geometry but the pyramid is available on the device.
   *
   * @tparam K Polynomial degree.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Backend>
  class H1<K, Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        H1<K, Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = ScalarType;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = H1Element<K, RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType =
        Rodin::Kokkos::Variational::LagrangeDOFs<K, typename MeshType::MemorySpace>;

      /**
       * @brief Constructs the space on a resident mesh.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Mesh without pyramids supporting the space, which must
       * outlive it.
       */
      H1([[maybe_unused]] std::integral_constant<size_t, K> degree, const MeshType& mesh)
        : m_mesh(mesh),
          m_mirror(std::make_shared<const MirrorType>(mesh, 0))
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_mirror(other.m_mirror)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_mirror(std::move(other.m_mirror))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_mirror = other.m_mirror;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_mirror = std::move(other.m_mirror);
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_mirror->getScalar().getFiniteElement(d, i);
      }

      /// @brief Gets the number of degrees of freedom. @return Size of the space.
      size_t getSize() const override
      {
        return m_mirror->getScalar().getSize();
      }

      /// @brief Gets the vector dimension. @return One, the space being scalar.
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Global indices in local order.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_mirror->getScalar().getDOFs(d, i);
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local index within the polytope.
       * @return Global index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        return m_mirror->getScalar().getGlobalIndex(idx, local);
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Scalar H1 space on @p mesh.
       */
      H1<K, RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return H1<K, RangeType, Geometry::Mesh<Context::Local>>(
          std::integral_constant<size_t, K>{}, mesh);
      }

      /**
       * @brief Gets the pullback of a callable on a polytope, as on the host.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] v Callable of the physical points.
       * @return Callable of the reference coordinates.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        using Pullback = typename H1<K, RangeType,
          Geometry::Mesh<Context::Local>>::template Pullback<Callable>;
        const auto& [d, i] = idx;
        return Pullback(*getMesh().getPolytope(d, i), std::forward<Callable>(v));
      }

      /// @brief Gets the device degrees of freedom. @return Device numbering and basis.
      const DOFsType& getDOFsView() const
      {
        return m_mirror->getDOFs();
      }

    private:
      using MirrorType = Rodin::Kokkos::Variational::Internal::LagrangeMirror<K, Backend>;

      std::reference_wrapper<const MeshType> m_mesh;
      std::shared_ptr<const MirrorType> m_mirror;
  };

  /**
   * @ingroup H1Specializations
   * @brief Vector valued H1 Lagrange space on a resident Kokkos mesh.
   *
   * The space numbers its degrees of freedom as the host vector space:
   * component @f$ r @f$ of the scalar degree of freedom @f$ q @f$ has the
   * global index @f$ q + r n @f$, where @f$ n @f$ is the size of the scalar
   * space, and the local index @f$ a\,d + r @f$ when @f$ q @f$ is the local
   * node @f$ a @f$ of a polytope and @f$ d @f$ the vector dimension.
   *
   * @tparam K Polynomial degree.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Backend>
  class H1<K, ::Rodin::Math::SpatialVector<Real>,
    Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>
    final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        H1<K, ::Rodin::Math::SpatialVector<Real>,
          Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = ::Rodin::Math::SpatialVector<ScalarType>;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = H1Element<K, RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType =
        Rodin::Kokkos::Variational::LagrangeDOFs<K, typename MeshType::MemorySpace>;

      /**
       * @brief Constructs the vector space on a resident mesh.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Mesh without pyramids supporting the space, which must
       * outlive it.
       * @param[in] vdim Vector dimension, at most three.
       */
      H1([[maybe_unused]] std::integral_constant<size_t, K> degree, const MeshType& mesh,
        size_t vdim)
        : m_mesh(mesh),
          m_vdim(vdim),
          m_mirror(std::make_shared<const MirrorType>(mesh, vdim))
      {
        assert(vdim > 0 && vdim <= RODIN_MAXIMAL_SPACE_DIMENSION);
      }

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_mirror(other.m_mirror)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(other.m_vdim),
          m_mirror(std::move(other.m_mirror))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
          m_mirror = other.m_mirror;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = other.m_vdim;
          m_mirror = std::move(other.m_mirror);
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Vector element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_mirror->getVector().getFiniteElement(d, i);
      }

      /// @brief Gets the number of degrees of freedom. @return Size of the space.
      size_t getSize() const override
      {
        return m_mirror->getVector().getSize();
      }

      /// @brief Gets the vector dimension. @return Number of components.
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Global indices in local order.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_mirror->getVector().getDOFs(d, i);
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local index within the polytope.
       * @return Global index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        return m_mirror->getVector().getGlobalIndex(idx, local);
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Vector H1 space on @p mesh.
       */
      H1<K, RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return H1<K, RangeType, Geometry::Mesh<Context::Local>>(
          std::integral_constant<size_t, K>{}, mesh, m_vdim);
      }

      /**
       * @brief Gets the pullback of a callable on a polytope, as on the host.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] v Callable of the physical points.
       * @return Callable of the reference coordinates.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        using Pullback = typename H1<K, RangeType,
          Geometry::Mesh<Context::Local>>::template Pullback<Callable>;
        const auto& [d, i] = idx;
        return Pullback(*getMesh().getPolytope(d, i), std::forward<Callable>(v));
      }

      /// @brief Gets the device degrees of freedom. @return Scalar device numbering and basis.
      const DOFsType& getDOFsView() const
      {
        return m_mirror->getDOFs();
      }

    private:
      using MirrorType = Rodin::Kokkos::Variational::Internal::LagrangeMirror<K, Backend>;

      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
      std::shared_ptr<const MirrorType> m_mirror;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief Quadratic H1 space on a resident mesh.
  template <class Range, class Backend>
  using KokkosH1Space = ::Rodin::Variational::H1<2, Range,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;

  /// @brief Real valued quadratic H1 space on a resident mesh.
  template <class Backend>
  using KokkosH1 = KokkosH1Space<Real, Backend>;

  /// @brief Vector valued quadratic H1 space on a resident mesh.
  template <class Backend>
  using KokkosVectorH1 = KokkosH1Space<::Rodin::Math::SpatialVector<Real>, Backend>;
}

#endif
