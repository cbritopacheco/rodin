/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_HPP
#define RODIN_VARIATIONAL_H1_H1ELEMENT_HPP

/**
 * @file
 * @brief Template implementation for H1Element class.
 *
 * This file contains the template member function implementations for
 * H1Element, including:
 * - Basis function evaluation using modal Dubiner basis
 * - First-order derivative computation via chain rule
 * - Vandermonde matrix inversion for nodal-to-modal conversion
 *
 * The implementation uses:
 * - **Segments**: Classical Lagrange interpolation with GLL nodes
 * - **Triangles**: Dubiner modal basis with Fekete nodes
 * - **Quadrilaterals**: Tensor product of 1D GLL Lagrange bases
 * - **Tetrahedra**: 3D Dubiner modal basis with Fekete nodes
 * - **Wedges**: Product of triangle Fekete basis and 1D GLL basis
 */

#include "Rodin/Math/Common.h"

#include "Rodin/QF/QuadratureFormula.h"

#include "H1Element.h"

#include "LagrangeBasis.h"
#include "Dubiner.h"
#include "Fekete.h"
#include "WarpBlend.h"
#include "GLL.h"
#include "LegendrePolynomial.h"

/// @cond RODIN_DOXYGEN_INTERNAL
namespace Rodin::Variational
{
  template <size_t K>
  struct PyramidIndex
  {
    static constexpr size_t Count = (K + 1) * (K + 2) * (2 * K + 3) / 6;

    struct IJ
    {
      size_t i;
      size_t j;
    };

    static constexpr size_t getLayerOffset(size_t layer)
    {
      size_t out = 0;
      for (size_t k = 0; k < layer; ++k)
      {
        const size_t n = K - k + 1;
        out += n * n;
      }
      return out;
    }

    static constexpr size_t getIndex(size_t i, size_t j, size_t k)
    {
      const size_t n = K - k + 1;
      return getLayerOffset(k) + j * n + i;
    }

    static constexpr void decode(size_t idx, size_t& i, size_t& j, size_t& k)
    {
      size_t rem = idx;
      for (k = 0; k <= K; ++k)
      {
        const size_t n = K - k + 1;
        const size_t layer = n * n;
        if (rem < layer)
        {
          j = rem / n;
          i = rem % n;
          return;
        }
        rem -= layer;
      }

      i = j = k = 0;
    }

    static constexpr IJ getTriangleLattice(size_t alpha)
    {
      size_t pos = 0;
      for (size_t j = 0; j <= K; ++j)
      {
        for (size_t i = 0; i <= K - j; ++i, ++pos)
        {
          if (pos == alpha)
            return IJ{ i, j };
        }
      }
      return IJ{ 0, 0 };
    }

    static constexpr size_t getSideIndex(size_t local, size_t alpha)
    {
      const auto ij = getTriangleLattice(alpha);
      const size_t i = ij.i;
      const size_t k = ij.j;
      const size_t n = K - k;

      switch (local)
      {
        case 1:
          return getIndex(i, 0, k);
        case 2:
          return getIndex(n, i, k);
        case 3:
          return getIndex(n - i, n, k);
        case 4:
          return getIndex(0, n - i, k);
        default:
          return 0;
      }
    }
  };

  template <size_t K>
  class BernsteinPyramid
  {
    public:
      static Real getBinomial(size_t n, size_t k)
      {
        if (k > n)
          return 0;
        if (k > n - k)
          k = n - k;

        Real out = 1;
        for (size_t i = 1; i <= k; ++i)
        {
          out *= static_cast<Real>(n + 1 - i);
          out /= static_cast<Real>(i);
        }
        return out;
      }

      static Real getBasis(size_t n, size_t i, Real x)
      {
        if (i > n)
          return 0;

        Real out = getBinomial(n, i);
        for (size_t p = 0; p < i; ++p)
          out *= x;
        for (size_t p = 0; p < n - i; ++p)
          out *= (Real(1) - x);
        return out;
      }

      static Real getDerivative(size_t n, size_t i, Real x)
      {
        if (n == 0)
          return 0;

        Real left = 0;
        if (i > 0)
          left = getBasis(n - 1, i - 1, x);

        Real right = 0;
        if (i < n)
          right = getBasis(n - 1, i, x);

        return static_cast<Real>(n) * (left - right);
      }
  };

  template <size_t K>
  class PyramidModal
  {
    public:
      static Real getBasis(size_t mode, const Math::SpatialPoint& r)
      {
        size_t i, j, k;
        PyramidIndex<K>::decode(mode, i, j, k);

        const Real z = r.z();
        const Real q = Real(1) - z;
        if (q <= RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE)
          return k == K ? Real(1) : Real(0);

        const size_t n = K - k;
        const Real a = r.x() / q;
        const Real b = r.y() / q;

        return BernsteinPyramid<K>::getBasis(n, i, a)
             * BernsteinPyramid<K>::getBasis(n, j, b)
             * BernsteinPyramid<K>::getBasis(K, k, z);
      }

      static Real getDerivative(size_t mode, size_t deriv, const Math::SpatialPoint& r)
      {
        size_t i, j, k;
        PyramidIndex<K>::decode(mode, i, j, k);

        const Real z = r.z();
        const Real q = Real(1) - z;
        if (q <= RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE)
          return 0;

        const size_t n = K - k;
        const Real a = r.x() / q;
        const Real b = r.y() / q;

        const Real Ba  = BernsteinPyramid<K>::getBasis(n, i, a);
        const Real Bb  = BernsteinPyramid<K>::getBasis(n, j, b);
        const Real Bz  = BernsteinPyramid<K>::getBasis(K, k, z);
        const Real dBa = BernsteinPyramid<K>::getDerivative(n, i, a);
        const Real dBb = BernsteinPyramid<K>::getDerivative(n, j, b);
        const Real dBz = BernsteinPyramid<K>::getDerivative(K, k, z);

        if (deriv == 0)
          return dBa * Bb * Bz / q;
        if (deriv == 1)
          return Ba * dBb * Bz / q;
        if (deriv == 2)
          return (dBa * (a / q) * Bb + Ba * dBb * (b / q)) * Bz
               + Ba * Bb * dBz;
        return 0;
      }
  };

  template <size_t K>
  class VandermondePyramid
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static const Math::Matrix<Real> s_vandermonde = [] {
          constexpr size_t N = PyramidIndex<K>::Count;
          const auto& nodes =
            H1Element<K, Real>::getNodes(Geometry::Polytope::Type::Pyramid);
          Math::Matrix<Real> vandermonde(N, N);

          for (size_t node = 0; node < N; ++node)
          {
            for (size_t mode = 0; mode < N; ++mode)
              vandermonde(node, mode) = PyramidModal<K>::getBasis(mode, nodes[node]);
          }
          return vandermonde;
        }();

        return s_vandermonde;
      }

      static const Math::Matrix<Real>& getInverse()
      {
        static const Math::Matrix<Real> s_inv = [] {
          const auto& V = getMatrix();
          Eigen::BDCSVD<Math::Matrix<Real>> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(V.rows(), V.cols());
          return svd.solve(I).eval();
        }();

        return s_inv;
      }
  };

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::Tabulation&
  H1Element<K, Scalar>::getTabulation(const QF::QuadratureFormulaBase& qf) const
  {
    struct CacheEntry
    {
      const QF::QuadratureFormulaBase* qf;
      Geometry::Polytope::Type g;
      size_t nqp;
      bool valid;
      Tabulation tab;

      CacheEntry()
        : qf(nullptr),
          g(Geometry::Polytope::Type::Point),
          nqp(0),
          valid(false),
          tab()
      {}
    };

    struct Cache
    {
      std::array<CacheEntry, 8> e;
      size_t next = 0; // eviction pointer
    };

    static thread_local Cache s_cache;

    const auto g   = this->getGeometry();
    const auto nqp = qf.getSize();

    // 1) lookup
    for (auto& ce : s_cache.e)
    {
      if (ce.valid && ce.qf == &qf && ce.g == g && ce.nqp == nqp)
        return ce.tab;
    }

    // 2) miss -> rebuild into an entry
    CacheEntry& ce = s_cache.e[s_cache.next];
    s_cache.next = (s_cache.next + 1) % s_cache.e.size();

    ce.valid = true;
    ce.qf = &qf;
    ce.g = g;
    ce.nqp = nqp;

    Tabulation& t = ce.tab;

    const size_t dim  = Geometry::Polytope::Traits(g).getDimension();
    const size_t ndof = this->getCount();

    t.nqp  = nqp;
    t.ndof = ndof;
    t.dim  = dim;

    t.phi.resize(nqp * ndof);
    t.dphi.resize(nqp * ndof * dim);

    // Cache basis handles (cheap, avoids repeated getBasis(a) calls)
    std::vector<std::reference_wrapper<const BasisFunction>> bf;
    bf.reserve(ndof);
    for (size_t a = 0; a < ndof; ++a)
      bf.emplace_back(this->getBasis(a));

    for (size_t qp = 0; qp < nqp; ++qp)
    {
      const auto& rc = qf.getPoint(qp);

      Scalar* phi_row  = t.phi.data()  + qp * ndof;
      Scalar* dphiRow = t.dphi.data() + qp * ndof * dim;

      for (size_t a = 0; a < ndof; ++a)
      {
        const auto& bfa = bf[a].get();

        phi_row[a] = bfa(rc);

        Scalar* d = dphiRow + a * dim;
        for (size_t i = 0; i < dim; ++i)
          d[i] = bfa.template getDerivative<1>(i)(rc);
      }
    }

    return t;
  }

  template <size_t K, class Scalar>
  constexpr
  const Math::SpatialPoint& H1Element<K, Scalar>::getNode(size_t i) const
  {
    return this->getNodes(this->getGeometry())[i];
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::LinearForm& H1Element<K, Scalar>::getLinearForm(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // These objects depend only on the compile-time order and reference
    // geometry. Sharing their immutable storage avoids per-thread replicas.
    const auto makeLinearForms = [this, g] {
      std::vector<LinearForm> forms;
      const size_t count = getCount();
      forms.reserve(count);
      for (size_t j = 0; j < count; ++j)
        forms.emplace_back(j, g);
      return forms;
    };

    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        static const std::vector<LinearForm> s_lfs = makeLinearForms();
        return s_lfs[i];
      }
    }

    // Fallback (should never happen)
    static const LinearForm s_null(0, Geometry::Polytope::Type::Point);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  const typename H1Element<K, Scalar>::BasisFunction& H1Element<K, Scalar>::getBasis(size_t i) const
  {
    const Geometry::Polytope::Type g = this->getGeometry();

    // Basis descriptors are immutable after construction and are safe to share
    // across concurrent element evaluations.
    const auto makeBasis = [this, g] {
      std::vector<BasisFunction> basis;
      const size_t count = getCount();
      basis.reserve(count);
      for (size_t j = 0; j < count; ++j)
        basis.emplace_back(j, g);
      return basis;
    };

    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Segment:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Triangle:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Wedge:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        static const std::vector<BasisFunction> s_bs = makeBasis();
        return s_bs[i];
      }
    }

    // Fallback (should never happen)
    static const BasisFunction s_null(0, Geometry::Polytope::Type::Point);
    assert(false);
    return s_null;
  }

  template <size_t K, class Scalar>
  Scalar H1Element<K, Scalar>::BasisFunction::operator()(
      const Math::SpatialPoint& r) const
  {
    switch (m_g)
    {
      case Geometry::Polytope::Type::Point:
      {
        return Scalar(1);
      }
      case Geometry::Polytope::Type::Segment:
      {
        // Canonical [0,1] segment with GLL01 nodes
        return LagrangeBasisSegment<K>::getBasis(m_local, r.value());
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // Modal Dubiner basis + Vandermonde on Fekete nodes
        const auto& inverse = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar result = Scalar(0);
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Real psi;
            DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
            result += inverse(modeIdx, m_local) * psi;
            ++modeIdx;
          });
        });

        return result;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // Tensor product GLL01 × GLL01 (Lagrange)
        const size_t jIdx = m_local / (K + 1);
        const size_t iIdx = m_local % (K + 1);
        return LagrangeBasisQuadrilateral<K>::getBasis(iIdx, jIdx, r.x(), r.y());
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // Modal Dubiner basis + Vandermonde on tetra Fekete nodes
        const auto& inverse = VandermondeTetrahedron<K>::getInverse();

        Real ac, bc, cc;
        DubinerTetrahedron<K>::getCollapsed(
            ac, bc, cc, r.x(), r.y(), r.z());

        Real result = 0;
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P - Q>([&](auto rIdx) {
              constexpr size_t R = rIdx.value;
              Real psi;
              DubinerTetrahedron<K>::template getBasis<P, Q, R>(psi, ac, bc, cc);
              result += inverse(modeIdx, m_local) * psi;
              ++modeIdx;
            });
          });
        });

        return result;
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        const auto& inverse = VandermondePyramid<K>::getInverse();
        constexpr size_t count = PyramidIndex<K>::Count;

        Scalar result = Scalar(0);
        for (size_t mode = 0; mode < count; ++mode)
          result += inverse(mode, m_local) * PyramidModal<K>::getBasis(mode, r);

        return result;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        // Number of triangle DOFs (Fekete triangle)
        constexpr size_t ntri = FeketeTriangle<K>::Count; // same as FeketeTriangle<K>::Count

        const size_t k     = m_local / ntri; // z-index, 0..K
        const size_t alpha = m_local % ntri; // triangle node index

        // --- triangle factor: nodal basis on Fekete nodes ---
        const auto& Vinv = VandermondeTriangle<K>::getInverse();

        Real rc, sc;
        DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

        Scalar triVal = Scalar(0);
        size_t modeIdx = 0;

        Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
          constexpr size_t P = pIdx.value;
          Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
            constexpr size_t Q = qIdx.value;
            Real psi;
            DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
            triVal += Vinv(modeIdx, alpha) * psi;
            ++modeIdx;
          });
        });

        // --- segment factor in z (GLL01, Lagrange) ---
        const Real z = r(2);
        const Real segVal = LagrangeBasisSegment<K>::getBasis(k, z);

        return triVal * segVal;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        // Tensor product GLL01 × GLL01 × GLL01
        // Node ordering matches getNodes(): i + (K+1)*(j + (K+1)*k)
        const size_t n1 = K + 1;
        const size_t k  = m_local / (n1 * n1);
        const size_t r2 = m_local % (n1 * n1);
        const size_t j  = r2 / n1;
        const size_t i  = r2 % n1;

        const Real x = r.x();
        const Real y = r.y();
        const Real z = r.z();

        const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
        const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
        const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

        return static_cast<Scalar>(lx * ly * lz);
      }
    }

    return Math::nan<Scalar>();
  }

  template <size_t K, class Scalar>
  template <size_t Order>
  Scalar H1Element<K, Scalar>::BasisFunction::DerivativeFunction<Order>::operator()(
      const Math::SpatialPoint& r) const
  {
    if constexpr (Order == 0)
    {
      return BasisFunction(m_local, m_g)(r);
    }
    else if constexpr (Order == 1)
    {
      switch (m_g)
      {
        case Geometry::Polytope::Type::Point:
        {
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Segment:
        {
          // \partialφ_i / \partialx on [0,1] with GLL01 nodes
          assert(m_i == 0);
          return LagrangeBasisSegment<K>::getDerivative(
              m_local, r.value());
        }

        case Geometry::Polytope::Type::Triangle:
        {
          // Dubiner modal gradients + chain rule (r,s) → (x,y)
          const auto& Vinv = VandermondeTriangle<K>::getInverse();

          Scalar rc, sc;
          DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

          const Scalar x = r.x();
          const Scalar y = r.y();
          const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

          Scalar result = Scalar(0);
          size_t modeIdx = 0;

          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;

              Scalar dpsi_dr = Scalar(0), dpsi_ds = Scalar(0);
              DubinerTriangle<K>::template getGradient<P, Q>(dpsi_dr, dpsi_ds, rc, sc);

              Scalar dpsi_dx = Scalar(0), dpsi_dy = Scalar(0);

              // r = 2x/(1-y) - 1, s = 2y - 1
              if (Math::abs(Scalar(1) - y) > eps)
              {
                const Scalar denom = Scalar(1) - y;
                const Scalar dr_dx = Scalar(2) / denom;
                const Scalar dr_dy = Scalar(2) * x / (denom * denom);
                const Scalar ds_dx = Scalar(0);
                const Scalar ds_dy = Scalar(2);

                dpsi_dx = dpsi_dr * dr_dx + dpsi_ds * ds_dx;
                dpsi_dy = dpsi_dr * dr_dy + dpsi_ds * ds_dy;
              }

              if (m_i == 0) // \partial/\partialx
                result += Vinv(modeIdx, m_local) * dpsi_dx;
              else if (m_i == 1) // \partial/\partialy
                result += Vinv(modeIdx, m_local) * dpsi_dy;

              ++modeIdx;
            });
          });

          return result;
        }

        case Geometry::Polytope::Type::Quadrilateral:
        {
          // Tensor product Lagrange on GLL01 × GLL01
          const size_t jIdx = m_local / (K + 1);
          const size_t iIdx = m_local % (K + 1);

          if (m_i == 0) // \partial/\partialx
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
              iIdx, jIdx, 0, r.x(), r.y());
          }
          else if (m_i == 1) // \partial/\partialy
          {
            return LagrangeBasisQuadrilateral<K>::getDerivative(
              iIdx, jIdx, 1, r.x(), r.y());
          }
          return Scalar(0);
        }

        case Geometry::Polytope::Type::Tetrahedron:
        {
          // Dubiner modal gradients + chain rule (a,b,c) → (x,y,z)
          const auto& Vinv = VandermondeTetrahedron<K>::getInverse();

          Scalar ac, bc, cc;
          DubinerTetrahedron<K>::getCollapsed(
              ac, bc, cc, r.x(), r.y(), r.z());

          const Scalar x = r.x();
          const Scalar y = r.y();
          const Scalar z = r.z();
          const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

          Scalar result = Scalar(0);
          size_t modeIdx = 0;

          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P - Q>([&](auto rIdx) {
                constexpr size_t R = rIdx.value;

                Scalar dpsi_da = Scalar(0);
                Scalar dpsi_db = Scalar(0);
                Scalar dpsi_dc = Scalar(0);
                DubinerTetrahedron<K>::template getGradient<P, Q, R>(
                  dpsi_da, dpsi_db, dpsi_dc, ac, bc, cc);

                Scalar dpsi_dx = Scalar(0);
                Scalar dpsi_dy = Scalar(0);
                Scalar dpsi_dz = Scalar(0);

                const Scalar denom2 = Scalar(1) - z; // 1 - z
                const Scalar denom3 = Scalar(1) - y - z; // 1 - y - z

                if (Math::abs(denom2) > eps && Math::abs(denom3) > eps)
                {
                  // a = 2x / (1 - y - z) - 1
                  const Scalar da_dx = Scalar(2) / denom3;
                  const Scalar da_dy = Scalar(2) * x / (denom3 * denom3);
                  const Scalar da_dz = da_dy;

                  // b = 2y / (1 - z) - 1
                  const Scalar db_dx = Scalar(0);
                  const Scalar db_dy = Scalar(2) / denom2;
                  const Scalar db_dz = Scalar(2) * y / (denom2 * denom2);

                  // c = 2z - 1
                  const Scalar dc_dx = Scalar(0);
                  const Scalar dc_dy = Scalar(0);
                  const Scalar dc_dz = Scalar(2);

                  dpsi_dx = dpsi_da * da_dx + dpsi_db * db_dx + dpsi_dc * dc_dx;

                  dpsi_dy = dpsi_da * da_dy + dpsi_db * db_dy + dpsi_dc * dc_dy;

                  dpsi_dz = dpsi_da * da_dz + dpsi_db * db_dz + dpsi_dc * dc_dz;
                }

                if (m_i == 0) // \partial/\partialx
                  result += Vinv(modeIdx, m_local) * dpsi_dx;
                else if (m_i == 1) // \partial/\partialy
                  result += Vinv(modeIdx, m_local) * dpsi_dy;
                else if (m_i == 2) // \partial/\partialz
                  result += Vinv(modeIdx, m_local) * dpsi_dz;

                ++modeIdx;
              });
            });
          });

          return result;
        }
        case Geometry::Polytope::Type::Pyramid:
        {
          const auto& inverse = VandermondePyramid<K>::getInverse();
          constexpr size_t count = PyramidIndex<K>::Count;

          Scalar result = Scalar(0);
          for (size_t mode = 0; mode < count; ++mode)
          {
            result += inverse(mode, m_local)
                    * PyramidModal<K>::getDerivative(mode, m_i, r);
          }

          return result;
        }

        case Geometry::Polytope::Type::Wedge:
        {
          constexpr size_t ntri = FeketeTriangle<K>::Count;

          const size_t k = m_local / ntri;
          const size_t alpha = m_local % ntri;

          const Real z = r(2);

          if (m_i < 2) // \partial/\partialx or \partial/\partialy
          {
            // --- triangle gradient (same as Triangle case, but index = alpha) ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();

            Scalar rc, sc;
            DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

            const Scalar x   = r.x();
            const Scalar y   = r.y();
            const Scalar eps = RODIN_VARIATIONAL_H1ELEMENT_TOLERANCE;

            Scalar triDeriv = Scalar(0);
            size_t modeIdx = 0;

            Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
              constexpr size_t P = pIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
                constexpr size_t Q = qIdx.value;

                Scalar dpsi_dr = Scalar(0), dpsi_ds = Scalar(0);
                DubinerTriangle<K>::template getGradient<P, Q>(dpsi_dr, dpsi_ds, rc, sc);

                Scalar dpsi_dx = Scalar(0), dpsi_dy = Scalar(0);

                // r = 2x/(1-y) - 1, s = 2y - 1
                if (Math::abs(Scalar(1) - y) > eps)
                {
                  const Scalar denom = Scalar(1) - y;
                  const Scalar dr_dx = Scalar(2) / denom;
                  const Scalar dr_dy = Scalar(2) * x / (denom * denom);
                  const Scalar ds_dx = Scalar(0);
                  const Scalar ds_dy = Scalar(2);

                  dpsi_dx = dpsi_dr * dr_dx + dpsi_ds * ds_dx;
                  dpsi_dy = dpsi_dr * dr_dy + dpsi_ds * ds_dy;
                }

                if (m_i == 0) // \partial/\partialx
                  triDeriv += Vinv(modeIdx, alpha) * dpsi_dx;
                else if (m_i == 1) // \partial/\partialy
                  triDeriv += Vinv(modeIdx, alpha) * dpsi_dy;

                ++modeIdx;
              });
            });

            // --- segment value in z ---
            const Real segVal = LagrangeBasisSegment<K>::getBasis(k, z);
            return triDeriv * segVal;
          }
          else // m_i == 2 → \partial/\partialz
          {
            // --- triangle value (same as in BasisFunction wedge case) ---
            const auto& Vinv = VandermondeTriangle<K>::getInverse();

            Real rc, sc;
            DubinerTriangle<K>::getCollapsed(rc, sc, r.x(), r.y());

            Scalar triVal = Scalar(0);
            size_t modeIdx = 0;

            Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
              constexpr size_t P = pIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
                constexpr size_t Q = qIdx.value;
                Real psi;
                DubinerTriangle<K>::template getBasis<P, Q>(psi, rc, sc);
                triVal += Vinv(modeIdx, alpha) * psi;
                ++modeIdx;
              });
            });

            // --- 1D derivative in z ---
            const Real dseg = LagrangeBasisSegment<K>::getDerivative(k, z);
            return triVal * dseg;
          }
        }
        case Geometry::Polytope::Type::Hexahedron:
        {
          // Tensor product derivative on [0,1]^3 with GLL nodes.
          // Node ordering: i + (K+1)*(j + (K+1)*k)
          const size_t n1 = K + 1;
          const size_t k  = m_local / (n1 * n1);
          const size_t r2 = m_local % (n1 * n1);
          const size_t j  = r2 / n1;
          const size_t i  = r2 % n1;

          const Real x = r.x();
          const Real y = r.y();
          const Real z = r.z();

          const Real lx = LagrangeBasisSegment<K>::getBasis(i, x);
          const Real ly = LagrangeBasisSegment<K>::getBasis(j, y);
          const Real lz = LagrangeBasisSegment<K>::getBasis(k, z);

          const Real dlx = LagrangeBasisSegment<K>::getDerivative(i, x);
          const Real dly = LagrangeBasisSegment<K>::getDerivative(j, y);
          const Real dlz = LagrangeBasisSegment<K>::getDerivative(k, z);

          Real val = 0;
          if (m_i == 0)      // \partial/\partialx
            val = dlx * ly * lz;
          else if (m_i == 1) // \partial/\partialy
            val = lx * dly * lz;
          else if (m_i == 2) // \partial/\partialz
            val = lx * ly * dlz;

          return static_cast<Scalar>(val);
        }
      }

      return Scalar(0);
    }
    else
    {
      // Higher-order derivatives not implemented
      return Scalar(0);
    }
  }
}

/// @endcond
#endif
