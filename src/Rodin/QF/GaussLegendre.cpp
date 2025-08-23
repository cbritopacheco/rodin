#include "GaussLegendre.h"

namespace Rodin::QF
{
  void GaussLegendre::gl_1d_unit(size_t n, std::vector<Real>& x, std::vector<Real>& w, size_t maxIt, Real tol)
  {
    assert(n >= 1);
    std::vector<Real> xi(n), wi(n);
    const size_t m = (n + 1) / 2;

    for (size_t k = 0; k < m; ++k)
    {
      const Real kk = static_cast<Real>(k) + 0.75;
      Real xk = std::cos(M_PI * kk / (static_cast<Real>(n) + 0.5));

      for (size_t it = 0; it < maxIt; ++it)
      {
        Real p0 = 1.0;
        Real p1 = xk;
        for (size_t j = 2; j <= n; ++j)
        {
          const Real jj = static_cast<Real>(j);
          const Real p = ((2 * jj - 1) * xk * p1 - (jj - 1) * p0) / jj;
          p0 = p1;
          p1 = p;
        }
        const Real pn = p1;
        const Real pnm1 = p0;
        const Real dpn = (static_cast<Real>(n) * (xk * pn - pnm1)) / (xk * xk - 1.0);
        const Real dx = pn / dpn;
        xk -= dx;
        if (std::abs(dx) < tol) break;
      }

      Real p0 = 1.0;
      Real p1 = xk;
      for (size_t j = 2; j <= n; ++j)
      {
        const Real jj = static_cast<Real>(j);
        const Real p = ((2 * jj - 1) * xk * p1 - (jj - 1) * p0) / jj;
        p0 = p1;
        p1 = p;
      }
      const Real pn = p1;
      const Real pnm1 = p0;
      const Real dpn = (static_cast<Real>(n) * (xk * pn - pnm1)) / (xk * xk - 1.0);
      const Real wv = 2.0 / ((1.0 - xk * xk) * dpn * dpn);

      xi[k]         = -xk;
      xi[n - 1 - k] =  xk;
      wi[k]         =  wv;
      wi[n - 1 - k] =  wv;
    }

    x.resize(n);
    w.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
      x[i] = 0.5 * (xi[i] + 1.0);
      w[i] = 0.5 * wi[i];
    }
  }

  void GaussLegendre::build()
  {
    switch (getGeometry())
    {
      case Geometry::Polytope::Type::Point:         build_point(); break;
      case Geometry::Polytope::Type::Segment:       build_segment(m_nx); break;
      case Geometry::Polytope::Type::Quadrilateral: build_quad(m_nx, m_ny); break;
      case Geometry::Polytope::Type::Triangle:      build_tri(m_nx, m_nx); break;
      case Geometry::Polytope::Type::Tetrahedron:   build_tet(m_nx, m_ny, m_nz); break;
      case Geometry::Polytope::Type::Wedge:         build_wedge(m_nx, m_ny); break;
    }
  }

  void GaussLegendre::build_point()
  {
    m_points.clear();
    m_weights.resize(1);
    Math::SpatialVector<Real> p;
    p.resize(1);
    p[0] = 0.0;
    m_points.push_back(std::move(p));
    m_weights[0] = 1.0;
  }

  void GaussLegendre::build_segment(size_t n)
  {
    std::vector<Real> x, w;
    gl_1d_unit(n, x, w);
    m_points.clear();
    m_points.reserve(n);
    m_weights.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
      Math::SpatialVector<Real> p;
      p.resize(1);
      p[0] = x[i];
      m_points.push_back(std::move(p));
      m_weights[i] = w[i];
    }
  }

  void GaussLegendre::build_quad(size_t nx, size_t ny)
  {
    std::vector<Real> x, wx, y, wy;
    gl_1d_unit(nx, x, wx);
    gl_1d_unit(ny, y, wy);

    const size_t N = nx * ny;
    m_points.clear();
    m_points.reserve(N);
    m_weights.resize(N);

    size_t k = 0;
    for (size_t j = 0; j < ny; ++j)
    {
      for (size_t i = 0; i < nx; ++i)
      {
        Math::SpatialVector<Real> p;
        p.resize(2);
        p[0] = x[i];
        p[1] = y[j];
        m_points.push_back(std::move(p));
        m_weights[k++] = wx[i] * wy[j];
      }
    }
  }

  void GaussLegendre::build_tri(size_t nu, size_t nv)
  {
    std::vector<Real> u, wu, v, wv;
    gl_1d_unit(nu, u, wu);
    gl_1d_unit(nv, v, wv);

    const size_t N = nu * nv;
    m_points.clear();
    m_points.reserve(N);
    m_weights.resize(N);

    size_t k = 0;
    for (size_t j = 0; j < nv; ++j)
    {
      for (size_t i = 0; i < nu; ++i)
      {
        const Real r = u[i];
        const Real s = (1.0 - u[i]) * v[j];
        Math::SpatialVector<Real> p;
        p.resize(2);
        p[0] = r;
        p[1] = s;
        m_points.push_back(std::move(p));
        m_weights[k++] = wu[i] * wv[j] * (1.0 - u[i]);
      }
    }
  }

  void GaussLegendre::build_tet(size_t nu, size_t nv, size_t nw)
  {
    std::vector<Real> u, wu, v, wv, w, ww;
    gl_1d_unit(nu, u, wu);
    gl_1d_unit(nv, v, wv);
    gl_1d_unit(nw, w, ww);

    const size_t N = nu * nv * nw;
    m_points.clear();
    m_points.reserve(N);
    m_weights.resize(N);

    size_t k = 0;
    for (size_t k3 = 0; k3 < nw; ++k3)
    {
      for (size_t k2 = 0; k2 < nv; ++k2)
      {
        for (size_t k1 = 0; k1 < nu; ++k1)
        {
          const Real r = u[k1];
          const Real s = (1.0 - u[k1]) * v[k2];
          const Real t = (1.0 - u[k1]) * (1.0 - v[k2]) * w[k3];

          Math::SpatialVector<Real> p;
          p.resize(3);
          p[0] = r;
          p[1] = s;
          p[2] = t;
          m_points.push_back(std::move(p));

          m_weights[k++] = wu[k1] * wv[k2] * ww[k3]
                         * (1.0 - u[k1]) * (1.0 - u[k1]) * (1.0 - v[k2]);
        }
      }
    }
  }

  void GaussLegendre::build_wedge(size_t ntri, size_t nz)
  {
    std::vector<Real> u, wu, v, wv, z, wz;
    gl_1d_unit(ntri, u, wu);
    gl_1d_unit(ntri, v, wv);
    gl_1d_unit(nz, z, wz);

    const size_t N = ntri * ntri * nz;
    m_points.clear();
    m_points.reserve(N);
    m_weights.resize(N);

    size_t k = 0;
    for (size_t kk = 0; kk < nz; ++kk)
    {
      for (size_t j = 0; j < ntri; ++j)
      {
        for (size_t i = 0; i < ntri; ++i)
        {
          const Real r = u[i];
          const Real s = (1.0 - u[i]) * v[j];
          const Real t = z[kk];

          Math::SpatialVector<Real> p;
          p.resize(3);
          p[0] = r;
          p[1] = s;
          p[2] = t;
          m_points.push_back(std::move(p));

          m_weights[k++] = wu[i] * wv[j] * (1.0 - u[i]) * wz[kk];
        }
      }
    }
  }
}
