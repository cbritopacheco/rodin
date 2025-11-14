/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Pk.h"

namespace Rodin::Variational
{
  // Explicit template instantiations for common cases
  // K=0 (P0)
  template class Pk<0, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<0, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<0, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<0, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=1 (P1)
  template class Pk<1, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<1, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<1, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<1, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=2 (P2)
  template class Pk<2, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<2, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<2, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<2, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=3 (P3)
  template class Pk<3, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<3, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<3, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<3, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=4 (P4)
  template class Pk<4, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<4, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<4, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<4, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=5 (P5)
  template class Pk<5, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<5, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<5, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<5, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;

  // K=6 (P6)
  template class Pk<6, Real, Geometry::Mesh<Context::Local>>;
  template class Pk<6, Complex, Geometry::Mesh<Context::Local>>;
  template class Pk<6, Math::Vector<Real>, Geometry::Mesh<Context::Local>>;
  template class Pk<6, Math::Vector<Complex>, Geometry::Mesh<Context::Local>>;
}
