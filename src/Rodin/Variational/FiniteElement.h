/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include "Rodin/Math/DenseMatrix.h"
#include "Rodin/Geometry/Simplex.h"

#include "ForwardDecls.h"

#include <mfem.hpp>

namespace Rodin::Variational
{
  class FiniteElement
  {
    public:
      FiniteElement(
          const Geometry::Simplex& simplex,
          const mfem::FiniteElement* handle)
        : m_simplex(simplex), m_handle(handle)
      {}

      Math::DenseMatrix getBasis(const Math::Vector& p) const;

      Math::DenseMatrix getJacobian(const Math::Vector& p) const;

      Math::DenseMatrix getDivergence(const Math::Vector& p) const;

      Math::DenseMatrix getCurl(const Math::Vector& p) const;

      Math::DenseMatrix getHessian(const Math::Vector& p) const;

      const Geometry::Simplex& getSimplex() const
      {
        return m_simplex.get();
      }

      const mfem::FiniteElement& getHandle() const
      {
        return *m_handle;
      }

    private:
      std::reference_wrapper<const Geometry::Simplex> m_simplex;
      const mfem::FiniteElement* m_handle;
  };
}

#endif

