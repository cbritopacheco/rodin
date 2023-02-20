/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Simplex.h"
#include "Rodin/Geometry/SimplexTransformation.h"

#include "ForwardDecls.h"

#include <mfem.hpp>

namespace Rodin::Variational
{
  class FiniteElement
  {
    public:
      FiniteElement(const Geometry::Simplex& simplex, const Geometry::Transformation& trans,
                    const mfem::FiniteElement* handle)
        : m_simplex(simplex), m_trans(trans), m_handle(handle)
      {}

      Math::Matrix getBasis(const Math::Vector& p) const;

      Math::Matrix getJacobian(const Math::Vector& p) const;

      Math::Matrix getDivergence(const Math::Vector& p) const;

      Math::Matrix getCurl(const Math::Vector& p) const;

      Math::Matrix getHessian(const Math::Vector& p) const;

      inline
      const Geometry::Simplex& getSimplex() const
      {
        return m_simplex.get();
      }

      inline
      const Geometry::Transformation& getTransformation() const
      {
        return m_trans.get();
      }

      inline
      const mfem::FiniteElement& getHandle() const
      {
        return *m_handle;
      }

    private:
      std::reference_wrapper<const Geometry::Simplex> m_simplex;
      std::reference_wrapper<const Geometry::Transformation> m_trans;
      const mfem::FiniteElement* m_handle;
  };
}

#endif

