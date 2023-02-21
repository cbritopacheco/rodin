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
  class FiniteElementOrder
  {
    public:
      explicit
      constexpr
      FiniteElementOrder(size_t v)
        : m_v(v)
      {}

      constexpr
      FiniteElementOrder(const FiniteElementOrder&) = default;

      constexpr
      FiniteElementOrder(FiniteElementOrder&&) = default;

      inline
      constexpr
      operator size_t() const
      {
        return m_v;
      }

    private:
      const size_t m_v;
  };


  class FiniteElement
  {
    public:
      FiniteElement(const Geometry::Simplex& simplex, const mfem::FiniteElement* handle)
        : m_simplex(simplex), m_handle(handle)
      {}

      Math::Matrix getBasis(const Math::Vector& p) const;

      Math::Matrix getGradient(const Math::Vector& p) const;

      Math::Matrix getDivergence(const Math::Vector& p) const;

      Math::Matrix getCurl(const Math::Vector& p) const;

      Math::Matrix getHessian(const Math::Vector& p) const;

      inline
      size_t getOrder() const
      {
        return m_handle->GetOrder();
      }

      inline
      size_t getDOFs() const
      {
        return m_handle->GetDof();
      }

      inline
      const Geometry::Simplex& getSimplex() const
      {
        return m_simplex.get();
      }

      inline
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

