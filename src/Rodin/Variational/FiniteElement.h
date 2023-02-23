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
#include "MFEM.h"

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

  template <class ... Ps>
  class FiniteElement<H1<Scalar, Ps...>>
  {
    public:
      using FES = H1<Scalar, Ps...>;

      FiniteElement(const Geometry::Simplex& simplex, const mfem::FiniteElement* handle)
        : m_simplex(simplex), m_handle(handle)
      {}

      /**
       * @brief Computes the basis functions.
       * @param[in] p Reference coordinates.
       *
       * Computes the vector @f$ \mathbf{\phi} (r) \in \mathbb{R}^{n} @f$ where
       * the i-th entry is equal to the value of the i-th basis function @f$
       * \phi_i(r) \in \mathbb{R} @f$ at the reference point @f$ r @f$, and @f$
       * n @f$ is the number of degrees of freedom.
       *
       * @returns Vector of dimension @f$ n @f$.
       */
      Math::Vector getBasis(const Math::Vector& p) const
      {
        Math::Vector basis(getDOFs());
        mfem::Vector tmp(basis.data(), basis.size());
        const mfem::IntegrationPoint ip = Internal::vec2ip(p);
        m_handle->CalcShape(ip, tmp);
        return basis;
      }

      /**
       * @brief Computes the gradient of each basis function.
       * @param[in] p Reference coordinates.
       *
       * Computes the matrix @f$ \nabla \mathbf{\phi} (r) \in \mathbb{R}^{s
       * \times n} @f$ defined by
       * @f[
       *  \nabla \phi (r) = \begin{bmatrix}
       *    \nabla \phi_1 (r) \ldots \nabla \phi_n (r)
       *  \end{bmatrix}
       * @f]
       * where @f$ s @f$ is the space dimension, @f$ n @f$ is the number of
       * degrees of freedom, and the column vector @f$ \nabla \phi_i (r) @f$ is
       * the gradient of the i-th the basis function at the reference point @f$
       * r @f$.
       *
       * @returns Matrix of dimension @f$ s \times n @f$.
       */
      Math::Matrix getGradient(const Math::Vector& p) const
      {
        Math::Matrix gradient(getDOFs(), getSimplex().getMesh().getSpaceDimension());
        mfem::DenseMatrix tmp;
        tmp.UseExternalData(gradient.data(), gradient.rows(), gradient.cols());
        const mfem::IntegrationPoint ip = Internal::vec2ip(p);
        m_handle->CalcDShape(ip, tmp);
        return gradient.transpose();
      }

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

