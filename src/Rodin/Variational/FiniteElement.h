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

#include "TensorBasis.h"

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
  class FiniteElement<H1<Scalar, Ps...>> final
  {
    public:
      using FES = H1<Scalar, Ps...>;

      constexpr
      FiniteElement(const Geometry::Simplex& simplex, const mfem::FiniteElement* handle)
        : m_simplex(simplex), m_handle(handle)
      {}

      constexpr
      FiniteElement(const FiniteElement&) = default;

      constexpr
      FiniteElement(FiniteElement&&) = default;

      /**
       * @f$ n \times 1 @f$
       */
      Math::Vector getBasis(const Math::Vector& r) const
      {
        assert(r.size() == getSimplex().getDimension());
        Math::Vector basis(getDOFs());
        mfem::Vector tmp(basis.data(), basis.size());
        m_handle->CalcShape(Internal::vec2ip(r), tmp);
        return basis;
      }

      /**
       * @f$ n \times r @f$
       */
      Math::Matrix getGradient(const Math::Vector& r) const
      {
        assert(r.size() == getSimplex().getDimension());
        Math::Matrix gradient(getDOFs(), getSimplex().getDimension());
        mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
        m_handle->CalcDShape(Internal::vec2ip(r), tmp);
        return gradient;
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

  template <class ... Ps>
  class FiniteElement<H1<Math::Vector, Ps...>> final
  {
    public:
      using FES = H1<Math::Vector, Ps...>;

      constexpr
      FiniteElement(size_t vdim, const Geometry::Simplex& simplex, const mfem::FiniteElement* handle)
        : m_vdim(vdim), m_simplex(simplex), m_handle(handle)
      {}

      constexpr
      FiniteElement(const FiniteElement&) = default;

      constexpr
      FiniteElement(FiniteElement&&) = default;

      /**
       * @f$ nd \times d @f$
       */
      Math::Matrix getBasis(const Math::Vector& r) const
      {
        assert(r.size() == getSimplex().getDimension());
        const size_t n = m_handle->GetDof();
        Math::Matrix basis = Math::Matrix::Zero(getDOFs(), m_vdim);
        Math::Vector shape(n);
        mfem::Vector tmp(shape.data(), shape.size());
        m_handle->CalcShape(Internal::vec2ip(r), tmp);
        for (size_t i = 0; i < m_vdim; i++)
          basis.block(i * n, i, n, 1) = shape;
        return basis;
      }

      /**
       * @f$ nd \times 1 @f$
       */
      Math::Vector getDivergence(const Math::Vector& r) const
      {
        assert(r.size() == getSimplex().getDimension());
        const size_t n = m_handle->GetDof();
        Math::Matrix gradient(n, r.size());
        mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
        m_handle->CalcDShape(Internal::vec2ip(r), tmp);
        return gradient.reshaped();
      }

      /**
       * @f$ nd \times d \times r @f$
       */
      Math::Tensor<3> getJacobian(const Math::Vector& r) const
      {
        assert(m_vdim == r.size());
        assert(r.size() == getSimplex().getDimension());
        const size_t n = m_handle->GetDof();
        const size_t rdim = r.size();
        Math::Matrix gradient(n, rdim);
        mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
        m_handle->CalcDShape(Internal::vec2ip(r), tmp);
        Math::Tensor<3> jacobian(getDOFs(), Math::Tensor<3>::Index(m_vdim), Math::Tensor<3>::Index(rdim));
        jacobian.setZero();
        for (size_t i = 0; i < n; i++)
        {
          for (size_t j = 0; j < m_vdim; j++)
          {
            Math::Slice<Math::Tensor<3>, 0>(jacobian, Eigen::Index(i + j * n)).row(Eigen::Index(j))
              = gradient.row(Eigen::Index(i));
          }
        }
        return jacobian;
      }

      inline
      size_t getOrder() const
      {
        return m_handle->GetOrder();
      }

      inline
      size_t getDOFs() const
      {
        return m_vdim * m_handle->GetDof();
      }

      inline
      const Geometry::Simplex& getSimplex() const
      {
        return m_simplex.get();
      }

      inline
      constexpr
      size_t getVectorDimension() const
      {
        return m_vdim;
      }

      const mfem::FiniteElement& getHandle() const
      {
        return *m_handle;
      }

    private:
      const size_t m_vdim;
      std::reference_wrapper<const Geometry::Simplex> m_simplex;
      const mfem::FiniteElement* m_handle;
  };
}

#endif

