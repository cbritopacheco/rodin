/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include <unordered_map>

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
      FiniteElement(const Geometry::Polytope& simplex, const mfem::FiniteElement* handle)
        : m_rdim(simplex.getDimension()), m_handle(handle)
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
        assert(r.size() == getDimension());
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
        assert(r.size() == getDimension());
        Math::Matrix gradient(getDOFs(), getDimension());
        mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
        m_handle->CalcDShape(Internal::vec2ip(r), tmp);
        return gradient;
      }

      inline
      constexpr
      size_t getDimension() const
      {
        return m_rdim;
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

      const mfem::FiniteElement& getHandle() const
      {
        return *m_handle;
      }

    private:
      const size_t m_rdim;
      const mfem::FiniteElement* m_handle;
  };

  template <class ... Ps>
  class FiniteElement<H1<Math::Vector, Ps...>> final
  {
    struct Cache
    {
      using Key = const Math::Vector*;
      std::unordered_map<Key, Math::Matrix> basis;
      std::unordered_map<Key, Math::Vector> divergence;
      std::unordered_map<Key, Math::Tensor<3>> jacobian;
    };

    public:
      using FES = H1<Math::Vector, Ps...>;

      constexpr
      FiniteElement(size_t vdim, const Geometry::Polytope& simplex, const mfem::FiniteElement* handle)
        : m_rdim(simplex.getDimension()), m_vdim(vdim), m_handle(handle)
      {}

      constexpr
      FiniteElement(const FiniteElement&) = default;

      constexpr
      FiniteElement(FiniteElement&&) = default;

      /**
       * @f$ nd \times d @f$
       */
      const Math::Matrix& getBasis(const Math::Vector& r) const
      {
        assert(r.size() == getDimension());
        auto search = m_cache.basis.find(&r);
        if (search != m_cache.basis.end())
        {
          return search->second;
        }
        else
        {
          const size_t n = m_handle->GetDof();
          auto [it, inserted] =
            m_cache.basis.emplace(std::make_pair(&r, Math::Matrix(getDOFs(), m_vdim)));
          assert(inserted);
          Math::Matrix& basis = it->second;
          basis.setZero();
          Math::Vector shape(n);
          mfem::Vector tmp(shape.data(), shape.size());
          m_handle->CalcShape(Internal::vec2ip(r), tmp);
          for (size_t i = 0; i < m_vdim; i++)
            basis.block(i * n, i, n, 1) = shape;
          return basis;
        }
      }

      /**
       * @f$ nd \times 1 @f$
       *
       * @f[
       * \begin{pmatrix}
       *  \partial_{x_1} \phi_1\\
       *  \vdots\\
       *  \partial_{x_1} \phi_n\\
       *  \vdots\\
       *  \partial_{x_r} \phi_1\\
       *  \vdots\\
       *  \partial_{x_r} \phi_n
       * \end{pmatrix}
       *
       * @f]
       */
      const Math::Vector& getDivergence(const Math::Vector& r) const
      {
        assert(r.size() == getDimension());
        auto search = m_cache.divergence.find(&r);
        if (search != m_cache.divergence.end())
        {
          return search->second;
        }
        else
        {
          const size_t n = m_handle->GetDof();
          Math::Matrix gradient(n, r.size());
          mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
          m_handle->CalcDShape(Internal::vec2ip(r), tmp);
          auto [it, inserted] =
            m_cache.divergence.emplace(std::make_pair(&r, Math::Vector(getDOFs())));
          assert(inserted);
          Math::Vector& divergence = it->second;
          divergence = gradient.reshaped();
          return divergence;
        }
      }

      /**
       * @f$ nd \times d \times r @f$
       */
      const Math::Tensor<3>& getJacobian(const Math::Vector& r) const
      {
        assert(m_vdim == r.size());
        assert(r.size() == getDimension());

        auto search = m_cache.jacobian.find(&r);
        if (search != m_cache.jacobian.end())
        {
          return search->second;
        }
        else
        {
          const size_t n = m_handle->GetDof();
          const size_t rdim = r.size();
          Math::Matrix gradient(n, rdim);
          mfem::DenseMatrix tmp(gradient.data(), gradient.rows(), gradient.cols());
          m_handle->CalcDShape(Internal::vec2ip(r), tmp);
          auto [it, inserted] =
            m_cache.jacobian.emplace(std::make_pair(
                  &r, Math::Tensor<3>(getDOFs(), Math::Tensor<3>::Index(m_vdim), Math::Tensor<3>::Index(rdim))));
          assert(inserted);
          Math::Tensor<3>& jacobian = it->second;
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
      }

      inline
      size_t getComponentDOFs() const
      {
        return m_handle->GetDof();
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
      size_t getDimension() const
      {
        return m_rdim;
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
      const size_t m_rdim;
      const size_t m_vdim;
      const mfem::FiniteElement* m_handle;
      mutable Cache m_cache;
  };
}

#endif

