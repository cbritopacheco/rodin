/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

#include <Eigen/QR>

#include "Rodin/Geometry/Simplex.h"

#include "SimplexTransformation.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Polytope isoparametric transformation.
   */
  template <class FE>
  class IsoparametricTransformation final : public PolytopeTransformation
  {
    static_assert(std::is_same_v<typename FE::RangeType, Scalar>,
        "Type of finite element must be scalar valued.");

    public:
      using Parent = PolytopeTransformation;

      IsoparametricTransformation(Math::Matrix&& pm, FE&& fe)
        : m_pm(std::move(pm)),
          m_sdim(m_pm.rows()),
          m_fe(std::move(fe)),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      /**
       * pm : sdim x dof
       */
      IsoparametricTransformation(const Math::Matrix& pm, const FE& fe)
        : m_pm(pm),
          m_sdim(m_pm.rows()),
          m_fe(fe),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(Math::Matrix&& pm, const FE& fe)
        : m_pm(std::move(pm)),
          m_sdim(m_pm.rows()),
          m_fe(fe),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const Math::Matrix& pm, FE&& fe)
        : m_pm(pm),
          m_sdim(m_pm.rows()),
          m_fe(std::move(fe)),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(const IsoparametricTransformation& other)
        : Parent(other),
          m_fe(other.m_fe),
          m_sdim(m_pm.rows()),
          m_pm(other.m_pm),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      IsoparametricTransformation(IsoparametricTransformation&& other)
        : Parent(std::move(other)),
          m_fe(std::move(other.m_fe)),
          m_sdim(m_pm.rows()),
          m_pm(std::move(other.m_pm)),
          m_rdim(Polytope::getGeometryDimension(m_fe.getGeometry()))
      {
        assert(m_pm.cols() >= 0);
        assert(static_cast<size_t>(m_pm.cols()) == m_fe.getCount());
      }

      inline
      Math::Vector transform(const Math::Vector& rc) const override
      {
        Math::Vector res = Math::Vector::Zero(m_sdim);
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          assert(res.size() == m_pm.col(local).size());
          res.noalias() += m_pm.col(local) * m_fe.getBasis(local)(rc);
        }
        return res;
      }

      inline
      const Math::Vector& transform(CacheResultType, const Math::Vector& rc) const override
      {
        Math::Vector res = Math::Vector::Zero(m_sdim);
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          assert(res.size() == m_pm.col(local).size());
          res.noalias() += m_pm.col(local) * m_fe.getBasis(local)(CacheResult, rc);
        }

        auto it = m_transform.find(&rc);
        if (it == m_transform.end())
        {
          auto rit = m_transform.insert(it, { &rc, std::move(res) });
          return rit->second;
        }
        else
        {
          return it->second;
        }
      }

      inline
      Math::Matrix jacobian(const Math::Vector& rc) const override
      {
        Math::Matrix res = Math::Matrix::Zero(m_sdim, m_rdim);
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          const auto& gradient = m_fe.getGradient(local)(rc);
          for (size_t i = 0; i < m_rdim; i++)
          {
            assert(res.col(i).size() == m_pm.col(local).size());
            res.col(i).noalias() += m_pm.col(local) * gradient.coeff(i);
          }
        }
        return res;
      }

      inline
      const Math::Matrix& jacobian(CacheResultType, const Math::Vector& rc) const override
      {
        Math::Matrix res = Math::Matrix::Zero(m_sdim, m_rdim);
        for (size_t local = 0; local < m_fe.getCount(); local++)
        {
          const auto& gradient = m_fe.getGradient(local)(CacheResult, rc);
          for (size_t i = 0; i < m_rdim; i++)
          {
            assert(res.col(i).size() == m_pm.col(local).size());
            res.col(i).noalias() += m_pm.col(local) * gradient.coeff(i);
          }
        }

        auto it = m_jacobian.find(&rc);
        if (it == m_jacobian.end())
        {
          auto rit = m_jacobian.insert(it, { &rc, std::move(res) });
          return rit->second;
        }
        else
        {
          return it->second;
        }
      }

      inline
      const Math::Matrix& getPointMatrix() const
      {
        return m_pm;
      }

    private:
      Math::Matrix m_pm;
      const size_t m_sdim;
      FE m_fe;
      const size_t m_rdim;
  };
}

#endif
