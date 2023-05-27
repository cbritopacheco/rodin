/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

#include <Eigen/QR>

#include "Rodin/Variational/P1.h"

#include "SimplexTransformation.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  template <class FE>
  class IsoparametricTransformation final : public PolytopeTransformation
  {
    public:
      IsoparametricTransformation(Math::Matrix&& pm, FE&& fe)
        : m_pm(std::move(pm)), m_fe(std::move(fe))
      {}

      /**
       * pm : sdim x dof
       */
      IsoparametricTransformation(const Math::Matrix& pm, const FE& fe)
        : m_pm(pm), m_fe(fe)
      {}

      IsoparametricTransformation(Math::Matrix&& pm, const FE& fe)
        : m_pm(std::move(pm)), m_fe(fe)
      {}

      IsoparametricTransformation(Math::Matrix&& pm)
        : m_pm(std::move(pm))
      {}

      inline
      Math::Vector transform(const Math::Vector& rc) const final override
      {
        return m_pm * m_fe.getBasis(rc);
      }

      inline
      Math::Matrix jacobian(const Math::Vector& rc) const final override
      {
        return m_pm * m_fe.getGradient(rc);
      }

    private:
      Math::Matrix m_pm;
      FE m_fe;
  };
}

#endif
