/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H
#define RODIN_GEOMETRY_ISOPARAMETRICTRANSFORMATION_H

#include "Rodin/Variational/MFEM.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "SimplexTransformation.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents a transformation between the reference space to the
   * physical space of a simplex.
   */
  class IsoparametricTransformation final : public SimplexTransformation
  {
    public:
      IsoparametricTransformation(mfem::IsoparametricTransformation* trans)
        : m_handle(trans)
      {}

      // inline
      // Math::Vector transform(const Math::Vector& rc) const final override
      // {
      //   const mfem::IntegrationPoint ip = Variational::Internal::vec2ip(rc);
      //   Math::Vector shape(m_handle->GetFE()->GetDof());
      //   mfem::Vector tmp(shape.data(), shape.size());
      //   m_handle->GetFE()->CalcShape(ip, tmp);
      //   return m_pm * shape;
      // }

      // inline
      // Math::Matrix jacobian(const Math::Vector& rc) const final override
      // {
      //   const mfem::IntegrationPoint ip = Variational::Internal::vec2ip(rc);
      //   Math::Matrix dshape(m_handle->GetFE()->GetDof(), m_handle->GetFE()->GetDim());
      //   mfem::DenseMatrix tmp;
      //   tmp.UseExternalData(dshape.data(), dshape.rows(), dshape.cols());
      //   m_handle->GetFE()->CalcDShape(ip, tmp);
      //   return m_pm * dshape;
      // }

      inline
      mfem::IsoparametricTransformation& getHandle() const final override
      {
        return *m_handle;
      }

    private:
      std::unique_ptr<mfem::IsoparametricTransformation> m_handle;
      // Eigen::Map<Math::Matrix> m_pm;
  };
}

#endif
