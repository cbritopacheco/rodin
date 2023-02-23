/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MFEM_H
#define RODIN_VARIATIONAL_MFEM_H

#include <mfem.hpp>

#include "Function.h"

namespace Rodin::Utility
{
  template <class T>
  inline
  mfem::Array<int> set2marker(const std::set<T>& s, size_t size)
  {
    mfem::Array<int> res(size);
    res = 0;
    for (const auto& v : s)
    {
      assert(v > 0);
      assert(v - 1 < size);
      res[static_cast<int>(v) - 1] = 1;
    }
    return res;
  }
}

namespace Rodin::Variational::Internal
{
  inline
  mfem::IntegrationPoint vec2ip(const Math::Vector& vec)
  {
    mfem::IntegrationPoint ip;
    ip.Set(vec.data(), vec.size());
    return ip;
  }

  inline
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Eigen::ColMajor, RODIN_MAXIMAL_SPACE_DIMENSION, 1>
  ip2vec(const mfem::IntegrationPoint& ip, const size_t dim)
  {
    if (dim == 1)
    {
      Math::FixedSizeVector<1> vec;
      vec << ip.x;
      return vec;
    }
    else if (dim == 2)
    {
      Math::FixedSizeVector<2> vec;
      vec << ip.x, ip.y;
      return vec;
    }
    else if (dim == 3)
    {
      Math::FixedSizeVector<3> vec;
      vec << ip.x, ip.y, ip.z;
      return vec;
    }
    else
    {
      assert(false);
      return Math::Vector::Zero(0);
    }
  }

  class MFEMElementTransformation final : public Geometry::SimplexTransformation
  {
    public:
      MFEMElementTransformation(mfem::ElementTransformation& trans)
        : m_handle(trans)
      {}

      Math::Vector transform(const Math::Vector& rc) const override
      {
        Math::Vector res(m_handle.get().GetSpaceDim());
        const mfem::IntegrationPoint ip = Internal::vec2ip(rc);
        m_handle.get().SetIntPoint(&ip);
        mfem::Vector tmp;
        tmp.SetDataAndSize(res.data(), res.size());
        assert(!tmp.OwnsData());
        m_handle.get().Transform(ip, tmp);
        return res;
      }

      Math::Matrix jacobian(const Math::Vector& rc) const override
      {
        Math::Matrix res(m_handle.get().GetSpaceDim(), m_handle.get().GetDimension());
        const mfem::IntegrationPoint ip = Internal::vec2ip(rc);
        m_handle.get().SetIntPoint(&ip);
        mfem::DenseMatrix tmp;
        tmp.UseExternalData(res.data(), res.rows(), res.cols());
        tmp = m_handle.get().Jacobian();
        return res;
      }

      mfem::ElementTransformation& getHandle() const override
      {
        return m_handle.get();
      }

    private:
      std::reference_wrapper<mfem::ElementTransformation> m_handle;
  };

  template <class Derived>
  class MFEMScalarCoefficient final : public mfem::Coefficient
  {
    public:
      MFEMScalarCoefficient(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
         : m_mesh(mesh), m_s(s)
      {}

      MFEMScalarCoefficient(const MFEMScalarCoefficient& other)
        : mfem::Coefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      MFEMScalarCoefficient(MFEMScalarCoefficient&& other)
        : mfem::Coefficient(std::move(other)),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      inline
      double Eval(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      override
      {
        Scalar res = 0;
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            res = m_s.get().getValue(Geometry::Point(
                  *m_mesh.get().getElement(trans.ElementNo),
                  MFEMElementTransformation(trans),
                  Internal::ip2vec(ip, trans.GetDimension())));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            res = m_s.get().getValue(Geometry::Point(
                  *m_mesh.get().getFace(faceIdx),
                  MFEMElementTransformation(trans),
                  Internal::ip2vec(ip, trans.GetDimension())));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            res = m_s.get().getValue(Geometry::Point(
                  *m_mesh.get().getFace(trans.ElementNo),
                  MFEMElementTransformation(trans),
                  Internal::ip2vec(ip, trans.GetDimension())));
            break;
          }
        }
        return res;
      }

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };

  template <class Derived>
  MFEMScalarCoefficient(const Geometry::MeshBase&, const FunctionBase<Derived>&)
    -> MFEMScalarCoefficient<Derived>;

  template <class Derived>
  class MFEMVectorCoefficient final : public mfem::VectorCoefficient
  {
    public:
      MFEMVectorCoefficient(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
        : mfem::VectorCoefficient(
            s.getRangeShape().height() == 1 ?
              s.getRangeShape().width() : s.getRangeShape().height()),
          m_mesh(mesh),
          m_s(s)
      {}

      MFEMVectorCoefficient(const MFEMVectorCoefficient& other)
        : mfem::VectorCoefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      MFEMVectorCoefficient(MFEMVectorCoefficient&& other)
        : mfem::VectorCoefficient(std::move(other)),
          m_mesh(other.mesh),
          m_s(other.m_s)
      {}

      inline
      void Eval(
          mfem::Vector& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      override
      {
        Math::Vector vec;
        assert(false);
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            // vec = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            // int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            // vec = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            // vec = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
            break;
          }
        }
        value.SetSize(vec.size());
        std::copy(vec.begin(), vec.end(), value.begin());
      }

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };

  template <class Derived>
  MFEMVectorCoefficient(const Geometry::MeshBase&, const FunctionBase<Derived>&)
    -> MFEMVectorCoefficient<Derived>;

  template <class Derived>
  class MFEMMatrixCoefficient final : public mfem::MatrixCoefficient
  {
    public:
      MFEMMatrixCoefficient(const Geometry::MeshBase& mesh, const FunctionBase<Derived>& s)
        : mfem::MatrixCoefficient(s.getRangeShape().height(), s.getRangeShape().width()),
          m_mesh(mesh),
          m_s(s)
      {}

      MFEMMatrixCoefficient(const MFEMMatrixCoefficient& other)
        : mfem::MatrixCoefficient(other),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      MFEMMatrixCoefficient(MFEMMatrixCoefficient&& other)
        : mfem::MatrixCoefficient(std::move(other)),
          m_mesh(other.m_mesh),
          m_s(other.m_s)
      {}

      inline
      void Eval(
          mfem::DenseMatrix& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      override
      {
        Math::Matrix mat;
        assert(false);
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            // mat = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            // int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            // mat = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            // mat = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
            break;
          }
        }
        value.SetSize(mat.rows(), mat.cols());
        value = mat.data();
      }
    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      std::reference_wrapper<const FunctionBase<Derived>> m_s;
  };

  template <class Derived>
  MFEMMatrixCoefficient(const Geometry::MeshBase&, const FunctionBase<Derived>&)
    -> MFEMMatrixCoefficient<Derived>;
}

#endif
