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
  inline
  mfem::Array<int> set2marker(const std::set<int>& s, int size)
  {
    mfem::Array<int> res(size);
    res = 0;
    for (const auto& v : s)
    {
      assert(v > 0);
      assert(v - 1 < size);
      res[v - 1] = 1;
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

  template <class Derived>
  class MFEMScalarCoefficient : public mfem::Coefficient
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
      final override
      {
        Scalar res = 0;
        assert(false);
        switch (trans.ElementType)
        {
          case mfem::ElementTransformation::ELEMENT:
          {
            // res = m_s.getValue(Geometry::Point(*m_mesh.get().getElement(trans.ElementNo), ip));
            break;
          }
          case mfem::ElementTransformation::BDR_ELEMENT:
          {
            // int faceIdx = m_mesh.get().getHandle().GetBdrFace(trans.ElementNo);
            // res = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(faceIdx), ip));
            break;
          }
          case mfem::ElementTransformation::FACE:
          {
            // res = m_s.getValue(Geometry::Point(*m_mesh.get().getFace(trans.ElementNo), ip));
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
  class MFEMVectorCoefficient : public mfem::VectorCoefficient
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
      final override
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
  class MFEMMatrixCoefficient : public mfem::MatrixCoefficient
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
      final override
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
}

#endif
