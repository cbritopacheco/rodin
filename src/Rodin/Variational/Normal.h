#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class Normal : public VectorFunctionBase
  {
    public:
      /**
       * @brief Constructs the outward unit normal.
       */
      Normal(size_t dimension)
        : m_dimension(dimension)
      {
        assert(dimension > 0);
      }

      Normal(const Normal& other)
        :  VectorFunctionBase(other),
          m_dimension(other.m_dimension)
      {}

      Normal(Normal&& other)
        :  VectorFunctionBase(std::move(other)),
          m_dimension(other.m_dimension)
      {}

      int getDimension() const override
      {
        return m_dimension;
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        FunctionValue::Vector value;
        const auto& simplex = p.getSimplex();
        const auto& mesh = simplex.getMesh();
        auto& trans = simplex.getTransformation();

        if ((mesh.getDimension() + 1 == mesh.getSpaceDimension()) &&
            simplex.getDimension() == mesh.getDimension())
        {
          // Or we are on an element of a d-mesh in (d + 1)-space.
          value.SetSize(m_dimension);
          mfem::CalcOrtho(trans.Jacobian(), value);
          const double norm = value.Norml2();
          assert(norm > 0.0);
          assert(std::isfinite(norm));
          value /= norm;
        }
        else if ((mesh.getDimension() == mesh.getSpaceDimension()) &&
            simplex.getDimension() == mesh.getDimension() - 1)
        {
          // We are on a face of a d-mesh in d-space
          if (mesh.isBoundary(simplex.getIndex()))
          {
            mfem::FaceElementTransformations* ft =
              const_cast<Geometry::MeshBase&>(mesh).getHandle()
              .GetFaceElementTransformations(simplex.getIndex());
            value.SetSize(m_dimension);
            ft->SetAllIntPoints(&p.getIntegrationPoint());
            assert(ft->Elem1);
            mfem::CalcOrtho(trans.Jacobian(), value);
            const double norm = value.Norml2();
            assert(norm > 0.0);
            assert(std::isfinite(norm));
            value /= norm;
          }
          else
          {
            mfem::FaceElementTransformations* ft =
              const_cast<Geometry::MeshBase&>(mesh).getHandle()
              .GetFaceElementTransformations(simplex.getIndex());
            value.SetSize(m_dimension);
            ft->SetAllIntPoints(&p.getIntegrationPoint());
            if (ft->Elem1 && getTraceDomain() == ft->Elem1->Attribute)
            {
              mfem::CalcOrtho(trans.Jacobian(), value);
              const double norm = value.Norml2();
              assert(norm > 0.0);
              assert(std::isfinite(norm));
              value /= norm;
            }
            else if (ft->Elem2 && getTraceDomain() == ft->Elem2->Attribute)
            {
              mfem::CalcOrtho(trans.Jacobian(), value);
              const double norm = value.Norml2();
              assert(norm > 0.0);
              assert(std::isfinite(norm));
              value /= -norm;
            }
            else
            {
              assert(false);
            }
          }
        }
        else
        {
          assert(false);
          return FunctionValue::Vector{0, 0};
        }
        return value;
      }

      Normal* copy() const noexcept override
      {
        return new Normal(*this);
      }
    private:
      const size_t m_dimension;
  };
}

#endif
