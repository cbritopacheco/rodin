#ifndef RODIN_VARIATIONAL_FACENORMAL_H
#define RODIN_VARIATIONAL_FACENORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/SimplexTransformation.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal on a face.
   */
  class FaceNormal : public VectorFunctionBase<FaceNormal>
  {
    public:
      using Parent = VectorFunctionBase<FaceNormal>;

      /**
       * @brief Constructs the outward unit on a face.
       */
      FaceNormal(const Geometry::MeshBase& surface)
        : m_dimension(surface.getSpaceDimension())
      {
        assert(m_dimension > 0);
      }

      FaceNormal(const FaceNormal& other)
        : Parent(other),
          m_dimension(other.m_dimension)
      {}

      FaceNormal(FaceNormal&& other)
        : Parent(std::move(other)),
          m_dimension(std::move(other.m_dimension))
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      Math::Vector getValue(const Geometry::Point& p) const
      {
        assert(p.getSimplex().getDimension() == p.getSimplex().getMesh().getSpaceDimension() - 1);
        const auto& simplex = p.getSimplex();
        const auto& mesh = simplex.getMesh();
        const auto& jacobian = p.getJacobian();

        Math::Vector value(m_dimension);
        if (jacobian.rows() == 2)
        {
          value(0) =  jacobian(1, 0);
          value(1) = -jacobian(0, 0);
        }
        else if (jacobian.rows() == 3)
        {
          value(0) = jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1);
          value(1) = jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1);
          value(2) = jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
        }
        else
        {
          assert(false);
          value.setConstant(NAN);
          return value;
        }

        // Or we are on a face of a d-mesh in d-space
        if (mesh.isBoundary(simplex.getIndex()))
        {
          return value.normalized();
        }
        else
        {
          // int el1 = -1;
          // int el2 = -1;
          assert(false);
          // const auto& meshHandle = simplex.getMesh().getHandle();
          // meshHandle.GetFaceElements(simplex.getIndex(), &el1, &el2);
          // if (el1 >= 0 && getTraceDomain().count(meshHandle.GetAttribute(el1)))
          // {
          //   return value.normalized();
          // }
          // else if (el2 >= 0 && getTraceDomain().count(meshHandle.GetAttribute(el2)))
          // {
          //   value = -1.0 * value;
          //   return value.normalized();
          // }
          // else
          // {
          //   assert(false);
          //   value.setConstant(NAN);
          //   return value;
          // }
        }
      }

      inline
      constexpr
      FaceNormal& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        return *this;
      }

      inline
      constexpr
      FaceNormal& traceOf(const std::set<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        return *this;
      }

      inline FaceNormal* copy() const noexcept override
      {
        return new FaceNormal(*this);
      }

    private:
      const size_t m_dimension;
  };
}

#endif

