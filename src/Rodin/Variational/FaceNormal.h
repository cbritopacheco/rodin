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
        const auto& polytope = p.getPolytope();
        const auto& mesh = polytope.getMesh();
        const auto& jacobian = p.getJacobian();
        assert(polytope.getDimension() == mesh.getDimension() - 1);
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
        return value.normalized();
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
      FaceNormal& traceOf(const FlatSet<Geometry::Attribute>& attrs)
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

