#ifndef RODIN_VARIATIONAL_BOUNDARYNORMAL_H
#define RODIN_VARIATIONAL_BOUNDARYNORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/SimplexTransformation.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class BoundaryNormal final : public VectorFunctionBase<BoundaryNormal>
  {
    public:
      using Parent = VectorFunctionBase<BoundaryNormal>;

      /**
       * @brief Constructs the outward unit normal.
       */
      constexpr
      BoundaryNormal(const Geometry::MeshBase& mesh)
        : m_dimension(mesh.getSpaceDimension())
      {
        assert(m_dimension > 0);
      }

      constexpr
      BoundaryNormal(const BoundaryNormal& other)
        : Parent(other),
          m_dimension(other.m_dimension)
      {}

      constexpr
      BoundaryNormal(BoundaryNormal&& other)
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
        const auto& jacobian = p.getJacobian();
        Math::Vector value(m_dimension);
        if (jacobian.rows() == 2)
        {
          value <<
            jacobian(1, 0), -jacobian(0, 0);
        }
        else if (jacobian.rows() == 3)
        {
          value <<
            jacobian(1, 0) * jacobian(2, 1) - jacobian(2, 0) * jacobian(1, 1),
            jacobian(2, 0) * jacobian(0, 1) - jacobian(0, 0) * jacobian(2, 1),
            jacobian(0, 0) * jacobian(1, 1) - jacobian(1, 0) * jacobian(0, 1);
        }
        else
        {
          assert(false);
          value.setConstant(NAN);
        }
        return value.normalized();
      }

      inline BoundaryNormal* copy() const noexcept override
      {
        return new BoundaryNormal(*this);
      }

    private:
      const size_t m_dimension;
  };
}

#endif

