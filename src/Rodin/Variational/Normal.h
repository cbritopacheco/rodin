#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

#include "ForwardDecls.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Outward unit normal.
   */
  class Normal final : public VectorFunctionBase<Normal>
  {
    public:
      using Parent = VectorFunctionBase<Normal>;

      /**
       * @brief Constructs the outward unit normal.
       */
      Normal(const Geometry::MeshBase& surface)
        : m_dimension(surface.getSpaceDimension())
      {
        assert(m_dimension > 0);
      }

      Normal(const Normal& other)
        : Parent(other),
          m_dimension(other.m_dimension)
      {}

      Normal(Normal&& other)
        : Parent(std::move(other)),
          m_dimension(std::move(other.m_dimension))
      {}

      inline
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      Math::SpatialVector<Real> getValue(const Geometry::Point& p) const
      {
        Math::SpatialVector<Real> res;
        // this->getValue(res, p);
        return res;
      }

      inline Normal* copy() const noexcept override
      {
        return new Normal(*this);
      }

    private:
      const size_t m_dimension;
  };
}

#endif
