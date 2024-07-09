/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_P0_GRIDFUNCTION_H

#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Geometry/SubMesh.h"

#include "P0.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range, class Mesh>
  struct Traits<Variational::GridFunction<Variational::P0<Range, Mesh>>>
  {
    using FESType = Variational::P0<Range, Mesh>;
    using MeshType = typename FormLanguage::Traits<FESType>::MeshType;
    using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
    using ContextType = typename FormLanguage::Traits<FESType>::ContextType;
    using ElementType = typename FormLanguage::Traits<FESType>::ElementType;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup GridFunctionSpecializations
   * @brief P0 GridFunction
   */
  template <class Number, class Mesh>
  class GridFunction<P0<Number, Mesh>> final
    : public GridFunctionBase<P0<Number, Mesh>, GridFunction<P0<Number, Mesh>>>
  {
    public:
      using FESType = P0<Number, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using MeshType = typename FormLanguage::Traits<FESType>::MeshType;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using ElementType = typename FormLanguage::Traits<FESType>::ElementType;

      using Parent = GridFunctionBase<FESType, GridFunction<FESType>>;

      using Parent::getValue;
      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;

      /**
       * @brief Constructs a grid function on a finite element space.
       * @param[in] fes Finite element space to which the function belongs
       * to.
       */
      GridFunction(const FESType& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copies the grid function.
       * @param[in] other Other grid function to copy.
       */
      GridFunction(const GridFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructs the grid function.
       * @param[in] other Other grid function to move.
       */
      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Move assignment operator.
       */
      inline
      constexpr
      GridFunction& operator=(GridFunction&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      GridFunction& operator=(const GridFunction&)  = delete;

      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Real>);
        const auto& fes = this->getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        const auto& polytope = p.getPolytope();
        assert(fesMesh == polytope.getMesh());
        const size_t d = fesMesh.getDimension();
        assert(d == polytope.getDimension());
        const Index i = polytope.getIndex();
        assert(fes.getFiniteElement(d, i).getCount() == 1);
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          res = getValue({ d, i }, 0);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          assert(false);
        }
        else
        {
          assert(false);
        }
      }

      GridFunction& setWeights()
      {
        auto& data = this->getData();
        if (!this->getWeights().has_value())
        {
          auto& weights = this->getWeights().emplace(this->getFiniteElementSpace().getSize());
          if constexpr (std::is_same_v<RangeType, Real>)
          {
            assert(data.rows() == 1);
            std::copy(data.data(), data.data() + data.size(), weights.data());
          }
          else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
          {
            const auto& fes = this->getFiniteElementSpace();
            const size_t vdim = fes.getVectorDimension();
            for (size_t i = 0; i < fes.getSize(); i++)
              weights.coeffRef(i) = data.col(i).coeff(i % vdim);
          }
          else
          {
            assert(false);
            weights.setConstant(NAN);
          }
        }
        return *this;
      }

      template <class Vector>
      GridFunction& setWeights(Vector&& weights)
      {
        assert(weights.size() >= 0);
        assert(static_cast<size_t>(weights.size()) == this->getFiniteElementSpace().getSize());
        auto& data = this->getData();
        const auto& w = this->getWeights().emplace(std::forward<Vector>(weights));
        if constexpr (std::is_same_v<RangeType, Real>)
        {
          assert(data.rows() == 1);
          std::copy(w.data(), w.data() + w.size(), data.data());
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
        {
          const size_t sz = w.size();
          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          const size_t count = mesh.getVertexCount();
          const size_t vdim = fes.getVectorDimension();
          data.setZero();
          for (size_t i = 0; i < sz; i++)
            for (size_t d = 0; d < vdim; d++)
              data.col(i).coeffRef(d) = w(i % count + d * count);
        }
        else
        {
          assert(false);
          data.setConstant(NAN);
        }
        return *this;
      }

  };

  template <class Number, class Mesh>
  GridFunction(const P0<Number, Mesh>&) -> GridFunction<P0<Number, Mesh>>;
}

#endif

