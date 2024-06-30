/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_P1_GRIDFUNCTION_H

#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Geometry/SubMesh.h"

#include "P1.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range, class Mesh>
  struct Traits<Variational::GridFunction<Variational::P1<Range, Mesh>>>
  {
    using FESType = Variational::P1<Range, Mesh>;
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
   * @brief P1 GridFunction
   */
  template <class Range, class Mesh>
  class GridFunction<P1<Range, Mesh>> final
    : public GridFunctionBase<GridFunction<P1<Range, Mesh>>>
  {
    public:
      /// Type of finite element space to which the GridFunction belongs to
      using FESType = P1<Range, Mesh>;

      using MeshType = typename FormLanguage::Traits<FESType>::MeshType;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using ElementType = typename FormLanguage::Traits<FESType>::ElementType;

      /// Parent class
      using Parent = GridFunctionBase<GridFunction<FESType>>;

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

      Scalar interpolate(const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Scalar>);
        const auto& fes = this->getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        const auto& polytope = p.getPolytope();
        assert(fesMesh == polytope.getMesh());
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& r = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        Scalar res = 0;
        for (Index local = 0; local < fe.getCount(); local++)
        {
          const auto& basis = fe.getBasis(local);
          res += getValue({d, i}, local) * basis(r);
        }
        return res;
      }

      void interpolate(Math::Vector<Scalar>& res, const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<Scalar>>);
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& r = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        const size_t vdim = fes.getVectorDimension();
        const size_t dofs = fe.getCount();
        res.resize(vdim);
        res.setZero();
        for (Index local = 0; local < dofs; local++)
        {
          const auto& basis = fe.getBasis(local);
          res += getValue({d, i}, local).coeff(local % vdim) * basis(r);
        }
      }

      GridFunction& setWeights()
      {
        auto& data = this->getData();
        auto& weights = this->getWeights().emplace(this->getFiniteElementSpace().getSize());
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(data.rows() == 1);
          std::copy(data.data(), data.data() + data.size(), weights.data());
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Scalar>>)
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
        return *this;
      }

      template <class Vector>
      GridFunction& setWeights(Vector&& weights)
      {
        assert(weights.size() >= 0);
        assert(static_cast<size_t>(weights.size()) == this->getFiniteElementSpace().getSize());
        auto& data = this->getData();
        const auto& w = this->getWeights().emplace(std::forward<Vector>(weights));
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(data.rows() == 1);
          std::copy(w.data(), w.data() + w.size(), data.data());
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Scalar>>)
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

  template <class Range, class Mesh>
  GridFunction(const P1<Range, Mesh>&) -> GridFunction<P1<Range, Mesh>>;
}

#endif
