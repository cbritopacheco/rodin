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
    : public GridFunctionBase<P1<Range, Mesh>, GridFunction<P1<Range, Mesh>>>
  {
    public:
      /// Type of finite element space to which the GridFunction belongs to
      using FESType = P1<Range, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using MeshType = typename FormLanguage::Traits<FESType>::MeshType;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ContextType = typename FormLanguage::Traits<FESType>::ContextType;

      using ElementType = typename FormLanguage::Traits<FESType>::ElementType;

      /// Parent class
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
        const auto& fes = this->getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        const auto& polytope = p.getPolytope();
        assert(fesMesh == polytope.getMesh());
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& r = p.getReferenceCoordinates();
        if constexpr (std::is_same_v<RangeType, Real>)
        {
          res = Real(0);
          for (Index local = 0; local < fe.getCount(); local++)
            res += getValue({d, i}, local) * fe.getBasis(local)(r);
        }
        else if constexpr (std::is_same_v<RangeType, Complex>)
        {
          res = Complex(0, 0);
          for (Index local = 0; local < fe.getCount(); local++)
            res += 0.5 * getValue({d, i}, local) * Complex(1, 1) * fe.getBasis(local)(r);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
        {
          const size_t vdim = fes.getVectorDimension();
          res.resize(vdim);
          res.setZero();
          RangeType basis;
          for (Index local = 0; local < fe.getCount(); local++)
          {
            fe.getBasis(local)(basis, r);
            res += getValue({d, i}, local).coeff(local % vdim) * basis;
          }
        }
        else
        {
          assert(false);
        }
      }

      GridFunction& setWeights()
      {
        auto& data = this->getData();
        auto& w = this->getWeights().emplace(this->getFiniteElementSpace().getSize());
        if constexpr (std::is_same_v<RangeType, Real>)
        {
          assert(data.rows() == 1);
          w = data.transpose();
        }
        else if constexpr (std::is_same_v<RangeType, Complex>)
        {
          assert(data.rows() == 1);
          w = data.transpose() / Complex(1, 1);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          const auto& fes = this->getFiniteElementSpace();
          const size_t vdim = fes.getVectorDimension();
          for (size_t i = 0; i < fes.getSize(); i++)
            w.coeffRef(i) = data.col(i).coeff(i % vdim);
        }
        else
        {
          assert(false);
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
          data = w.transpose();
        }
        else if constexpr (std::is_same_v<RangeType, Complex>)
        {
          assert(data.rows() == 1);
          data = w.adjoint() * Complex(1, -1);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          const size_t sz = w.size();
          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          const size_t count = mesh.getVertexCount();
          const size_t vdim = fes.getVectorDimension();
          for (size_t i = 0; i < sz; i++)
            for (size_t d = 0; d < vdim; d++)
              data.col(i).coeffRef(d) = w(i % count + d * count);
        }
        else
        {
          assert(false);
        }
        return *this;
      }

  };

  template <class Range, class Mesh>
  GridFunction(const P1<Range, Mesh>&) -> GridFunction<P1<Range, Mesh>>;
}

#endif
