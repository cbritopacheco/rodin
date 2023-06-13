/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_P1_GRIDFUNCTION_H

#include "Rodin/Variational/GridFunction.h"
#include "P1.h"

namespace Rodin::Variational
{
  /**
   * @ingroup GridFunctionSpecializations
   * @brief P1 GridFunction
   */
  template <class ... Ts>
  class GridFunction<P1<Ts...>> final : public GridFunctionBase<GridFunction<P1<Ts...>>, P1<Ts...>>
  {
    public:
      /// Type of finite element space to which the GridFunction belongs to
      using FES = P1<Ts...>;

      /// Type of finite element
      using Element = typename FES::Element;

      /// Parent class
      using Parent = GridFunctionBase<GridFunction<P1<Ts...>>, P1<Ts...>>;

      /// Range type of value
      using RangeType = typename FES::RangeType;

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
      GridFunction(const FES& fes)
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

      inline
      auto getValue(const Geometry::Point& p) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const auto& r = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          Scalar res = 0;
          for (Index local = 0; local < fe.getCount(); local++)
          {
            const auto& basis = fe.getBasis(local);
            res += getValue({d, i}, local) * basis(r);
          }
          return res;
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          const size_t vdim = fes.getVectorDimension();
          const size_t dofs = fe.getCount();
          Math::Vector res = Math::Vector::Zero(vdim);
          for (Index local = 0; local < dofs; local++)
          {
            const auto& basis = fe.getBasis(local);
            res += getValue({d, i}, local).coeff(local % vdim) * basis(r);
          }
          return res;
        }
        else
        {
          assert(false);
          return void();
        }
      }
  };

  template <class ... Ts>
  GridFunction(const P1<Ts...>&) -> GridFunction<P1<Ts...>>;

}
#endif
