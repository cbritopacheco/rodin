#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "FiniteElementSpace.h"
#include "GridFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Function representing the Derivative of a GridFunction in H1 space
   *
   * Given a GridFunction @f$ u : \mathbb{R}^s \rightarrow \mathbb{R}^d @f$,
   * this class represents the derivative in the @f$ i @f$-th direction of the
   * @f$ j @f$-th component:
   * @f[
   *   \dfrac{\partial u_j}{ \partial x_i } ,
   * @f]
   * where @f$ \quad 0 \leq i < s @f$ and @f$ \ 0 \leq j < d @f$.
   */
  // template <class ... Ts>
  // class Derivative<GridFunction<H1<Ts...>>> final
  //   : public ScalarFunctionBase<Derivative<GridFunction<H1<Ts...>>>>
  // {
  //   public:
  //     using Operand = GridFunction<H1<Ts...>>;
  //     using Parent = ScalarFunctionBase<Derivative<Operand>>;

  //     /**
  //      * @brief Constructs the derivative in the i-th direction of the j-th
  //      * component
  //      * @param[in] direction Spatial direction @f$ x_i @f$
  //      * @param[in] component Component @f$ u_j @f$ to differentiate
  //      * @param[in] u GridFunction in H1 space
  //      */
  //     Derivative(size_t direction, size_t component, const Operand& u)
  //       : m_direction(direction),
  //         m_component(component),
  //         m_u(u)
  //     {}

  //     Derivative(const Derivative& other)
  //       : Parent(other),
  //         m_direction(other.m_direction),
  //         m_component(other.m_component),
  //         m_u(other.m_u)
  //     {}

  //     Derivative(Derivative&& other)
  //       : Parent(std::move(other)),
  //         m_direction(std::move(other.m_direction)),
  //         m_component(std::move(other.m_component)),
  //         m_u(std::move(other.m_u))
  //     {}

  //     inline
  //     constexpr
  //     Scalar getValue(const Geometry::Point&) const
  //     {
  //       assert(false);
  //       return Scalar(NAN);
  //     }

  //   private:
  //     const size_t m_direction;
  //     const size_t m_component;
  //     std::reference_wrapper<const Operand> m_u;
  // };

  /**
   * @brief %Utility function for computing @f$ \partial_x u @f$
   * @param[in] u GridFunction in H1 space
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ x @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial x}
   * @f$
   */
  template <class FES>
  auto Dx(GridFunction<FES>& u)
  {
    assert(u.getRangeType() == RangeType::Scalar);
    return Derivative(0, 0, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_y u @f$
   * @param[in] u GridFunction in H1 space
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial y}
   * @f$
   */
  template <class FES>
  auto Dy(GridFunction<FES>& u)
  {
    assert(u.getRangeType() == RangeType::Scalar);
    return Derivative(1, 0, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_z u @f$
   * @param[in] u GridFunction in H1 space
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial z}
   * @f$
   */
  template <class FES>
  auto Dz(GridFunction<FES>& u)
  {
    assert(u.getRangeType() == RangeType::Scalar);
    return Derivative(2, 0, u);
  }
}

#endif
