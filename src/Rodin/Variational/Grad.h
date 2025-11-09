/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRAD_H
#define RODIN_VARIATIONAL_GRAD_H

/**
 * @file
 * @brief Gradient operator for scalar and vector functions.
 *
 * This file defines the gradient operator @f$ \nabla @f$ which computes spatial
 * derivatives of functions. The gradient is fundamental to many variational
 * formulations and appears in diffusion, advection, and mechanics problems.
 *
 * # Mathematical Foundation
 *
 * ## Scalar Functions
 * For a scalar function @f$ u: \Omega \rightarrow \mathbb{R} @f$, the gradient
 * is a vector field:
 * @f[
 *   \nabla u = \left(\frac{\partial u}{\partial x_1}, \frac{\partial u}{\partial x_2}, \ldots, \frac{\partial u}{\partial x_d}\right)
 * @f]
 * where @f$ d @f$ is the spatial dimension.
 *
 * ## Vector Functions
 * For a vector function @f$ \mathbf{u}: \Omega \rightarrow \mathbb{R}^n @f$, the gradient
 * is the Jacobian matrix:
 * @f[
 *   \nabla \mathbf{u} = \begin{pmatrix}
 *     \frac{\partial u_1}{\partial x_1} & \cdots & \frac{\partial u_1}{\partial x_d} \\
 *     \vdots & \ddots & \vdots \\
 *     \frac{\partial u_n}{\partial x_1} & \cdots & \frac{\partial u_n}{\partial x_d}
 *   \end{pmatrix}
 * @f]
 *
 * ## Properties
 * - **Linearity**: @f$ \nabla(au + bv) = a\nabla u + b\nabla v @f$
 * - **Product rule**: @f$ \nabla(uv) = u\nabla v + v\nabla u @f$
 * - **Chain rule**: @f$ \nabla(f(u)) = f'(u)\nabla u @f$
 *
 * # Applications
 * - **Diffusion**: @f$ -\nabla \cdot (k\nabla u) = f @f$
 * - **Advection**: @f$ \mathbf{v} \cdot \nabla u @f$
 * - **Elasticity**: Strain tensor @f$ \varepsilon = \frac{1}{2}(\nabla u + \nabla u^T) @f$
 * - **Gradient flow**: @f$ \frac{\partial u}{\partial t} = -\nabla E(u) @f$
 *
 * # Usage Examples
 *
 * ## Poisson Equation (Laplace Operator)
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * P1 Vh(mesh);
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Weak form: (grad u, grad v) = (f, v)
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))  // Diffusion term
 *         - Integral(f * v)
 *         + DirichletBC(u, 0.0);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Advection-Diffusion
 * @code{.cpp}
 * // Advection-diffusion: -nu*Laplace(u) + v·grad(u) = f
 * double nu = 0.01;  // Diffusion coefficient
 *
 * auto velocity = VectorFunction([](auto& x) {
 *   return Math::Vector2D(1.0, 0.0);  // Flow to the right
 * });
 *
 * Problem problem(u, v);
 * problem = Integral(nu * Grad(u), Grad(v))  // Diffusion
 *         + Integral(Dot(velocity, Grad(u)) * v)  // Advection
 *         - Integral(f * v);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Linear Elasticity (Strain Tensor)
 * @code{.cpp}
 * // 2D elasticity
 * P1 Vh(mesh, 2);  // Vector displacement
 *
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Strain tensor: epsilon = (grad u + grad u^T) / 2
 * auto epsilon = [](auto& w) {
 *   return 0.5 * (Grad(w) + Transpose(Grad(w)));
 * };
 *
 * // Stress-strain relation
 * double lambda = 1.0, mu = 1.0;
 *
 * Problem problem(u, v);
 * problem = Integral(lambda * Trace(epsilon(u)) * Trace(epsilon(v))
 *                   + 2 * mu * Dot(epsilon(u), epsilon(v)));
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Computing H1 Seminorm
 * @code{.cpp}
 * P1 Vh(mesh);
 * GridFunction<P1> u(Vh);
 *
 * // |u|_H^1 = ||grad u||_L^2
 * auto grad_u = Grad(u);
 * double h1_seminorm = std::sqrt(Integral(Dot(grad_u, grad_u)).compute());
 * @endcode
 *
 * ## Anisotropic Diffusion
 * @code{.cpp}
 * // Diffusion with tensor coefficient: -div(K·grad u) = f
 * auto K = MatrixFunction([](auto& x) {
 *   // Anisotropic diffusion tensor
 *   Math::Matrix2D k;
 *   k << 2.0, 0.5,
 *        0.5, 1.0;
 *   return k;
 * });
 *
 * Problem problem(u, v);
 * problem = Integral(Dot(K * Grad(u), Grad(v)))
 *         - Integral(f * v);
 *
 * problem.solve(solver);
 * @endcode
 *
 * @see Div, Jacobian, Derivative
 * @see Integral, Problem, TrialFunction
 */

#include "ForwardDecls.h"

#include "VectorFunction.h"

namespace Rodin::FormLanguage
{
  template <class FES, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
    using FESType = FES;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };

  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   */
  template <class FES, class Data, class Derived>
  class GradBase<GridFunction<FES, Data>, Derived>
    : public VectorFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, GradBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = VectorFunctionBase<ScalarType, GradBase<OperandType, Derived>>;

      GradBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      GradBase(const GradBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      GradBase(GradBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local SpatialVectorType s_out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          this->interpolate(s_out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(s_out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(s_out, *restriction);
        }
        else
        {
          assert(false);
        }
        return s_out;
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      GradBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a GridFunction
   */
  template <class FES, class Data>
  Grad(const GridFunction<FES, Data>&) -> Grad<GridFunction<FES, Data>>;

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a ShapeFunction
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, FES, Space>&)
    -> Grad<ShapeFunction<NestedDerived, FES, Space>>;
}

#endif
