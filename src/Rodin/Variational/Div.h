/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

/**
 * @file
 * @brief Divergence operator for vector-valued functions.
 *
 * This file defines the divergence operator @f$ \nabla \cdot @f$ which maps
 * vector fields to scalar fields. The divergence measures the "outflow" or
 * "expansion" of a vector field at each point.
 *
 * # Mathematical Foundation
 *
 * For a vector field @f$ \mathbf{u} = (u_1, u_2, \ldots, u_d): \Omega \rightarrow \mathbb{R}^d @f$,
 * the divergence is the scalar field:
 * @f[
 *   \nabla \cdot \mathbf{u} = \frac{\partial u_1}{\partial x_1} +
 *                              \frac{\partial u_2}{\partial x_2} + \cdots +
 *                              \frac{\partial u_d}{\partial x_d} = \mathrm{tr}(\nabla \mathbf{u})
 * @f]
 *
 * The divergence is the trace of the Jacobian matrix.
 *
 * ## Properties
 * - **Linearity**: @f$ \nabla \cdot (a\mathbf{u} + b\mathbf{v}) = a\nabla \cdot \mathbf{u} + b\nabla \cdot \mathbf{v} @f$
 * - **Product rule**: @f$ \nabla \cdot (f\mathbf{u}) = f\nabla \cdot \mathbf{u} + \nabla f \cdot \mathbf{u} @f$
 * - **Divergence theorem**: @f$ \int_\Omega \nabla \cdot \mathbf{u} \, dx = \int_{\partial\Omega} \mathbf{u} \cdot \mathbf{n} \, ds @f$
 *
 * # Applications
 * - **Incompressibility constraint**: @f$ \nabla \cdot \mathbf{u} = 0 @f$ (Stokes, Navier-Stokes)
 * - **Conservation laws**: Mass, momentum, energy conservation
 * - **Elasticity**: Volumetric strain @f$ \varepsilon_v = \nabla \cdot \mathbf{u} @f$
 * - **Maxwell equations**: @f$ \nabla \cdot \mathbf{E} = \rho / \varepsilon_0 @f$
 *
 * # Usage Examples
 *
 * ## Incompressible Flow (Stokes/Navier-Stokes)
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * // Velocity space (vector) and pressure space (scalar)
 * P1 Vh(mesh, 2);  // Velocity
 * P0 Qh(mesh);     // Pressure
 *
 * TrialFunction u(Vh), p(Qh);  // Velocity and pressure
 * TestFunction  v(Vh), q(Qh);
 *
 * // Incompressibility constraint: div(u) = 0
 * Problem problem(u, v, p, q);
 * problem = Integral(Grad(u), Grad(v))      // Viscous term
 *         - Integral(p * Div(v))             // Pressure gradient
 *         - Integral(q * Div(u));            // Incompressibility
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Mixed Poisson (Darcy Flow)
 * @code{.cpp}
 * // Darcy's law: u = -K·grad p, div(u) = f
 * P0 Vh(mesh, 2);  // Flux (vector)
 * P0 Qh(mesh);     // Pressure (scalar)
 *
 * TrialFunction u(Vh), p(Qh);
 * TestFunction  v(Vh), q(Qh);
 *
 * // Permeability
 * double K = 1.0;
 *
 * Problem problem(u, v, p, q);
 * problem = Integral((1.0/K) * Dot(u, v))    // Constitutive law
 *         + Integral(Div(v) * p)              // Flux-pressure coupling
 *         + Integral(q * Div(u))              // Mass conservation
 *         - Integral(q * f);                  // Source term
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Volumetric Strain in Elasticity
 * @code{.cpp}
 * // Linear elasticity with volumetric/deviatoric decomposition
 * P1 Vh(mesh, 3);  // 3D displacement
 *
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Volumetric strain: epsilon_v = div(u)
 * double lambda = 1.0, mu = 1.0;
 *
 * Problem problem(u, v);
 * problem = Integral(lambda * Div(u) * Div(v))  // Volumetric energy
 *         + Integral(/* deviatoric terms */);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Verifying Divergence-Free Solution
 * @code{.cpp}
 * // Check incompressibility of computed velocity field
 * P1 Vh(mesh, 2);
 * GridFunction<P1> u(Vh);  // Velocity solution
 *
 * // Compute L2 norm of div(u)
 * auto div_u = Div(u);
 * double div_norm = std::sqrt(Integral(div_u * div_u).compute());
 *
 * std::cout << "||div(u)||_L2 = " << div_norm << std::endl;
 * // Should be near zero for incompressible flow
 * @endcode
 *
 * ## Divergence Theorem Verification
 * @code{.cpp}
 * // Verify: int(div u) dx = int(u·n) ds
 * P1 Vh(mesh, 2);
 * GridFunction<P1> u(Vh);
 *
 * // Volume integral of divergence
 * double volume_integral = Integral(Div(u)).compute();
 *
 * // Boundary integral of normal flux
 * auto n = BoundaryNormal();
 * double boundary_integral = BoundaryIntegral(Dot(u, n)).compute();
 *
 * std::cout << "Volume: " << volume_integral << std::endl;
 * std::cout << "Boundary: " << boundary_integral << std::endl;
 * // Should be approximately equal
 * @endcode
 *
 * @see Grad, Trace, Jacobian
 * @see Integral, BoundaryIntegral, Problem
 */

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DivSpecializations Div Template Specializations
   * @brief Template specializations of the Div class.
   * @see Div
   */

  /**
   * @ingroup RodinVariational
   * @brief Base class for divergence operators.
   *
   * @tparam Operand Type of the operand (vector-valued function)
   * @tparam Derived Derived class type (CRTP pattern)
   */
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a GridFunction.
   *
   * Computes the divergence @f$ \nabla \cdot \mathbf{u} @f$ of a vector-valued
   * grid function @f$ \mathbf{u} @f$. The result is a scalar-valued function.
   *
   * ## Mathematical Definition
   * For a vector field @f$ \mathbf{u}: \Omega \rightarrow \mathbb{R}^d @f$:
   * @f[
   *   \text{div}(\mathbf{u}) = \nabla \cdot \mathbf{u} =
   *   \sum_{i=1}^d \frac{\partial u_i}{\partial x_i}
   * @f]
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @tparam Derived Derived class type
   */
  template <class FES, class Data, class Derived>
  class DivBase<GridFunction<FES, Data>, Derived>
    : public ScalarFunctionBase<typename FormLanguage::Traits<FES>::ScalarType, DivBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = GridFunction<FES, Data>;

      /// Parent class
      using Parent = ScalarFunctionBase<ScalarType, DivBase<OperandType, Derived>>;

      /**
       * @brief Constructs the Div of a @f$ \mathbb{P}_1 @f$ function @f$ u
       * @f$.
       * @param[in] u P1 GridFunction
       */
      DivBase(const OperandType& u)
        : m_u(u)
      {}

      /**
       * @brief Copy constructor
       */
      DivBase(const DivBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      DivBase(DivBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      ScalarType getValue(const Geometry::Point& p) const
      {
        static thread_local ScalarType s_out;
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

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      DivBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
