/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_WORKLOADS_COUPLEDREACTIONDIFFUSIONWORKLOAD_H
#define RODIN_WORKLOADS_COUPLEDREACTIONDIFFUSIONWORKLOAD_H

#include <benchmark/benchmark.h>
#include "convergence/ReactionDiffusion.h"
#include "QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Shared two-field reaction--diffusion forms and manufactured data.
   *
   * On the unit box, @f$-\kappa_i\Delta u_i+u_i+\alpha u_{1-i}=f_i@f$
   * uses @f$\kappa=(1,2)@f$ and @f$\alpha=0.2@f$. The reaction matrix
   * has eigenvalues 0.8 and 1.2, so full Dirichlet elimination preserves
   * positive definiteness. The affine/quadratic fields and analytic sources
   * are supplied by @ref Convergence::ReactionDiffusionData.
   * Complex fields have distinct constant component multipliers; their sources
   * are shifted analytically, not obtained from an assembled matrix.
   * No assembly backend or solver policy is selected here.
   *
   * @tparam K Degree of both scalar H1 spaces.
   * @tparam Scalar Real or complex coefficient field.
   */
  template <size_t K, class Scalar>
  class CoupledReactionDiffusionWorkload
  {
    public:
    /// @brief Matches the polynomial patch quadrature in coupled certification.
    static constexpr size_t QuadratureOrder = 16;

    /**
     * @brief Returns the fixed component multiplier for the coefficient field.
     * @param component Zero-based component index, smaller than two.
     * @returns Real unit multiplier or a distinct complex component multiplier.
     */
    static Scalar getMultiplier(size_t component)
    {
      assert(component < 2);
      if constexpr (std::is_same_v<Scalar, Complex>)
        return component == 0 ? Scalar(1, 1) : Scalar(2, -1);
      else
        return 1;
    }

    /**
     * @brief Returns the manufactured physical field for one component.
     * @param data Certified real manufactured data.
     * @param component Zero-based component index, smaller than two.
     * @returns Function taking values in the requested coefficient field.
     */
    static auto getSolution(
      const Convergence::ReactionDiffusionData& data, size_t component)
    {
      const auto exact = data.getSolution(component);
      const Scalar multiplier = getMultiplier(component);
      const auto evaluate = [exact, multiplier](const Geometry::Point& point) {
        return multiplier * exact(point);
      };
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Variational::ComplexFunction(evaluate);
      else
        return Variational::RealFunction(evaluate);
    }

    /**
     * @brief Returns the analytically scaled coupled source.
     * @param data Certified real manufactured data.
     * @param component Zero-based component index, smaller than two.
     * @returns Source for the original nonzero coupling coefficient.
     */
    static auto getSource(
      const Convergence::ReactionDiffusionData& data, size_t component)
    {
      const auto source = data.getSource(component);
      const auto other = data.getSolution(1 - component);
      const Scalar multiplier = getMultiplier(component),
                   otherMultiplier = getMultiplier(1 - component);
      const auto evaluate = [source, other, multiplier, otherMultiplier](
                              const Geometry::Point& point) {
        return multiplier * source(point) +
          Convergence::ReactionDiffusionData::Coupling * (otherMultiplier - multiplier) *
          other(point);
      };
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Variational::ComplexFunction(evaluate);
      else
        return Variational::RealFunction(evaluate);
    }

    /**
     * @brief Defines the physical block problem or its omitted-coupling control.
     *
     * Original reference integrals replace only local operator kernels.
     * Omitting coupling leaves analytic sources and prescribed traces unchanged.
     * Both exact fields are interpolated through their own DOF functionals.
     *
     * @tparam Reference Whether original local operator kernels are requested.
     * @tparam TestWeighted Whether the real second diffusion coefficient is
     * placed on the test operand to exercise generic gradient expressions.
     * @tparam U First trial function type.
     * @tparam W Second trial function type.
     * @tparam V First test function type.
     * @tparam Z Second test function type.
     * @tparam ProblemType Two-field problem representation.
     * @param u First trial function.
     * @param w Second trial function.
     * @param v First test function.
     * @param z Second test function.
     * @param problem Destination formulation.
     * @param data Manufactured real fields and physical sources.
     * @param omitCoupling Whether both off-diagonal operator terms are zeroed.
     */
    template <bool Reference = false, bool TestWeighted = false, class U, class W,
      class V, class Z, class ProblemType>
    static void configure(U& u, W& w, V& v, Z& z, ProblemType& problem,
      const Convergence::ReactionDiffusionData& data, bool omitCoupling = false)
    {
      using namespace Variational;
      const auto exactU = getSolution(data, 0), exactW = getSolution(data, 1);
      const Real coupling =
        omitCoupling ? Real(0) : Convergence::ReactionDiffusionData::Coupling;
      auto aUU = Integral(Grad(u), Grad(v));
      auto aWW = [&]() {
        if constexpr (TestWeighted)
          return Integral(Grad(w), Real(2) * Grad(z));
        else
          return Integral(Real(2) * Grad(w), Grad(z));
      }();
      auto rUU = Integral(u, v), rWW = Integral(w, z);
      auto rUW = Integral(coupling * w, v), rWU = Integral(coupling * u, z);
      auto bU = Integral(getSource(data, 0), v), bW = Integral(getSource(data, 1), z);
      aUU.setOrder(QuadratureOrder);
      aWW.setOrder(QuadratureOrder);
      rUU.setOrder(QuadratureOrder);
      rWW.setOrder(QuadratureOrder);
      rUW.setOrder(QuadratureOrder);
      rWU.setOrder(QuadratureOrder);
      bU.setOrder(QuadratureOrder);
      bW.setOrder(QuadratureOrder);
      if constexpr (Reference)
      {
        problem = ReferenceIntegral(aUU) + ReferenceIntegral(rUU) +
          ReferenceIntegral(rUW) - bU + ReferenceIntegral(aWW) + ReferenceIntegral(rWW) +
          ReferenceIntegral(rWU) - bW + DirichletBC(u, exactU) + DirichletBC(w, exactW);
      }
      else
      {
        problem = aUU + rUU + rUW - bU + aWW + rWW + rWU - bW + DirichletBC(u, exactU) +
          DirichletBC(w, exactW);
      }
      u.getSolution() = exactU;
      w.getSolution() = exactW;
    }

    /**
     * @brief Registers both physical fields on seven geometries at three sizes.
     * @tparam Callback Callable accepting state and reference geometry.
     * @param scope Owning operation, such as Assembly or CG.
     * @param callback Operation to measure.
     */
    template <class Callback>
    static void registerCases(const char* scope, Callback callback)
    {
      using Geometry::Polytope;
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      {
        const std::string name = std::string("CoupledReactionDiffusion/") + scope +
          "/H1/P" + std::to_string(K) + "/" +
          (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(
          name.c_str(), [callback, geometry](auto& state) { callback(state, geometry); })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->UseRealTime();
      }
    }
  };
}
#endif
