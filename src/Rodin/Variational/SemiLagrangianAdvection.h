/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SEMILAGRANGIANADVECTION_H
#define RODIN_VARIATIONAL_SEMILAGRANGIANADVECTION_H

#include <vector>
#include <functional>
#include <cmath>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeIterator.h"

#include "Rodin/QF/QuadratureFormula.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  /**
   * @brief Finite element space-agnostic semi-Lagrangian variational advection step for scalar fields.
   *
   * This class implements a semi-Lagrangian method for advecting scalar fields that works with
   * any scalar finite element space (Lagrange Pk, bubbles, DG, modal) on both volume and surface meshes.
   *
   * The method solves the advection equation:
   * @f[
   * \frac{\partial u}{\partial t} + \beta \cdot \nabla u = 0
   * @f]
   * 
   * **Algorithm Overview:**
   * 1. For each quadrature point, backtrace along characteristic curves
   * 2. Use RK2/RK3 integration with no-cross substep strategy
   * 3. Evaluate old solution at departure points
   * 4. Assemble right-hand side with proper metric factors
   * 5. Solve mass matrix system: M c^{n+1} = b
   *
   * **Key Features:**
   * - Works with any scalar FE space via template design
   * - Supports both volume meshes (metric factor = det(J)) and surface meshes (metric factor = sqrt(det(J^T J)))
   * - Conservative and non-conservative forms
   * - Handles inflow boundary conditions
   * - Element-agnostic through Rodin FE infrastructure
   *
   * @tparam FES Finite element space type
   * @tparam ScalarType Numerical scalar type (Real, Complex, etc.)
   */
  template <class FES, class ScalarType = Real>
  class SemiLagrangianAdvection
  {
    public:
      using FESType = FES;
      using MeshType = typename FormLanguage::Traits<FES>::MeshType;
      using GridFunctionType = GridFunction<FES>;
      using VectorFieldType = GridFunction<FES>;
      using VectorType = Math::Vector<ScalarType>;
      using SpatialVectorType = Math::SpatialVector<ScalarType>;
      using PointType = Geometry::Point;

      /**
       * @brief Runge-Kutta integrator type for characteristic tracing
       */
      enum class IntegratorType
      {
        RK2,  ///< Second-order Runge-Kutta
        RK3   ///< Third-order Runge-Kutta
      };

      /**
       * @brief Mass matrix type for final solve
       */
      enum class MassMatrixType
      {
        Lumped,     ///< Lumped mass matrix (diagonal)
        Consistent  ///< Consistent mass matrix (full)
      };

      /**
       * @brief Constructs a semi-Lagrangian advection solver.
       *
       * @param fes Finite element space
       * @param mesh Mesh with element connectivity
       * @param oldSolution Previous time step solution coefficients c^n
       * @param velocityField Velocity field β(x,t)
       * @param timestep Time step Δt
       */
      SemiLagrangianAdvection(
          const FES& fes,
          const MeshType& mesh,
          const GridFunctionType& oldSolution,
          const VectorFieldType& velocityField,
          ScalarType timestep);

      /**
       * @brief Sets the mass matrix to use for the final solve.
       * @param massMatrix Precomputed mass matrix M
       * @param type Type of mass matrix (lumped or consistent)
       * @return Reference to this object for method chaining
       */
      SemiLagrangianAdvection& setMassMatrix(
          const Math::Matrix<ScalarType>& massMatrix,
          MassMatrixType type = MassMatrixType::Consistent);

      /**
       * @brief Sets the Runge-Kutta integrator type for characteristic tracing.
       * @param type Integrator type (RK2 or RK3)
       * @return Reference to this object for method chaining
       */
      SemiLagrangianAdvection& setIntegratorType(IntegratorType type);

      /**
       * @brief Enables conservative form with divergence compensation.
       * @param enable Whether to enable conservative form
       * @return Reference to this object for method chaining
       *
       * When enabled, multiplies by exp(-∫ div βτ dt) using midpoint rule along substep.
       */
      SemiLagrangianAdvection& setConservativeForm(bool enable = true);

      /**
       * @brief Sets inflow boundary data for boundary conditions.
       * @param boundaryData Function g(x,t) for inflow boundaries
       * @return Reference to this object for method chaining
       */
      template <class BoundaryFunction>
      SemiLagrangianAdvection& setInflowBoundaryData(const BoundaryFunction& boundaryData);

      /**
       * @brief Sets the exit time factor for no-cross substep strategy.
       * @param theta Safety factor θ ∈ (0,1) for limiting substep size
       * @return Reference to this object for method chaining
       *
       * Substep size is limited by: δt = min(Δt_rem, θ * t_exit)
       */
      SemiLagrangianAdvection& setExitTimeFactor(ScalarType theta);

      /**
       * @brief Performs one semi-Lagrangian advection step.
       * @param[out] newSolution New solution coefficients c^{n+1}
       *
       * This method implements the complete semi-Lagrangian algorithm:
       * 1. Assembles RHS vector b by backtracing and evaluating old solution
       * 2. Solves mass matrix system M c^{n+1} = b
       */
      void step(GridFunctionType& newSolution);

    private:
      const FES& m_fes;
      const MeshType& m_mesh;
      const GridFunctionType& m_oldSolution;
      const VectorFieldType& m_velocityField;
      ScalarType m_timestep;

      // Algorithm parameters
      IntegratorType m_integratorType = IntegratorType::RK2;
      MassMatrixType m_massMatrixType = MassMatrixType::Consistent;
      bool m_conservativeForm = false;
      ScalarType m_exitTimeFactor = ScalarType(0.9);

      // Mass matrix
      std::optional<std::reference_wrapper<const Math::Matrix<ScalarType>>> m_massMatrix;

      // Boundary data
      std::function<ScalarType(const PointType&, ScalarType)> m_inflowBoundaryData;

      /**
       * @brief Backtraces from physical point along characteristic curve.
       * @param startPoint Starting physical point x_q
       * @param time Current time t
       * @return Departure point x* after backtracing over timestep
       */
      PointType backtrace(const PointType& startPoint, ScalarType time);

      /**
       * @brief Performs one RK2 substep for characteristic tracing.
       * @param point Current point in reference coordinates
       * @param element Current element index
       * @param dt Substep size
       * @param time Current time
       * @return New point and element after substep
       */
      std::pair<SpatialVectorType, size_t> rk2Substep(
          const SpatialVectorType& point,
          size_t element,
          ScalarType dt,
          ScalarType time);

      /**
       * @brief Performs one RK3 substep for characteristic tracing.
       * @param point Current point in reference coordinates
       * @param element Current element index
       * @param dt Substep size
       * @param time Current time
       * @return New point and element after substep
       */
      std::pair<SpatialVectorType, size_t> rk3Substep(
          const SpatialVectorType& point,
          size_t element,
          ScalarType dt,
          ScalarType time);

      /**
       * @brief Computes exit time from current reference element.
       * @param refPoint Point in reference coordinates
       * @param velocity Velocity in reference coordinates
       * @param element Current element
       * @return Exit time and exit face index
       */
      std::pair<ScalarType, size_t> computeExitTime(
          const SpatialVectorType& refPoint,
          const SpatialVectorType& velocity,
          size_t element);

      /**
       * @brief Hops to neighboring element through face connectivity.
       * @param refPoint Point on exit face in reference coordinates
       * @param currentElement Current element index
       * @param exitFace Exit face index
       * @return New element index and transformed reference point
       */
      std::pair<size_t, SpatialVectorType> hopToNeighbor(
          const SpatialVectorType& refPoint,
          size_t currentElement,
          size_t exitFace);

      /**
       * @brief Evaluates old solution at departure point.
       * @param departurePoint Physical coordinates of departure point
       * @return Value of u_h^n(x*)
       */
      ScalarType evaluateOldSolution(const PointType& departurePoint);

      /**
       * @brief Computes metric factor for assembly.
       * @param element Element index
       * @param refPoint Reference coordinates
       * @return Metric factor (det(J) for volume, sqrt(det(J^T J)) for surface)
       */
      ScalarType computeMetricFactor(size_t element, const SpatialVectorType& refPoint);

      /**
       * @brief Projects velocity to tangent space for surface meshes.
       * @param velocity Original velocity vector
       * @param normal Surface normal vector
       * @return Tangential velocity βτ = (I - n⊗n)β
       */
      SpatialVectorType projectToTangentSpace(
          const SpatialVectorType& velocity,
          const SpatialVectorType& normal);

      /**
       * @brief Computes divergence compensation factor for conservative form.
       * @param startPoint Starting point of characteristic
       * @param endPoint Ending point of characteristic
       * @param time Current time
       * @return exp(-∫ div βτ dt)
       */
      ScalarType computeConservativeFactor(
          const PointType& startPoint,
          const PointType& endPoint,
          ScalarType time);

      /**
       * @brief Checks if point hits inflow boundary during backtracing.
       * @param point Point to check
       * @return True if point is on inflow boundary
       */
      bool isInflowBoundary(const PointType& point);

      /**
       * @brief Assembles the right-hand side vector b.
       * @param[out] rhs Right-hand side vector to fill
       */
      void assembleRHS(VectorType& rhs);

      /**
       * @brief Solves the mass matrix system M c^{n+1} = b.
       * @param rhs Right-hand side vector b
       * @param[out] solution Solution vector c^{n+1}
       */
      void solveMassMatrixSystem(const VectorType& rhs, VectorType& solution);
  };
}

#include "SemiLagrangianAdvection.hpp"

#endif