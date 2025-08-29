/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SEMILAGRANGIANADVECTION_HPP
#define RODIN_VARIATIONAL_SEMILAGRANGIANADVECTION_HPP

#include "SemiLagrangianAdvection.h"

#include "Rodin/Assembly/Sequential.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/LinearAlgebra.h"

namespace Rodin::Variational
{
  template <class FES, class ScalarType>
  SemiLagrangianAdvection<FES, ScalarType>::SemiLagrangianAdvection(
      const FES& fes,
      const MeshType& mesh,
      const GridFunctionType& oldSolution,
      const VectorFieldType& velocityField,
      ScalarType timestep)
    : m_fes(fes),
      m_mesh(mesh),
      m_oldSolution(oldSolution),
      m_velocityField(velocityField),
      m_timestep(timestep)
  {
    // Validate inputs
    assert(timestep > ScalarType(0));
    assert(&fes.getMesh() == &mesh);
    assert(&oldSolution.getFiniteElementSpace().getMesh() == &mesh);
    assert(&velocityField.getFiniteElementSpace().getMesh() == &mesh);
  }

  template <class FES, class ScalarType>
  SemiLagrangianAdvection<FES, ScalarType>& 
  SemiLagrangianAdvection<FES, ScalarType>::setMassMatrix(
      const Math::Matrix<ScalarType>& massMatrix,
      MassMatrixType type)
  {
    m_massMatrix = std::cref(massMatrix);
    m_massMatrixType = type;
    return *this;
  }

  template <class FES, class ScalarType>
  SemiLagrangianAdvection<FES, ScalarType>& 
  SemiLagrangianAdvection<FES, ScalarType>::setIntegratorType(IntegratorType type)
  {
    m_integratorType = type;
    return *this;
  }

  template <class FES, class ScalarType>
  SemiLagrangianAdvection<FES, ScalarType>& 
  SemiLagrangianAdvection<FES, ScalarType>::setConservativeForm(bool enable)
  {
    m_conservativeForm = enable;
    return *this;
  }

  template <class FES, class ScalarType>
  template <class BoundaryFunction>
  SemiLagrangianAdvection<FES, ScalarType>& 
  SemiLagrangianAdvection<FES, ScalarType>::setInflowBoundaryData(const BoundaryFunction& boundaryData)
  {
    m_inflowBoundaryData = [boundaryData](const PointType& p, ScalarType t) -> ScalarType {
      return boundaryData(p, t);
    };
    return *this;
  }

  template <class FES, class ScalarType>
  SemiLagrangianAdvection<FES, ScalarType>& 
  SemiLagrangianAdvection<FES, ScalarType>::setExitTimeFactor(ScalarType theta)
  {
    assert(theta > ScalarType(0) && theta < ScalarType(1));
    m_exitTimeFactor = theta;
    return *this;
  }

  template <class FES, class ScalarType>
  void SemiLagrangianAdvection<FES, ScalarType>::step(GridFunctionType& newSolution)
  {
    // Get the size of the finite element space
    const size_t dofCount = m_fes.getCount();
    
    // Initialize RHS vector
    VectorType rhs(dofCount);
    rhs.setZero();

    // Assemble the right-hand side
    assembleRHS(rhs);

    // Solve the mass matrix system
    VectorType solutionVector(dofCount);
    solveMassMatrixSystem(rhs, solutionVector);

    // Set the new solution
    newSolution.setLeaf(std::move(solutionVector));
  }

  template <class FES, class ScalarType>
  void SemiLagrangianAdvection<FES, ScalarType>::assembleRHS(VectorType& rhs)
  {
    // Get mesh dimension for iteration
    const size_t meshDim = m_mesh.getMeshDimension();
    
    // Iterate over all elements of maximum dimension
    for (auto it = m_mesh.getPolytopeIterator(meshDim); it; ++it)
    {
      const auto& element = *it;
      const size_t elementIndex = element.getIndex();
      
      // Get element DOFs
      const auto& elementDOFs = m_fes.getDOFs(meshDim, elementIndex);
      
      // Get quadrature rule for this element
      const auto qf = QF::getQuadratureRule(element.getGeometry(), 2 * m_fes.getPolynomialDegree());
      
      // Loop over quadrature points
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        // Get quadrature point and weight in reference coordinates
        const auto refPoint = qf.getPoint(q);
        const auto weight = qf.getWeight(q);
        
        // Map to physical coordinates
        const auto physPoint = element.getTransformation().asFunction()(refPoint);
        
        // Backtrace to find departure point
        const auto departurePoint = backtrace(physPoint, ScalarType(0)); // Assuming current time = 0
        
        // Evaluate old solution at departure point
        ScalarType oldValue = evaluateOldSolution(departurePoint);
        
        // Apply conservative form correction if enabled
        if (m_conservativeForm)
        {
          const auto conservativeFactor = computeConservativeFactor(physPoint, departurePoint, ScalarType(0));
          oldValue *= conservativeFactor;
        }
        
        // Check for inflow boundary
        if (isInflowBoundary(departurePoint) && m_inflowBoundaryData)
        {
          oldValue = m_inflowBoundaryData(departurePoint, ScalarType(0));
        }
        
        // Compute metric factor
        const auto metricFactor = computeMetricFactor(elementIndex, refPoint);
        
        // Evaluate basis functions at quadrature point
        for (size_t i = 0; i < elementDOFs.size(); ++i)
        {
          const auto basisValue = m_fes.getBasisFunction(elementDOFs[i])(refPoint);
          const auto globalDOF = elementDOFs[i];
          
          // Accumulate RHS contribution
          rhs[globalDOF] += weight * metricFactor * oldValue * basisValue;
        }
      }
    }
  }

  template <class FES, class ScalarType>
  typename SemiLagrangianAdvection<FES, ScalarType>::PointType 
  SemiLagrangianAdvection<FES, ScalarType>::backtrace(const PointType& startPoint, ScalarType time)
  {
    // Initialize current state
    SpatialVectorType currentRefPoint = startPoint.getCoordinates(); // Convert to reference coordinates
    size_t currentElement = 0; // Find element containing start point
    ScalarType remainingTime = m_timestep;
    ScalarType currentTime = time;
    
    // Find initial element containing the start point
    for (auto it = m_mesh.getPolytopeIterator(m_mesh.getMeshDimension()); it; ++it)
    {
      if (it->contains(startPoint))
      {
        currentElement = it->getIndex();
        // Transform to reference coordinates
        currentRefPoint = it->getTransformation().inverse()(startPoint.getCoordinates());
        break;
      }
    }
    
    // Backtrace until time is consumed
    while (remainingTime > ScalarType(1e-12))
    {
      // Get current velocity in reference coordinates
      const auto physPoint = m_mesh.getPolytope(m_mesh.getMeshDimension(), currentElement)
                              ->getTransformation().asFunction()(currentRefPoint);
      auto velocity = m_velocityField(PointType(physPoint));
      
      // Project to tangent space if surface mesh
      if (m_mesh.isSurface())
      {
        const auto normal = m_mesh.getNormal(currentElement, currentRefPoint);
        velocity = projectToTangentSpace(velocity, normal);
      }
      
      // Transform velocity to reference coordinates
      const auto& jacobianInverse = m_mesh.getPolytope(m_mesh.getMeshDimension(), currentElement)
                                    ->getTransformation().getJacobianInverse(currentRefPoint);
      const auto refVelocity = jacobianInverse * velocity;
      
      // Compute exit time from current element
      const auto [exitTime, exitFace] = computeExitTime(currentRefPoint, refVelocity, currentElement);
      
      // Limit substep by exit time and remaining time
      const auto substepTime = std::min({remainingTime, m_exitTimeFactor * exitTime});
      
      // Perform integration substep
      std::pair<SpatialVectorType, size_t> result;
      if (m_integratorType == IntegratorType::RK2)
      {
        result = rk2Substep(currentRefPoint, currentElement, substepTime, currentTime);
      }
      else
      {
        result = rk3Substep(currentRefPoint, currentElement, substepTime, currentTime);
      }
      
      currentRefPoint = result.first;
      
      // Check if we need to hop to neighboring element
      if (substepTime >= m_exitTimeFactor * exitTime && exitTime > ScalarType(1e-12))
      {
        const auto [newElement, newRefPoint] = hopToNeighbor(currentRefPoint, currentElement, exitFace);
        currentElement = newElement;
        currentRefPoint = newRefPoint;
      }
      
      // Update time
      remainingTime -= substepTime;
      currentTime -= substepTime; // Backtracing goes backward in time
    }
    
    // Convert back to physical coordinates
    const auto finalPhysPoint = m_mesh.getPolytope(m_mesh.getMeshDimension(), currentElement)
                                ->getTransformation().asFunction()(currentRefPoint);
    return PointType(finalPhysPoint);
  }

  template <class FES, class ScalarType>
  std::pair<typename SemiLagrangianAdvection<FES, ScalarType>::SpatialVectorType, size_t>
  SemiLagrangianAdvection<FES, ScalarType>::rk2Substep(
      const SpatialVectorType& point,
      size_t element,
      ScalarType dt,
      ScalarType time)
  {
    // RK2 implementation for characteristic tracing
    // k1 = f(t, y)
    const auto physPoint1 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                             ->getTransformation().asFunction()(point);
    auto velocity1 = m_velocityField(PointType(physPoint1));
    
    if (m_mesh.isSurface())
    {
      const auto normal = m_mesh.getNormal(element, point);
      velocity1 = projectToTangentSpace(velocity1, normal);
    }
    
    const auto& jacobianInverse1 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                                   ->getTransformation().getJacobianInverse(point);
    const auto k1 = jacobianInverse1 * velocity1;
    
    // k2 = f(t + dt, y + dt*k1)
    const auto midPoint = point - dt * k1; // Negative because we're backtracing
    const auto physPoint2 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                             ->getTransformation().asFunction()(midPoint);
    auto velocity2 = m_velocityField(PointType(physPoint2));
    
    if (m_mesh.isSurface())
    {
      const auto normal = m_mesh.getNormal(element, midPoint);
      velocity2 = projectToTangentSpace(velocity2, normal);
    }
    
    const auto& jacobianInverse2 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                                   ->getTransformation().getJacobianInverse(midPoint);
    const auto k2 = jacobianInverse2 * velocity2;
    
    // Final update: y_{n+1} = y_n + dt/2 * (k1 + k2)
    const auto newPoint = point - dt * ScalarType(0.5) * (k1 + k2);
    
    return {newPoint, element};
  }

  template <class FES, class ScalarType>
  std::pair<typename SemiLagrangianAdvection<FES, ScalarType>::SpatialVectorType, size_t>
  SemiLagrangianAdvection<FES, ScalarType>::rk3Substep(
      const SpatialVectorType& point,
      size_t element,
      ScalarType dt,
      ScalarType time)
  {
    // RK3 implementation for characteristic tracing
    // k1 = f(t, y)
    const auto physPoint1 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                             ->getTransformation().asFunction()(point);
    auto velocity1 = m_velocityField(PointType(physPoint1));
    
    if (m_mesh.isSurface())
    {
      const auto normal = m_mesh.getNormal(element, point);
      velocity1 = projectToTangentSpace(velocity1, normal);
    }
    
    const auto& jacobianInverse1 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                                   ->getTransformation().getJacobianInverse(point);
    const auto k1 = jacobianInverse1 * velocity1;
    
    // k2 = f(t + dt/2, y + dt/2*k1)
    const auto midPoint1 = point - dt * ScalarType(0.5) * k1;
    const auto physPoint2 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                             ->getTransformation().asFunction()(midPoint1);
    auto velocity2 = m_velocityField(PointType(physPoint2));
    
    if (m_mesh.isSurface())
    {
      const auto normal = m_mesh.getNormal(element, midPoint1);
      velocity2 = projectToTangentSpace(velocity2, normal);
    }
    
    const auto& jacobianInverse2 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                                   ->getTransformation().getJacobianInverse(midPoint1);
    const auto k2 = jacobianInverse2 * velocity2;
    
    // k3 = f(t + dt, y - dt*k1 + 2*dt*k2)
    const auto midPoint2 = point - dt * k1 + ScalarType(2) * dt * k2;
    const auto physPoint3 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                             ->getTransformation().asFunction()(midPoint2);
    auto velocity3 = m_velocityField(PointType(physPoint3));
    
    if (m_mesh.isSurface())
    {
      const auto normal = m_mesh.getNormal(element, midPoint2);
      velocity3 = projectToTangentSpace(velocity3, normal);
    }
    
    const auto& jacobianInverse3 = m_mesh.getPolytope(m_mesh.getMeshDimension(), element)
                                   ->getTransformation().getJacobianInverse(midPoint2);
    const auto k3 = jacobianInverse3 * velocity3;
    
    // Final update: y_{n+1} = y_n + dt/6 * (k1 + 4*k2 + k3)
    const auto newPoint = point - dt / ScalarType(6) * (k1 + ScalarType(4) * k2 + k3);
    
    return {newPoint, element};
  }

  template <class FES, class ScalarType>
  std::pair<ScalarType, size_t>
  SemiLagrangianAdvection<FES, ScalarType>::computeExitTime(
      const SpatialVectorType& refPoint,
      const SpatialVectorType& velocity,
      size_t element)
  {
    // For simplicity, implement for simplex elements first
    // This computes when a ray from refPoint in direction -velocity 
    // exits the reference element using barycentric coordinates
    
    const auto polytope = m_mesh.getPolytope(m_mesh.getMeshDimension(), element);
    const auto geometry = polytope->getGeometry();
    
    ScalarType minExitTime = std::numeric_limits<ScalarType>::max();
    size_t exitFace = 0;
    
    if (geometry == Geometry::Polytope::Type::Triangle || 
        geometry == Geometry::Polytope::Type::Tetrahedron)
    {
      // For simplices, use barycentric coordinates
      // Each face is defined by one barycentric coordinate = 0
      const size_t dim = (geometry == Geometry::Polytope::Type::Triangle) ? 2 : 3;
      
      for (size_t face = 0; face <= dim; ++face)
      {
        // Check intersection with face where barycentric coordinate face = 0
        if (std::abs(velocity[face]) > ScalarType(1e-12))
        {
          const ScalarType t = -refPoint[face] / velocity[face];
          if (t > ScalarType(0) && t < minExitTime)
          {
            minExitTime = t;
            exitFace = face;
          }
        }
      }
    }
    else
    {
      // For tensor product elements (quad, hex), check each coordinate bound
      const size_t dim = (geometry == Geometry::Polytope::Type::Quadrilateral) ? 2 : 3;
      
      for (size_t d = 0; d < dim; ++d)
      {
        // Check exit through lower bound (coord = -1)
        if (velocity[d] < -ScalarType(1e-12))
        {
          const ScalarType t = (-ScalarType(1) - refPoint[d]) / velocity[d];
          if (t > ScalarType(0) && t < minExitTime)
          {
            minExitTime = t;
            exitFace = 2 * d; // Lower face
          }
        }
        
        // Check exit through upper bound (coord = +1)
        if (velocity[d] > ScalarType(1e-12))
        {
          const ScalarType t = (ScalarType(1) - refPoint[d]) / velocity[d];
          if (t > ScalarType(0) && t < minExitTime)
          {
            minExitTime = t;
            exitFace = 2 * d + 1; // Upper face
          }
        }
      }
    }
    
    return {minExitTime, exitFace};
  }

  template <class FES, class ScalarType>
  std::pair<size_t, typename SemiLagrangianAdvection<FES, ScalarType>::SpatialVectorType>
  SemiLagrangianAdvection<FES, ScalarType>::hopToNeighbor(
      const SpatialVectorType& refPoint,
      size_t currentElement,
      size_t exitFace)
  {
    // Find neighboring element through face connectivity
    const auto& connectivity = m_mesh.getConnectivity();
    const size_t meshDim = m_mesh.getMeshDimension();
    const size_t faceDim = meshDim - 1;
    
    // Get face indices for current element
    const auto& elementFaces = connectivity.getIncidence({meshDim, faceDim}, currentElement);
    
    if (exitFace < elementFaces.size())
    {
      const auto faceIndex = elementFaces[exitFace];
      
      // Find elements sharing this face
      const auto& faceElements = connectivity.getIncidence({faceDim, meshDim}, faceIndex);
      
      // Find the other element (not the current one)
      for (const auto& neighborIndex : faceElements)
      {
        if (neighborIndex != currentElement)
        {
          // Transform point to neighbor's reference coordinates
          // This is a simplified version - in practice would need proper face mapping
          return {neighborIndex, refPoint};
        }
      }
    }
    
    // If no neighbor found, stay in current element
    return {currentElement, refPoint};
  }

  template <class FES, class ScalarType>
  ScalarType SemiLagrangianAdvection<FES, ScalarType>::evaluateOldSolution(const PointType& departurePoint)
  {
    // Evaluate the old solution u_h^n(x*) = sum_j c_j^n * φ_j(x*)
    return m_oldSolution(departurePoint);
  }

  template <class FES, class ScalarType>
  ScalarType SemiLagrangianAdvection<FES, ScalarType>::computeMetricFactor(
      size_t element,
      const SpatialVectorType& refPoint)
  {
    const auto polytope = m_mesh.getPolytope(m_mesh.getMeshDimension(), element);
    const auto& transformation = polytope->getTransformation();
    const auto jacobian = transformation.getJacobian(refPoint);
    
    if (m_mesh.isSurface())
    {
      // For surface meshes: sqrt(det(J^T J))
      const auto gramMatrix = jacobian.transpose() * jacobian;
      return std::sqrt(gramMatrix.determinant());
    }
    else
    {
      // For volume meshes: det(J)
      return std::abs(jacobian.determinant());
    }
  }

  template <class FES, class ScalarType>
  typename SemiLagrangianAdvection<FES, ScalarType>::SpatialVectorType
  SemiLagrangianAdvection<FES, ScalarType>::projectToTangentSpace(
      const SpatialVectorType& velocity,
      const SpatialVectorType& normal)
  {
    // Tangential projection: βτ = (I - n⊗n)β = β - (β·n)n
    const auto normalComponent = velocity.dot(normal);
    return velocity - normalComponent * normal;
  }

  template <class FES, class ScalarType>
  ScalarType SemiLagrangianAdvection<FES, ScalarType>::computeConservativeFactor(
      const PointType& startPoint,
      const PointType& endPoint,
      ScalarType time)
  {
    // Simplified implementation: compute divergence at midpoint
    const auto midPoint = PointType((startPoint.getCoordinates() + endPoint.getCoordinates()) * ScalarType(0.5));
    
    // Compute div(βτ) at midpoint - this would need proper implementation
    // For now, return 1.0 (no divergence correction)
    return ScalarType(1.0);
  }

  template <class FES, class ScalarType>
  bool SemiLagrangianAdvection<FES, ScalarType>::isInflowBoundary(const PointType& point)
  {
    // Check if point is on mesh boundary and has inflow velocity
    // This is a simplified implementation
    return false; // For now, assume no inflow boundaries
  }

  template <class FES, class ScalarType>
  void SemiLagrangianAdvection<FES, ScalarType>::solveMassMatrixSystem(
      const VectorType& rhs,
      VectorType& solution)
  {
    if (!m_massMatrix.has_value())
    {
      // If no mass matrix provided, assume identity (lumped with unity diagonal)
      solution = rhs;
      return;
    }
    
    const auto& M = m_massMatrix.value().get();
    
    if (m_massMatrixType == MassMatrixType::Lumped)
    {
      // For lumped mass matrix, divide element-wise by diagonal
      solution.resize(rhs.size());
      for (size_t i = 0; i < rhs.size(); ++i)
      {
        solution[i] = rhs[i] / M(i, i);
      }
    }
    else
    {
      // For consistent mass matrix, solve linear system
      // This would typically use a sparse solver
      solution = M.lu().solve(rhs);
    }
  }
}

#endif