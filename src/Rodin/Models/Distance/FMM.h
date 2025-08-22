/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_FMM_H
#define RODIN_MODELS_DISTANCE_FMM_H

#include <vector>
#include <atomic>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>

#include "Rodin/Solver/CG.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Multithreaded Fast Marching Method for computing distance functions on FEM meshes.
   * 
   * This implementation uses a bucketed priority queue approach for parallelization
   * and solves local Eikonal equations on each polytope using least squares gradient
   * estimation in physical coordinates.
   */
  class FMM
  {
    public:
      /**
       * @brief DOF status in Fast Marching Method
       */
      enum class Status : std::uint8_t
      {
        FAR,      ///< DOF not yet processed
        NARROW,   ///< DOF in narrow band (active front)
        ACCEPTED  ///< DOF with final arrival time
      };

    private:
      // Constants for numerical stability
      static constexpr Real eps = 1e-14;
      static constexpr Real tol = 1e-12; 
      static constexpr Real epsF = 1e-14;

      /**
       * @brief Get physical coordinates of a DOF
       */
      template <class FES>
      Math::SpatialVector<Real> getDOFCoord(const FES& fes, Index dof) const
      {
        const auto& mesh = fes.getMesh();
        const auto& vertices = mesh.getVertices();
        const size_t spaceDim = mesh.getSpaceDimension();
        
        // For P1 elements, DOF coordinates are vertex coordinates
        // TODO: For higher order elements, would need shape function evaluation
        Math::SpatialVector<Real> coord(spaceDim);
        if (dof < static_cast<Index>(vertices.cols()))
        {
          for (size_t d = 0; d < spaceDim; ++d)
            coord[d] = vertices(d, dof);
        }
        else
        {
          // Handle edge/face DOFs for higher order elements
          // For now, just return zero coordinates
          coord.setZero();
        }
        return coord;
      }

      /**
       * @brief Get neighboring DOFs for a given DOF
       */
      template <class FES>
      std::vector<Index> getDOFNeighbors(const FES& fes, Index dof) const
      {
        std::vector<Index> neighbors;
        const auto& mesh = fes.getMesh();
        const size_t meshDim = mesh.getDimension();
        
        // Get elements containing this DOF
        const auto elements = getDOFElements(fes, dof);
        
        FlatSet<Index> neighborSet;
        for (Index elemIdx : elements)
        {
          const auto elemDOFs = getElementDOFs(fes, meshDim, elemIdx);
          for (Index neighborDOF : elemDOFs)
          {
            if (neighborDOF != dof)
              neighborSet.insert(neighborDOF);
          }
        }
        
        neighbors.assign(neighborSet.begin(), neighborSet.end());
        return neighbors;
      }

      /**
       * @brief Get elements containing a given DOF
       */
      template <class FES>
      std::vector<Index> getDOFElements(const FES& fes, Index dof) const
      {
        std::vector<Index> elements;
        const auto& mesh = fes.getMesh();
        const size_t meshDim = mesh.getDimension();
        
        // For P1 elements, search through all cells
        const size_t numCells = mesh.getCellCount();
        for (Index cellIdx = 0; cellIdx < numCells; ++cellIdx)
        {
          const auto elemDOFs = getElementDOFs(fes, meshDim, cellIdx);
          if (std::find(elemDOFs.begin(), elemDOFs.end(), dof) != elemDOFs.end())
          {
            elements.push_back(cellIdx);
          }
        }
        
        return elements;
      }

      /**
       * @brief Get DOF indices for an element
       */
      template <class FES>
      std::vector<Index> getElementDOFs(const FES& fes, size_t dim, Index elemIdx) const
      {
        std::vector<Index> dofs;
        const auto& fe = fes.getFiniteElement(dim, elemIdx);
        const size_t numDOFs = fe.getCount();
        
        for (size_t local = 0; local < numDOFs; ++local)
        {
          const Index global = fes.getGlobalIndex({dim, elemIdx}, local);
          dofs.push_back(global);
        }
        
        return dofs;
      }

      /**
       * @brief Get speed field value at DOF (default to 1.0)
       */
      Real getSpeed(Index dof) const
      {
        // TODO: Allow user-specified speed field
        return 1.0;
      }

      /**
       * @brief Quantize arrival time to bucket index
       */
      size_t quantize(Real t, Real dt) const
      {
        if (!std::isfinite(t) || t < 0)
          return SIZE_MAX;
        return static_cast<size_t>(std::floor(t / dt));
      }

      /**
       * @brief Solve local Eikonal equation on a polytope
       */
      template <class FES>
      Real solveLocalEikonal(const FES& fes, Index elemIdx, Index targetDOF,
                            const std::vector<Real>& T, 
                            const std::vector<Status>& S) const
      {
        const auto& mesh = fes.getMesh();
        const size_t meshDim = mesh.getDimension();
        const size_t spaceDim = mesh.getSpaceDimension();
        
        // Get element DOFs and their properties
        const auto elemDOFs = getElementDOFs(fes, meshDim, elemIdx);
        const size_t numDOFs = elemDOFs.size();
        
        // Find target DOF in element
        auto it = std::find(elemDOFs.begin(), elemDOFs.end(), targetDOF);
        if (it == elemDOFs.end())
          return std::numeric_limits<Real>::infinity();
        
        const Math::SpatialVector<Real> x0 = getDOFCoord(fes, targetDOF);
        const Real T0 = T[targetDOF];
        const Real F0 = getSpeed(targetDOF);
        
        // Collect accepted DOFs in this element
        std::vector<Index> acceptedDOFs;
        std::vector<Math::SpatialVector<Real>> acceptedCoords;
        std::vector<Real> acceptedTimes;
        Real maxAcceptedTime = -std::numeric_limits<Real>::infinity();
        
        for (Index dof : elemDOFs)
        {
          if (dof != targetDOF && S[dof] == Status::ACCEPTED)
          {
            acceptedDOFs.push_back(dof);
            acceptedCoords.push_back(getDOFCoord(fes, dof));
            acceptedTimes.push_back(T[dof]);
            maxAcceptedTime = std::max(maxAcceptedTime, T[dof]);
          }
        }
        
        if (acceptedDOFs.size() < spaceDim)
          return std::numeric_limits<Real>::infinity();
        
        // Least squares gradient estimation
        const size_t n = acceptedDOFs.size();
        if (n == 0)
          return std::numeric_limits<Real>::infinity();
          
        Eigen::MatrixXd A(n, spaceDim);
        Eigen::VectorXd b(n);
        
        for (size_t i = 0; i < n; ++i)
        {
          const Math::SpatialVector<Real> dx = acceptedCoords[i] - x0;
          for (size_t d = 0; d < spaceDim; ++d)
            A(i, d) = dx[d];
          b(i) = acceptedTimes[i] - T0;
        }
        
        // Solve least squares problem
        const auto AtA = A.transpose() * A;
        const auto Atb = A.transpose() * b;
        
        // Check if AtA is well-conditioned
        const Real conditionThreshold = 1e-12;
        Eigen::LDLT<Eigen::MatrixXd> solver(AtA);
        if (solver.info() != Eigen::Success)
          return std::numeric_limits<Real>::infinity();
        
        // Additional condition number check
        const auto AtA_diag = AtA.diagonal();
        const Real minDiag = AtA_diag.minCoeff();
        const Real maxDiag = AtA_diag.maxCoeff();
        if (minDiag < conditionThreshold || maxDiag / minDiag > 1e12)
          return std::numeric_limits<Real>::infinity();
          
        const Eigen::VectorXd g = solver.solve(Atb);
        const Real gNorm = g.norm();
        
        if (!std::isfinite(gNorm) || gNorm < eps)
          return std::numeric_limits<Real>::infinity();
        
        // Construct candidate arrival times using upwind projection
        Real tentativeTime = std::numeric_limits<Real>::infinity();
        
        for (size_t i = 0; i < n; ++i)
        {
          const Math::SpatialVector<Real> dx = acceptedCoords[i] - x0;
          const Real proj = dx.dot(g) / gNorm;
          const Real candidate = acceptedTimes[i] + std::abs(proj) / std::max(F0, epsF);
          
          // Check causality
          if (candidate + tol >= acceptedTimes[i] && candidate + tol >= maxAcceptedTime)
          {
            tentativeTime = std::min(tentativeTime, candidate);
          }
        }
        
        return tentativeTime;
      }

      /**
       * @brief Compute local update for a DOF
       */
      template <class FES>
      Real localUpdate(const FES& fes, Index dof,
                      const std::vector<Real>& T,
                      const std::vector<Status>& S) const
      {
        const auto elements = getDOFElements(fes, dof);
        Real minTime = std::numeric_limits<Real>::infinity();
        
        for (Index elemIdx : elements)
        {
          const Real candidateTime = solveLocalEikonal(fes, elemIdx, dof, T, S);
          if (std::isfinite(candidateTime))
            minTime = std::min(minTime, candidateTime);
        }
        
        return minTime;
      }

    public:
      // Legacy enum for backward compatibility
      enum class Label
      {
        Far = static_cast<int>(Status::FAR),
        Considered = static_cast<int>(Status::NARROW), 
        Accepted = static_cast<int>(Status::ACCEPTED)
      };

      template <class FES>
      auto operator()(
          Geometry::Attribute interface,
          Geometry::Attribute region,
          const FES& fes) const
      {
        return operator()(
            FlatSet<Geometry::Attribute>{interface},
            FlatSet<Geometry::Attribute>{region},
            fes);
      }

      template <class FES>
      auto operator()(
          const FlatSet<Geometry::Attribute>& interface,
          const FlatSet<Geometry::Attribute>& domain,
          const FES& fes) const
      {
        using RangeType = typename FormLanguage::Traits<FES>::RangeType;
        static_assert(std::is_same_v<RangeType, Real>);
        static_assert(std::numeric_limits<Real>::has_infinity);
        
        const size_t numDOFs = fes.getSize();
        const auto& mesh = fes.getMesh();
        const auto& meshDim = mesh.getDimension();
        
        // Initialize data structures
        std::vector<std::atomic<Real>> T(numDOFs);
        std::vector<std::atomic<Status>> S(numDOFs);
        
        // Initialize all DOFs
        for (size_t i = 0; i < numDOFs; ++i)
        {
          T[i].store(std::numeric_limits<Real>::infinity());
          S[i].store(Status::FAR);
        }
        
        // Set boundary conditions on interface
        std::vector<Index> initialDOFs;
        for (Index i = 0; i < mesh.getFaceCount(); i++)
        {
          if (interface.count(mesh.getAttribute(meshDim - 1, i)))
          {
            const auto& fe = fes.getFiniteElement(meshDim - 1, i);
            for (size_t local = 0; local < fe.getCount(); local++)
            {
              const Index global = fes.getGlobalIndex({ meshDim - 1, i }, local);
              S[global].store(Status::ACCEPTED);
              T[global].store(0.0);
              initialDOFs.push_back(global);
            }
          }
        }
        
        // Estimate discretization parameter for bucketing
        Real minEdgeLength = std::numeric_limits<Real>::max();
        Real minSpeed = std::numeric_limits<Real>::max();
        
        // Simple estimation based on mesh size
        if (numDOFs > 1)
        {
          const auto& vertices = mesh.getVertices();
          const size_t spaceDim = mesh.getSpaceDimension();
          const size_t numVertices = std::min(numDOFs, static_cast<Index>(vertices.cols()));
          
          // Sample a few vertices to estimate typical edge length
          const size_t maxSamples = std::min(static_cast<size_t>(100), numVertices);
          for (size_t i = 0; i < maxSamples && i < numVertices; ++i)
          {
            for (size_t j = i + 1; j < maxSamples && j < numVertices; ++j)
            {
              Real dist = 0.0;
              for (size_t d = 0; d < spaceDim; ++d)
              {
                const Real diff = vertices(d, i) - vertices(d, j);
                dist += diff * diff;
              }
              dist = std::sqrt(dist);
              if (dist > 1e-10)  // Avoid degenerate cases
                minEdgeLength = std::min(minEdgeLength, dist);
            }
            minSpeed = std::min(minSpeed, getSpeed(i));
          }
        }
        
        // Set reasonable defaults if estimation fails
        if (minEdgeLength == std::numeric_limits<Real>::max() || minEdgeLength < 1e-10)
          minEdgeLength = 1.0;
        if (minSpeed == std::numeric_limits<Real>::max() || minSpeed < 1e-10)
          minSpeed = 1.0;
        
        const Real dt = minEdgeLength / (4.0 * minSpeed);  // Conservative discretization
        
        // Bucketed priority queues
        std::vector<std::vector<Index>> buckets;
        size_t maxBucket = 0;
        
        // Create non-atomic copies for local computations
        std::vector<Real> T_local(numDOFs);
        std::vector<Status> S_local(numDOFs);
        
        for (size_t i = 0; i < numDOFs; ++i)
        {
          T_local[i] = T[i].load();
          S_local[i] = S[i].load();
        }
        
        // Initialize narrow band by updating neighbors of accepted DOFs
        std::set<Index> narrowDOFs;
        for (Index acceptedDOF : initialDOFs)
        {
          const auto neighbors = getDOFNeighbors(fes, acceptedDOF);
          for (Index neighbor : neighbors)
          {
            if (S_local[neighbor] == Status::FAR)
            {
              const Real tentTime = localUpdate(fes, neighbor, T_local, S_local);
              if (std::isfinite(tentTime))
              {
                T_local[neighbor] = tentTime;
                S_local[neighbor] = Status::NARROW;
                narrowDOFs.insert(neighbor);
                
                const size_t bucket = quantize(tentTime, dt);
                if (bucket != SIZE_MAX)
                {
                  maxBucket = std::max(maxBucket, bucket);
                  if (bucket >= buckets.size())
                    buckets.resize(bucket + 1);
                  buckets[bucket].push_back(neighbor);
                }
              }
            }
          }
        }
        
        // Update atomic arrays with initial narrow band
        for (size_t i = 0; i < numDOFs; ++i)
        {
          T[i].store(T_local[i]);
          S[i].store(S_local[i]);
        }
        
        // Main FMM loop over buckets
        const size_t maxIterations = numDOFs * 10;  // Safety limit
        size_t iterations = 0;
        
        for (size_t b = 0; b <= maxBucket && iterations < maxIterations; ++b)
        {
          if (b >= buckets.size() || buckets[b].empty())
            continue;
            
          // Take snapshot of current bucket
          std::vector<Index> currentNodes = std::move(buckets[b]);
          buckets[b].clear();
          
          // Thread-local storage for newly relaxed nodes
          std::vector<std::vector<Index>> threadLocalRelaxed;
          
#ifdef _OPENMP
          const int numThreads = omp_get_max_threads();
          threadLocalRelaxed.resize(numThreads);
          
          #pragma omp parallel for
#else
          threadLocalRelaxed.resize(1);
#endif
          for (size_t idx = 0; idx < currentNodes.size(); ++idx)
          {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            const Index u = currentNodes[idx];
            
            // Accept current node if not already accepted
            Status currentStatus = S[u].load();
            if (currentStatus == Status::ACCEPTED)
              continue;
              
            S[u].store(Status::ACCEPTED);
            
            // Update neighbors
            const auto neighbors = getDOFNeighbors(fes, u);
            for (Index v : neighbors)
            {
              if (S[v].load() != Status::ACCEPTED)
              {
                // Create local copies for computation
                std::vector<Real> T_thread(numDOFs);
                std::vector<Status> S_thread(numDOFs);
                
                for (size_t i = 0; i < numDOFs; ++i)
                {
                  T_thread[i] = T[i].load();
                  S_thread[i] = S[i].load();
                }
                
                const Real tentTime = localUpdate(fes, v, T_thread, S_thread);
                const Real currentTime = T[v].load();
                
                if (std::isfinite(tentTime) && tentTime < currentTime)
                {
                  // Atomic relaxation
                  Real expected = currentTime;
                  while (tentTime < expected && 
                         !T[v].compare_exchange_weak(expected, tentTime)) 
                  {
                    // Retry if another thread updated in the meantime
                  }
                  
                  S[v].store(Status::NARROW);
                  threadLocalRelaxed[tid].push_back(v);
                }
              }
            }
          }
          
          iterations++;
          
          // Merge thread-local relaxed nodes and re-bucket
          for (const auto& relaxedNodes : threadLocalRelaxed)
          {
            for (Index v : relaxedNodes)
            {
              const Real newTime = T[v].load();
              const size_t newBucket = quantize(newTime, dt);
              
              if (newBucket != SIZE_MAX && newBucket < numDOFs * 100)  // Sanity check
              {
                maxBucket = std::max(maxBucket, newBucket);
                if (newBucket >= buckets.size())
                  buckets.resize(newBucket + 1);
                buckets[newBucket].push_back(v);
              }
            }
          }
        }
        
        // Create result GridFunction
        Variational::GridFunction u(fes);
        for (size_t i = 0; i < numDOFs; ++i)
        {
          u[i] = T[i].load();
        }
        
        return u;
      }
  };
}

#endif



