/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_EIKONAL_FMM_H
#define RODIN_MODELS_EIKONAL_FMM_H

#include <queue>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>

#include "Rodin/Context/Local.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Models::Eikonal
{
  template <class Solution, class SpeedFunction, class Context>
  class FMM;

  /**
   * @brief Fast Marching Method implementation for solving the Eikonal equation.
   * 
   * This class implements the Fast Marching Method (FMM) algorithm for solving
   * the Eikonal equation |∇u| = F where F is the speed function and u represents
   * the arrival time or distance function.
   * 
   * The algorithm works by:
   * 1. Initializing all points as "Far" with infinite distance
   * 2. Setting interface/boundary points as "Accepted" with distance 0
   * 3. Adding immediate neighbors to "Considered" set with computed tentative values
   * 4. Iteratively processing the "Considered" point with minimum value:
   *    - Mark it as "Accepted"
   *    - Update its neighbors with better values if possible
   * 
   * @tparam Solution Type of the solution grid function (usually GridFunction)
   * @tparam SpeedFunction Type of the speed function callable
   */
  template <class Solution, class SpeedFunction>
  class FMM<Solution, SpeedFunction, Context::Local>
  {
    public:
      using ScalarType = Real;
      using Context = Context::Local;
      using Mesh = Geometry::Mesh<Context>;
      using FES = Variational::P1<ScalarType, Mesh>;

      using SolutionType = Solution;
      using SpeedFunctionType = SpeedFunction;

      enum class Label
      {
        Far,
        Considered,
        Accepted
      };

    private:
      struct PriorityQueueItem
      {
        Index nodeIndex;
        Real value;
        
        bool operator>(const PriorityQueueItem& other) const
        {
          return value > other.value;
        }
      };

      using PriorityQueue = std::priority_queue<PriorityQueueItem, 
                                                std::vector<PriorityQueueItem>, 
                                                std::greater<PriorityQueueItem>>;

    public:

      /**
       * @brief Constructor for the FMM solver.
       * @param u Reference to the solution grid function
       * @param speed Speed function that takes a spatial point and returns the speed value
       */
      template <class Callable>
      FMM(SolutionType& u, Callable&& speed)
        : m_u(u), m_speed(std::forward<Callable>(speed)), m_interfacePredicate(nullptr)
      {}

      /**
       * @brief Set the interface condition using a predicate function.
       * @param p Predicate function that takes a node index and returns true if it's on the interface
       * @return Reference to this FMM object for method chaining
       */
      template <class Pred>
      FMM& setInterface(Pred&& p)
      {
        m_interfacePredicate = std::forward<Pred>(p);
        return *this;
      }

      /**
       * @brief Solve the Eikonal equation using the Fast Marching Method.
       * 
       * This method implements the complete FMM algorithm:
       * - Initializes all nodes as Far with infinite distance
       * - Sets interface nodes as Accepted with value 0
       * - Iteratively processes nodes in order of increasing distance
       * - Updates neighbors using local Eikonal solver
       */
      void solve()
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& meshDim = mesh.getDimension();
        
        // Initialize all nodes as Far with infinite distance
        m_labels.resize(fes.getSize(), Label::Far);
        u = std::numeric_limits<Real>::infinity();
        
        // Set up priority queue for considered nodes
        PriorityQueue considered;
        
        // Initialize interface nodes as Accepted with value 0
        if (m_interfacePredicate)
        {
          for (Index i = 0; i < fes.getSize(); ++i)
          {
            if (m_interfacePredicate(i))
            {
              m_labels[i] = Label::Accepted;
              u[i] = 0.0;
              
              // Add neighbors of interface nodes to considered set
              addNeighborsToConsidered(i, u, fes, mesh, considered);
            }
          }
        }
        
        // Main Fast Marching loop
        while (!considered.empty())
        {
          // Get node with minimum tentative value
          PriorityQueueItem current = considered.top();
          considered.pop();
          
          Index nodeIdx = current.nodeIndex;
          
          // Skip if already processed (can happen due to multiple updates)
          if (m_labels[nodeIdx] == Label::Accepted)
            continue;
            
          // Mark as accepted
          m_labels[nodeIdx] = Label::Accepted;
          
          // Update neighbors
          addNeighborsToConsidered(nodeIdx, u, fes, mesh, considered);
        }
      }

    private:
      void addNeighborsToConsidered(Index nodeIdx, SolutionType& u, 
                                    const FES& fes, const Mesh& mesh, 
                                    PriorityQueue& considered)
      {
        // Get all elements connected to this vertex
        const auto& meshDim = mesh.getDimension();
        
        // For each element connected to this vertex
        for (auto elemIt = mesh.getPolytope(meshDim); elemIt; ++elemIt)
        {
          const auto& element = *elemIt;
          const auto& vertices = element.getVertices();
          
          // Check if this element contains our vertex
          bool containsVertex = false;
          for (auto v : vertices)
          {
            if (v == nodeIdx)
            {
              containsVertex = true;
              break;
            }
          }
          
          if (!containsVertex)
            continue;
            
          // Update all vertices of this element
          for (auto vertexIdx : vertices)
          {
            if (m_labels[vertexIdx] == Label::Far)
            {
              // Compute tentative value using local Eikonal solver
              Real tentativeValue = computeTentativeValue(vertexIdx, element, u, fes, mesh);
              
              if (tentativeValue < u[vertexIdx])
              {
                u[vertexIdx] = tentativeValue;
                m_labels[vertexIdx] = Label::Considered;
                considered.push({vertexIdx, tentativeValue});
              }
            }
            else if (m_labels[vertexIdx] == Label::Considered)
            {
              // Update existing considered node if we found a better value
              Real tentativeValue = computeTentativeValue(vertexIdx, element, u, fes, mesh);
              
              if (tentativeValue < u[vertexIdx])
              {
                u[vertexIdx] = tentativeValue;
                considered.push({vertexIdx, tentativeValue});
              }
            }
          }
        }
      }
      
      template<class Element>
      Real computeTentativeValue(Index nodeIdx, const Element& element, 
                                 const SolutionType& u, const FES& fes, 
                                 const Mesh& mesh)
      {
        const auto& vertices = element.getVertices();
        const Real speedValue = m_speed(mesh.getVertexCoordinates(nodeIdx));
        
        if (speedValue <= 0.0)
          return std::numeric_limits<Real>::infinity();
          
        Real minValue = std::numeric_limits<Real>::infinity();
        
        // Find minimum value among accepted neighbors in this element
        for (auto vertexIdx : vertices)
        {
          if (vertexIdx != nodeIdx && m_labels[vertexIdx] == Label::Accepted)
          {
            const auto& coord1 = mesh.getVertexCoordinates(nodeIdx);
            const auto& coord2 = mesh.getVertexCoordinates(vertexIdx);
            
            // Simple upwind scheme: distance = neighbor_value + h/speed
            Real distance = (coord1 - coord2).norm();
            Real tentative = u[vertexIdx] + distance / speedValue;
            minValue = std::min(minValue, tentative);
          }
        }
        
        return minValue;
      }

    private:
      std::reference_wrapper<SolutionType> m_u;
      SpeedFunctionType m_speed;
      std::vector<Label> m_labels;
      std::function<bool(Index)> m_interfacePredicate;
  };
}

#endif




