/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODININTEGRATION_MMG_SOLUTION_H
#define RODIN_RODININTEGRATION_MMG_SOLUTION_H

#include <cassert>
#include <boost/filesystem.hpp>

#include "Common.h"

namespace Rodin::External::MMG
{
  /**
   * @brief Base class for MMG solution objects.
   * @tparam Dimension Mesh dimension
   *
   * This class wraps the functionality of the type MMG5_Sol.
   */
  class SolutionBase
  {
    public:
      SolutionBase(MMG5_pSol sol = nullptr);

      SolutionBase(SolutionBase&& other);

      SolutionBase(const SolutionBase& other);

      SolutionBase& operator=(const SolutionBase& other);

      SolutionBase& operator=(SolutionBase&& other);

      virtual ~SolutionBase();

      /**
       * @brief Gets the dimension of the solution support.
       */
      int getDimension() const
      {
         return getHandle()->dim;
      }

      /**
       * @returns Number of points of the solution.
       */
      int count() const
      {
         return getHandle()->np;
      }

      /**
       * @returns Number of solutions per entity.
       */
      int size() const
      {
         return getHandle()->size;
      }

      /**
       * @internal
       * @returns Reference to underlying solution handle.
       */
      MMG5_pSol& getHandle()
      {
         return m_sol;
      }

      /**
       * @internal
       * @returns Constant reference to underlying solution handle.
       */
      const MMG5_pSol& getHandle() const
      {
         return m_sol;
      }

      virtual void save(const boost::filesystem::path& filename) = 0;
    private:
       MMG5_pSol m_sol;
  };

  class IncompleteSolutionBase
  {
    public:
      IncompleteSolutionBase();

      IncompleteSolutionBase(const IncompleteSolutionBase& other);

      IncompleteSolutionBase(IncompleteSolutionBase&& other);

      virtual ~IncompleteSolutionBase();

      MMG5_pSol release()
      {
         m_isOwner = false;
         auto sol = m_sol;
         m_sol = nullptr;
         return sol;
      }

      bool isOwner() const
      {
         return m_isOwner;
      }

      MMG5_pSol& getHandle()
      {
         return m_sol;
      }

      const MMG5_pSol& getHandle() const
      {
         return m_sol;
      }
    private:
      MMG5_pSol m_sol;
      bool m_isOwner;
  };
}

#endif
