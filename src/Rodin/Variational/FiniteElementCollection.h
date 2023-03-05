/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTCOLLECTION_H
#define RODIN_VARIATIONAL_FINITEELEMENTCOLLECTION_H

#include <mfem.hpp>

namespace Rodin::Variational
{
  /**
   * @brief Abstract base class for finite element collections.
   * @see L2
   * @see H1
   */
  class FiniteElementCollectionBase
  {
    public:
      FiniteElementCollectionBase() = default;

      FiniteElementCollectionBase(const FiniteElementCollectionBase&) = default;

      FiniteElementCollectionBase(FiniteElementCollectionBase&&) = default;

      FiniteElementCollectionBase& operator=(FiniteElementCollectionBase&&) = default;

      virtual ~FiniteElementCollectionBase() = default;

      inline
      size_t getOrder() const
      {
        return getHandle().GetOrder();
      }

      virtual mfem::FiniteElementCollection& getHandle() = 0;

      virtual const mfem::FiniteElementCollection& getHandle() const = 0;
  };
}

#endif
