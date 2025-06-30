/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H

#include <petsc.h>
#include <petscvec.h>

#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Variational
{
  template <class FES>
  class GridFunction<FES, ::Vec> final
    : public GridFunctionBase<GridFunction<FES, ::Vec>>
  {
    public:
      using FESType = FES;

      using ScalarType = PetscScalar;

      using DataType = ::Vec;

      using Parent = GridFunctionBase<GridFunction<FESType, DataType>>;
      using Parent::operator=;
      using Parent::min;
      using Parent::max;

      GridFunction(const FESType& fes, DataType& data)
        : Parent(fes, data)
      {
        const auto& mesh = fes.getMesh();
        const auto& comm = mesh.getCommunicator();
        size_t begin, end;
        fes.getOwnershipRange(begin, end);
        size_t ownedSize = end - begin;
        size_t ghostSize = fes.getShard().getSize() - ownedSize;
        // VecCreateGhost(comm, fes.getShard().getOwned().size(), fes.getSize(), );
      }

      GridFunction(const FESType& fes, DataType&& _data)
        : Parent(fes, std::move(_data))
      {
        auto& data = this->getData();
        data.resize(fes.getSize());
        data.setZero();
      }

      constexpr
      ScalarType min(Index& idx) const
      {
        return this->getData().minCoeff(&idx);
      }

      constexpr
      ScalarType max(Index& idx) const
      {
        return this->getData().maxCoeff(&idx);
      }

      ScalarType& operator[](Index global)
      {
        return this->getData()[global];
      }

      const ScalarType& operator[](Index global) const
      {
        return this->getData()[global];
      }

      GridFunction& operator+=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->getData() += rhs;
        return *this;
      }

      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->getData() -= rhs;
        return *this;
      }

      GridFunction& operator*=(const ScalarType& rhs)
      {
        this->getData() *= rhs;
        return *this;
      }

      GridFunction& operator/=(const ScalarType& rhs)
      {
        auto& data = this->getData();
        data = data.array() / rhs;
        return static_cast<GridFunction&>(*this);
      }

      GridFunction& operator+=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->getData().array() += rhs.getData().array();
        return *this;
      }

      GridFunction& operator-=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->getData().array() -= rhs.getData().array();
        return *this;
      }

      GridFunction& operator*=(const GridFunction& rhs)
      {
        this->getData().array() *= rhs.getData().array();
        return *this;
      }

      GridFunction& operator/=(const GridFunction& rhs)
      {
        this->getData().array() /= rhs.getData().array();
        return *this;
      }

      GridFunction& setData(const DataType& data, size_t offset = 0)
      {
        const auto sz = this->getFiniteElementSpace().getSize();
        assert(offset + sz <= data.size());
        this->getData() = data.segment(offset, sz);
        return *this;
      }
  };
}

#endif
