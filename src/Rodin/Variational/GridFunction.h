/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_GRIDFUNCTION_H

#include <cmath>
#include <utility>
#include <fstream>
#include <functional>
#include <boost/filesystem.hpp>
#include <type_traits>

#include <mfem.hpp>

#include "Rodin/Math.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/IO/ForwardDecls.h"

#include "ForwardDecls.h"

#include "H1.h"
#include "L2.h"
#include "MFEM.h"
#include "Function.h"
#include "Component.h"
#include "Restriction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{

  /**
   * @defgroup GridFunctionSpecializations GridFunction Template Specializations
   * @brief Template specializations of the GridFunction class.
   * @see GridFunction
   */

  template <class Derived>
  class GridFunctionBase : public FunctionBase<GridFunctionBase<Derived>>
  {
    public:
      using Parent = FunctionBase<GridFunctionBase<Derived>>;

      constexpr
      GridFunctionBase() = default;

      constexpr
      GridFunctionBase(const GridFunctionBase& other)
        : Parent(other)
      {}

      constexpr
      GridFunctionBase(GridFunctionBase&& other)
        : Parent(std::move(other))
      {}

      GridFunctionBase& operator=(GridFunctionBase&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      GridFunctionBase& operator=(const GridFunctionBase&) = delete;

      /**
       * @brief Searches for the maximum value in the grid function data.
       * @returns Maximum value in grid function.
       *
       * This function will compute the maximum value in the grid function
       * data array.
       */
      inline
      constexpr
      Scalar max() const
      {
        return static_cast<const Derived&>(*this).max();
      }

      /**
       * @brief Searches the minimum value in the grid function data.
       * @returns Minimum value in grid function.
       *
       * This function will compute the minimum value in the grid function
       * data array.
       */
      inline
      constexpr
      Scalar min() const
      {
        return static_cast<const Derived&>(*this).max();
      }

      inline
      constexpr
      size_t getDimension() const
      {
        return getFiniteElementSpace().getVectorDimension();
      }

      inline
      constexpr
      auto x() const
      {
        assert(getFiniteElementSpace().getVectorDimension() >= 1);
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      auto y() const
      {
        assert(getFiniteElementSpace().getVectorDimension() >= 2);
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      auto z() const
      {
        assert(getFiniteElementSpace().getVectorDimension() >= 3);
        return static_cast<const Derived&>(*this);
      }

      inline
      constexpr
      Math::Vector& getData()
      {
        return static_cast<Derived&>(*this).getData();
      }

      inline
      constexpr
      const Math::Vector& getData() const
      {
        return static_cast<const Derived&>(*this).getData();
      }

      inline
      constexpr
      size_t getSize() const
      {
        return getHandle().Size();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getFiniteElementSpace().getVectorDimension(), 1 };
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      /**
       * @brief Bulk assigns the value to the whole data array.
       */
      inline
      Derived& operator=(Scalar v)
      {
        getData().setConstant(v);
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator=(std::function<Scalar(const Geometry::Point&)> fn)
      {
        using RangeType = typename FormLanguage::Traits<FunctionBase<GridFunctionBase<Derived>>>::RangeType;
        static_assert(std::is_same_v<RangeType, Scalar>);
        assert(getFiniteElementSpace().getVectorDimension() == 1);
        return project(ScalarFunction(fn));
      }

      /**
       * @brief Projection of a function.
       */
      template <class NestedDerived>
      inline
      Derived& operator=(const FunctionBase<NestedDerived>& fn)
      {
        return project(fn);
      }

      /**
       * @brief Addition of a scalar value.
       */
      inline
      Derived& operator+=(Scalar rhs)
      {
        getData() = getData().array() + rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Substraction of a scalar value.
       */
      inline
      Derived& operator-=(Scalar rhs)
      {
        getData() = getData().array() - rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Multiplication by a scalar value.
       */
      inline
      Derived& operator*=(Scalar rhs)
      {
        getData() = getData().array() * rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Division by a scalar value.
       */
      inline
      Derived& operator/=(Scalar rhs)
      {
        getData() = getData().array() / rhs;
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator+=(const GridFunctionBase& rhs)
      {
        if (this == &rhs)
        {
          operator*=(Scalar(2));
        }
        else
        {
          assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
          getData() = getData().array() + rhs.getData().array();
        }
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator-=(const GridFunctionBase& rhs)
      {
        if (this == &rhs)
        {
          operator=(Scalar(0));
        }
        else
        {
          assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
          getData() = getData().array() - rhs.getData().array();
        }
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator*=(const GridFunctionBase& rhs)
      {
        if (this == &rhs)
        {
          getData() = getData().array() * getData().array();
        }
        else
        {
          assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
          getData() = getData().array() * rhs.getData().array();
        }
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator/=(const GridFunctionBase& rhs)
      {
        if (this == &rhs)
        {
          operator=(Scalar(1));
        }
        else
        {
          assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
          getData() = getData().array() / rhs.getData().array();
        }
        return static_cast<Derived&>(*this);
      }

      inline
      auto& project(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return project(fn, std::set<Geometry::Attribute>{attr});
      }

      inline
      auto& project(std::function<Scalar(const Geometry::Point&)> fn, const std::set<Geometry::Attribute>& attrs = {})
      {
        return project(ScalarFunction(fn), attrs);
      }

      /**
       * @brief Projects a FunctionBase instance
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attribute.
       *
       * It is a convenience function to call
       * project(const FunctionBase&, const std::set<int>&) with one
       * attribute.
       */
      template <class NestedDerived>
      Derived& project(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return project(fn, std::set<Geometry::Attribute>{attr});
      }

      /**
       * @brief Projects a FunctionBase instance on the grid function.
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attributes. If the attribute set is
       * empty, this function will project over all elements in the mesh.
       */
      template <class NestedDerived>
      Derived& project(const FunctionBase<NestedDerived>& fn, const std::set<Geometry::Attribute>& attrs = {})
      {
        using Value = FunctionBase<NestedDerived>;
        using ValueRangeType = typename FormLanguage::Traits<Value>::RangeType;
        using RangeType = typename FormLanguage::Traits<FunctionBase<GridFunctionBase<Derived>>>::RangeType;
        static_assert(std::is_same_v<RangeType, ValueRangeType>);
        if constexpr (std::is_same_v<ValueRangeType, Scalar>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          Internal::MFEMScalarCoefficient sc(getFiniteElementSpace().getMesh(), fn);
          if (attrs.size() == 0)
          {
            getHandle().ProjectCoefficient(sc);
          }
          else
          {
            mfem::Array<int> vdofs;
            const auto& fes = getFiniteElementSpace().getHandle();
            for (int i = 0; i < fes.GetNE(); i++)
            {
              if (attrs.count(fes.GetAttribute(i)) > 0)
              {
                fes.GetElementVDofs(i, vdofs);
                getHandle().ProjectCoefficient(sc, vdofs);
              }
            }
          }
          return static_cast<Derived&>(*this);
        }
        else if constexpr (std::is_same_v<ValueRangeType, Math::Vector>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == fn.getRangeShape().height());
          Internal::MFEMVectorCoefficient vc(getFiniteElementSpace().getMesh(), fn);
          if (attrs.size() == 0)
          {
            getHandle().ProjectCoefficient(vc);
          }
          else
          {
            mfem::Array<int> vdofs;
            const auto& fes = getFiniteElementSpace().getHandle();
            for (int i = 0; i < fes.GetNE(); i++)
            {
              if (attrs.count(fes.GetAttribute(i)) > 0)
              {
                fes.GetElementVDofs(i, vdofs);
                getHandle().ProjectCoefficient(vc, vdofs);
              }
            }
          }
          return static_cast<Derived&>(*this);
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      inline
      void save(const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM,
          size_t precision = RODIN_DEFAULT_GRIDFUNCTION_SAVE_PRECISION) const
      {
        return static_cast<const Derived&>(*this).save(filename, fmt, precision);
      }

      inline
      Derived& load(const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        return static_cast<Derived&>(*this).save(filename, fmt);
      }

      inline
      constexpr
      auto& getFiniteElementSpace()
      {
        return static_cast<Derived&>(*this).getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return static_cast<const Derived&>(*this).getFiniteElementSpace();
      }

      virtual GridFunctionBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

      /**
       * @internal
       * @brief Gets the underlying handle to the mfem::GridFunction object.
       * @returns Reference to the underlying object.
       */
      mfem::GridFunction& getHandle()
      {
        return static_cast<Derived&>(*this).getHandle();
      }

      /**
       * @internal
       * @brief Gets the underlying handle to the mfem::GridFunction object.
       * @returns Constant reference to the underlying object.
       */
      const mfem::GridFunction& getHandle() const
      {
        return static_cast<const Derived&>(*this).getHandle();
      }
  };

  template <class FESType>
  class FESGridFunction : public GridFunctionBase<FESGridFunction<FESType>>
  {
    public:
      using FES = FESType;
      using Parent = GridFunctionBase<FESGridFunction<FES>>;

      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;

      constexpr
      FESGridFunction(const FES& fes)
        : m_fes(fes),
          m_data(fes.getHandle().GetVSize()),
          m_gf(new mfem::GridFunction(&m_fes.get().getHandle(), m_data.data()))
      {
        assert(!m_gf->OwnsData());
      }

      constexpr
      FESGridFunction(const FESGridFunction& other)
        : Parent(other),
          m_fes(other.m_fes),
          m_data(other.m_data),
          m_gf(new mfem::GridFunction(*other.m_gf))
      {}

      constexpr
      FESGridFunction(FESGridFunction&& other)
        : Parent(std::move(other)),
          m_fes(std::move(other.m_fes)),
          m_data(std::move(other.m_data)),
          m_gf(std::move(other.m_gf))
      {}

      inline
      auto getValue(const Geometry::Point& p) const
      {
        using RangeType = typename FES::RangeType;
        static_assert(std::is_same_v<RangeType, Scalar> || std::is_same_v<RangeType, Math::Vector>);
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          return Scalar(m_gf->GetValue(p.getTransformation().getHandle(), p.getIntegrationPoint()));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          Math::Vector res(getFiniteElementSpace().getVectorDimension());
          mfem::Vector tmp(res.data(), res.size());
          m_gf->GetVectorValue(p.getTransformation().getHandle(), p.getIntegrationPoint(), tmp);
          return res;
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      inline
      constexpr
      Math::Vector& getData()
      {
        return m_data;
      }

      inline
      constexpr
      const Math::Vector& getData() const
      {
        return m_data;
      }

      virtual FESGridFunction* copy() const noexcept override
      {
        return new FESGridFunction(*this);
      }

      inline
      mfem::GridFunction& getHandle() const
      {
        assert(m_gf);
        return *m_gf;
      }

    private:
      std::reference_wrapper<const FES> m_fes;
      Math::Vector m_data;
      std::unique_ptr<mfem::GridFunction> m_gf;
  };

  /**
   * @ingroup GridFunctionSpecializations
   * @brief Represents a GridFunction which belongs to an L2 finite element
   * space.
   */
  template <class ... Ts>
  class GridFunction<L2<Ts ...>> final : public FESGridFunction<L2<Ts ...>>
  {
    public:
      using FES = L2<Ts ...>;
      using Parent = FESGridFunction<FES>;

      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;

      /**
       * @brief Constructs a grid function on an L2 finite element space.
       * @param[in] fes Finite element space to which the function belongs
       * to.
       */
      constexpr
      GridFunction(const FES& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copies the grid function.
       * @param[in] other Other grid function to copy.
       */
      constexpr
      GridFunction(const GridFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructs the grid function.
       * @param[in] other Other grid function to move.
       */
      constexpr
      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
      {}

      GridFunction& operator=(const GridFunction&)  = delete;

      inline GridFunction* copy() const noexcept override
      {
        return new GridFunction(*this);
      }
  };

  template <class ... Ts>
  GridFunction(const L2<Ts ...>&) -> GridFunction<L2<Ts ...>>;

  /**
   * @ingroup GridFunctionSpecializations
   * @brief Represents a GridFunction which belongs to an H1 finite element
   * space.
   */
  template <class ... Ts>
  class GridFunction<H1<Ts...>> final : public FESGridFunction<H1<Ts...>>
  {
    public:
      using FES = H1<Ts...>;
      using Parent = FESGridFunction<FES>;

      using Parent::operator=;
      using Parent::operator+=;
      using Parent::operator-=;
      using Parent::operator*=;
      using Parent::operator/=;

      /**
       * @brief Constructs a grid function on a finite element space.
       * @param[in] fes Finite element space to which the function belongs
       * to.
       */
      constexpr
      GridFunction(const FES& fes)
        : Parent(fes)
      {}

      /**
       * @brief Copies the grid function.
       * @param[in] other Other grid function to copy.
       */
      constexpr
      GridFunction(const GridFunction& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructs the grid function.
       * @param[in] other Other grid function to move.
       */
      constexpr
      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Move assignment operator.
       */
      inline
      constexpr
      GridFunction& operator=(GridFunction&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      GridFunction& operator=(const GridFunction&)  = delete;

      template <class NestedDerived>
      inline
      GridFunction& projectOnBoundary(const FunctionBase<NestedDerived>& fn,
                                      Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, std::set<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      GridFunction& projectOnBoundary(const FunctionBase<NestedDerived>& fn,
                                      const std::set<Geometry::Attribute>& attrs = {})
      {
        using Value = FunctionBase<NestedDerived>;
        using ValueRangeType = typename FormLanguage::Traits<Value>::RangeType;
        if constexpr (std::is_same_v<ValueRangeType, Scalar>)
        {
          int maxBdrAttr = this->getFiniteElementSpace()
                                .getMesh()
                                .getHandle().bdr_attributes.Max();
          mfem::Array<int> marker(maxBdrAttr);
          if (attrs.size() == 0)
          {
            marker = 1;
            Internal::MFEMScalarCoefficient sc(this->getFiniteElementSpace().getMesh(), fn);
            this->getHandle().ProjectBdrCoefficient(sc, marker);
            return *this;
          }
          else
          {
            marker = 0;
            for (const auto& attr : attrs)
            {
              assert(attr - 1 < maxBdrAttr);
              marker[attr - 1] = 1;
            }
            Internal::MFEMScalarCoefficient sc(this->getFiniteElementSpace().getMesh(), fn);
            this->getHandle().ProjectBdrCoefficient(sc, marker);
            return *this;
          }
        }
        else if constexpr (std::is_same_v<ValueRangeType, Math::Vector>)
        {
          int maxBdrAttr = this->getFiniteElementSpace()
                                  .getMesh()
                                  .getHandle().bdr_attributes.Max();
          mfem::Array<int> marker(maxBdrAttr);
          if (attrs.size() == 0)
          {
            marker = 1;
            Internal::MFEMVectorCoefficient vc(this->getFiniteElementSpace().getMesh(), fn);
            this->getHandle().ProjectBdrCoefficient(vc, marker);
            return *this;
          }
          else
          {
            marker = 0;
            for (const auto& attr : attrs)
            {
              assert(attr - 1 < maxBdrAttr);
              marker[attr - 1] = 1;
            }
            Internal::MFEMVectorCoefficient vc(this->getFiniteElementSpace().getMesh(), fn);
            this->getHandle().ProjectBdrCoefficient(vc, marker);
            return *this;
          }
        }
        else
        {
          assert(false);
          return *this;
        }
      }

      inline GridFunction* copy() const noexcept final override
      {
        return new GridFunction(*this);
      }
  };

  template <class ... Ts>
  GridFunction(const H1<Ts...>&) -> GridFunction<H1<Ts...>>;

}

#include "GridFunction.hpp"

#endif
