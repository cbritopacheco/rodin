/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_GRIDFUNCTION_H

#include <utility>
#include <fstream>
#include <functional>
#include <boost/filesystem.hpp>
#include <type_traits>

#include "Rodin/Math.h"

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/SubMesh.h"

#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/MEDIT.h"
#include "Rodin/IO/EnSight6.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "Component.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::FormLanguage
{
  template <class Derived, class FES, class Data>
  struct Traits<Variational::GridFunctionBase<Derived, FES, Data>>
  {
    using FESType = FES;
    using DataType = Data;
  };

  template <class FES, class Data>
  struct Traits<Variational::GridFunction<FES, Data>>
  {
    using FESType = FES;
    using DataType = Data;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GridFunctionSpecializations GridFunction Template Specializations
   * @brief Template specializations of the GridFunction class.
   * @see GridFunction
   */

  template <class StrictType>
  class GridFunctionBaseReference
    : public FunctionBase<GridFunctionBaseReference<StrictType>>
  {
    public:
      using Parent = FunctionBase<GridFunctionBaseReference<StrictType>>;

      /**
       * @brief R-Values are not allowed.
       */
      GridFunctionBaseReference(StrictType&&) = delete;

      /**
       * @brief Prevent implicit copies.
       */
      GridFunctionBaseReference(const StrictType& ref) = delete;

      /**
       * @brief Constructs the LazyEvaluator object from a constant reference
       * the data-full object.
       */
      explicit
      constexpr
      GridFunctionBaseReference(std::reference_wrapper<const StrictType> ref)
        : m_ref(ref)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      GridFunctionBaseReference(const GridFunctionBaseReference& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr
      GridFunctionBaseReference(GridFunctionBaseReference&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      GridFunctionBaseReference& operator=(const GridFunctionBaseReference&) = delete;

      GridFunctionBaseReference& operator=(GridFunctionBaseReference&&) = delete;

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return m_ref.get().getValue(p);
      }

      template <class T>
      constexpr
      void getValue(T& res, const Geometry::Point& p) const
      {
        m_ref.get().getValue(res, p);
      }

      GridFunctionBaseReference* copy() const noexcept final override
      {
        return new GridFunctionBaseReference(*this);
      }

    private:
      std::reference_wrapper<const StrictType> m_ref;
  };

  /**
   * @brief Abstract base class for GridFunction objects.
   *
   * This class contains the common routines for the behaviour of a
   * GridFunction object. It provides a common interface for the manipulation
   * of its data and weights, as well as projection utilities and convenience
   * functions.
   */
  template <
    class Derived,
    class FES = typename FormLanguage::Traits<Derived>::FESType,
    class Data = typename FormLanguage::Traits<Derived>::DataType>
  class GridFunctionBase
    : public GridFunctionBaseReference<GridFunctionBase<Derived, FES, Data>>
  {
    public:
      using FESType = FES;

      using DataType = Data;

      /// Range type of value
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Local>;

      /// Represents the Context of the P1 space
      using ContextType = Context::Local;

      /// Type of finite element
      using ElementType = typename FormLanguage::Traits<FESType>::ElementType;

      /// Parent class
      using Parent =
        GridFunctionBaseReference<GridFunctionBase<Derived, FESType, Data>>;

      static_assert(
          std::is_same_v<RangeType, ScalarType> ||
          std::is_same_v<RangeType, Math::Vector<ScalarType>>);

      GridFunctionBase(const FES& fes)
        : Parent(std::cref(*this)),
          m_fes(std::cref(fes))
      {}

      GridFunctionBase(const GridFunctionBase& other)
        : Parent(std::cref(*this)),
          m_fes(other.m_fes)
      {}

      GridFunctionBase(GridFunctionBase&& other)
        : Parent(std::cref(*this)),
          m_fes(std::move(other.m_fes))
      {}

      virtual ~GridFunctionBase() = default;

      GridFunctionBase& operator=(GridFunctionBase&& other)
      {
        m_fes = std::move(other.m_fes);
        return *this;
      }

      GridFunctionBase& operator=(const GridFunctionBase&) = delete;

      constexpr
      auto x() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(static_cast<Derived&>(*this), 0);
      }

      constexpr
      auto y() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(static_cast<Derived&>(*this), 1);
      }

      constexpr
      auto z() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(static_cast<Derived&>(*this), 2);
      }

      constexpr
      Derived& setData(const DataType& data, size_t offset = 0)
      {
        return static_cast<Derived&>(*this).setData(data, offset);
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      constexpr
      auto& getData()
      {
        return static_cast<Derived&>(*this).getData();
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      constexpr
      const DataType& getData() const
      {
        return static_cast<const Derived&>(*this).getData();
      }

      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      constexpr
      size_t getSize() const
      {
        return getFiniteElementSpace().getSize();
      }

      constexpr
      size_t getDimension() const
      {
        return getFiniteElementSpace().getVectorDimension();
      }

      Derived& load(
          const boost::filesystem::path& filename,
          IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        std::ifstream input(filename.c_str());
        if (!input)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open input file stream." << Alert::NewLine
            << "Filename: \"" << filename << "\"" << Alert::NewLine
            << "Please check if the file exists and is accessible."
            << Alert::Raise;
        }

        switch (fmt)
        {
          case IO::FileFormat::MFEM:
          {
            IO::GridFunctionLoader<IO::FileFormat::MFEM, FESType, DataType>(
              static_cast<Derived&>(*this)).load(input);
            break;
          }
          case IO::FileFormat::MEDIT:
          {
            IO::GridFunctionLoader<IO::FileFormat::MEDIT, FES, DataType>(
              static_cast<Derived&>(*this)).load(input);
            break;
          }
          default:
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Unsupported file format for loading GridFunction." << Alert::NewLine
              << "Format: \"" << fmt << "\""
              << Alert::Raise;
          }
        }
        return static_cast<Derived&>(*this);
      }

      void save(
          const boost::filesystem::path& filename,
          IO::FileFormat fmt = IO::FileFormat::MFEM) const
      {
        std::ofstream output(filename.c_str());
        if (!output)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open output file stream." << Alert::NewLine
            << "Filename: \"" << filename << "\"" << Alert::NewLine
            << "Please check if the path is valid and writable."
            << Alert::Raise;
        }

        switch (fmt)
        {
          case IO::FileFormat::MFEM:
          {
            IO::GridFunctionPrinter<IO::FileFormat::MFEM, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(output);
            break;
          }
          case IO::FileFormat::MEDIT:
          {
            IO::GridFunctionPrinter<IO::FileFormat::MEDIT, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(output);
            break;
          }
          case IO::FileFormat::ENSIGHT6:
          {
            IO::GridFunctionPrinter<IO::FileFormat::ENSIGHT6, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(output);
            break;
          }
          default:
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Unsupported file format for saving GridFunction." << Alert::NewLine
              << "Format: \"" << fmt << "\""
              << Alert::Raise;
          }
        }
        output.close();
      }

      constexpr
      RangeType getValue(const Geometry::Point& p) const
      {
        RangeType res;
        getValue(res, p);
        return res;
      }

      /**
       * @brief Gets the interpolated value at the point.
       */
      constexpr
      void getValue(RangeType& res, const Geometry::Point& p) const
      {
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& fes = m_fes.get();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          static_cast<const Derived&>(*this).interpolate(res, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          static_cast<const Derived&>(*this).interpolate(res, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          if (restriction)
          {
            static_cast<const Derived&>(*this).interpolate(res, *restriction);
          }
          else
          {
            assert(false);
          }
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @brief Interpolates the GridFunction at the given point.
       *
       * @note Can be overriden.
       */
      constexpr
      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(res, p);
      }

      template <class NestedDerived>
      void project(const FunctionBase<NestedDerived>& fn, const std::pair<size_t, Index>& p)
      {
        static_cast<Derived&>(*this).project(fn, p);
      }

      template <class NestedDerived>
      Derived& operator=(const FunctionBase<NestedDerived>& fn)
      {
        return projectOnCells(fn);
      }

      Derived& operator=(std::function<RangeType(const Geometry::Point&)> fn)
      {
        return projectOnCells(fn);
      }

      Derived& operator=(std::function<void(RangeType&, const Geometry::Point&)> fn)
      {
        return projectOnCells(fn);
      }

      Derived& operator=(const RangeType& v)
      {
        return projectOnCells([&](RangeType& res, const Geometry::Point&) { res = v; });
      }

      /**
       * @brief Projects a scalar valued function on the region of the mesh
       * with the given attribute.
       * @param[in] fn Scalar valued function
       * @param[in] attr Attribute
       */
      auto& projectOnCells(
          std::function<RangeType(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnCells(fn, FlatSet<Geometry::Attribute>{ attr });
      }

      auto& projectOnCells(
          std::function<void(RangeType&, const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnCells(fn, FlatSet<Geometry::Attribute>{ attr });
      }

      auto& projectOnCells(
          std::function<RangeType(const Geometry::Point&)> fn,
          const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return projectOnCells(Function(fn), attrs);
      }

      auto& projectOnCells(
          std::function<void(RangeType&, const Geometry::Point&)> fn,
          const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return projectOnCells(Function(fn), attrs);
      }

      template <class NestedDerived>
      Derived& projectOnCells(const FunctionBase<NestedDerived>& fn)
      {
        return projectOnCells(fn, FlatSet<Geometry::Attribute>{});
      }

      /**
       * @brief Projects a FunctionBase instance
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attribute.
       *
       * It is a convenience function to call
       * projectOnCells(const FunctionBase&, const FlatSet<Geometry::Atribute>&) with one
       * attribute.
       */
      template <class NestedDerived>
      Derived& projectOnCells(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnCells(fn, FlatSet<Geometry::Attribute>{attr});
      }

      /**
       * @brief Projects a FunctionBase instance on the grid function.
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attributes. If the attribute set is
       * empty, this function will project over all elements in the mesh.
       */
      template <class NestedDerived>
      Derived& projectOnCells(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        return static_cast<Derived&>(*this).projectOnCells(fn, attrs);
      }

      auto& projectOnBoundary(
          std::function<RangeType(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      }

      auto& projectOnBoundary(
          std::function<void(RangeType&, const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      }

      auto& projectOnBoundary(
          std::function<RangeType(const Geometry::Point&)> fn,
          const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnBoundary(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnBoundary(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      auto& projectOnBoundary(
          std::function<void(RangeType&, const Geometry::Point&)> fn,
          const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnBoundary(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnBoundary(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      template <class NestedDerived>
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        for (auto it = mesh.getBoundary(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
          {
            const auto& polytope = *it;
            if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
              project(fn, { polytope.getDimension(), polytope.getIndex() });
          }
        }
        return static_cast<Derived&>(*this);
      }

      auto& projectOnFaces(
          std::function<RangeType(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      auto& projectOnFaces(
          std::function<void(RangeType&, const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      auto& projectOnFaces(
          std::function<RangeType(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnFaces(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnFaces(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      auto& projectOnFaces(
          std::function<void(RangeType&, const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnFaces(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnFaces(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      template <class NestedDerived>
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        for (auto it = mesh.getFace(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
            project(fn, { polytope.getDimension(), polytope.getIndex() });
        }
        return static_cast<Derived&>(*this);
      }

      auto& projectOnInterfaces(
          std::function<RangeType(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      auto& projectOnInterfaces(
          std::function<void(RangeType&, const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{ attr });
      }

      auto& projectOnInterfaces(
          std::function<RangeType(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnInterfaces(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnInterfaces(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      auto& projectOnInterfaces(
          std::function<void(RangeType&, const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return projectOnInterfaces(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
        {
          return projectOnInterfaces(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
      }

      template <class NestedDerived>
      Derived& projectOnInterfaces(const FunctionBase<NestedDerived>& fn)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      Derived& projectOnInterfaces(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnInterfaces(
          const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        for (auto it = mesh.getInterface(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
            project(fn, { polytope.getDimension(), polytope.getIndex() });
        }
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Searches the minimum value in the grid function data.
       * @returns Minimum value in grid function.
       *
       * This function will compute the minimum value in the grid function
       * data array.
       *
       * @section Complexity
       * The operation is linear in the size of the number of entries in the
       * underlying matrix.
       */
      constexpr
      ScalarType min() const
      {
        Index _unused;
        return static_cast<const Derived&>(*this).min(_unused);
      }

      /**
       * @brief Searches for the maximum value in the grid function data.
       * @returns Maximum value in grid function.
       *
       * This function will compute the maximum value in the grid function
       * data array.
       *
       * @section Complexity
       * The operation is linear in the size of the number of entries in the
       * underlying matrix.
       */
      constexpr
      ScalarType max() const
      {
        Index _unused;
        return static_cast<const Derived&>(*this).max(_unused);
      }

      constexpr
      Index argmin() const
      {
        Index idx = 0;
        static_cast<const Derived&>(*this).min(idx);
        return idx;
      }

      constexpr
      Index argmax() const
      {
        Index idx;
        static_cast<const Derived&>(*this).max(idx);
        return idx;
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      constexpr
      ScalarType min(Index& idx) const
      {
        return static_cast<const Derived&>(*this).min(idx);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      constexpr
      ScalarType max(Index& idx) const
      {
        return static_cast<const Derived&>(*this).max(idx);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      ScalarType& operator[](Index global)
      {
        return static_cast<Derived&>(*this).operator[](global);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      const ScalarType& operator[](Index global) const
      {
        return static_cast<const Derived&>(*this).operator[](global);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator+=(const ScalarType& rhs)
      {
        return static_cast<Derived&>(*this).operator+=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator-=(const ScalarType& rhs)
      {
        return static_cast<Derived&>(*this).operator-=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator*=(const ScalarType& rhs)
      {
        return static_cast<Derived&>(*this).operator*=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator/=(const ScalarType& rhs)
      {
        return static_cast<Derived&>(*this).operator/=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator+=(const GridFunctionBase& rhs)
      {
        return static_cast<Derived&>(*this).operator+=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator-=(const GridFunctionBase& rhs)
      {
        return static_cast<Derived&>(*this).operator-=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator*=(const GridFunctionBase& rhs)
      {
        return static_cast<Derived&>(*this).operator*=(rhs);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      Derived& operator/=(const GridFunctionBase& rhs)
      {
        return static_cast<Derived&>(*this).operator/=(rhs);
      }

    private:
      std::reference_wrapper<const FESType> m_fes;
  };

  template <class FES>
  class GridFunction<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>> final
    : public GridFunctionBase<
        GridFunction<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
  {
    public:
      using FESType = FES;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DataType = Math::Vector<ScalarType>;

      using Parent = GridFunctionBase<GridFunction<FESType, DataType>>;

      using Parent::operator=;
      using Parent::min;
      using Parent::max;
      using Parent::projectOnCells;

      GridFunction(const FESType& fes)
        : Parent(fes)
      {
        auto& data = this->getData();
        data.resize(fes.getSize());
        data.setZero();
      }

      GridFunction(const GridFunction& other)
        : Parent(other)
      {}

      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
      {}

      GridFunction& operator=(GridFunction&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      virtual ~GridFunction() = default;

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

      /**
       * @brief Interpolates the GridFunction at the given point.
       *
       * @note Can be overriden.
       */
      constexpr
      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        RangeType v;
        for (Index local = 0; local < count; ++local)
        {
          const auto mapping = fes.getInverseMapping({ d, i }, fe.getBasis(local));
          mapping(v, p);
          const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * v;
          if (local == 0)
            res = k; // Initializes the result (resizes)
          else
            res += k; // Accumulates the result (does not resize)
        }
      }

      template <class NestedDerived>
      GridFunction& projectOnCells(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        for (auto it = mesh.getCell(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
            project(fn, { polytope.getDimension(), polytope.getIndex() });
        }
        return *this;
      }

      template <class NestedDerived>
      void project(const FunctionBase<NestedDerived>& fn, const std::pair<size_t, Index>& p)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& [d, i] = p;
        const auto& fe = fes.getFiniteElement(d, i);
        const auto mapping =
          fes.getMapping({ d, i }, fn.template cast<RangeType>());
        for (Index local = 0; local < fe.getCount(); local++)
        {
          const Index global = fes.getGlobalIndex({ d, i }, local);
          this->operator[](global) = fe.getLinearForm(local)(mapping);
        }
      }

      GridFunction& setData(const DataType& data, size_t offset = 0)
      {
        const auto sz = this->getFiniteElementSpace().getSize();
        assert(offset + sz <= data.size());
        this->getData() = data.segment(offset, sz);
        return *this;
      }

      constexpr
      auto& getData()
      {
        return m_data;
      }

      constexpr
      const DataType& getData() const
      {
        return m_data;
      }

    private:
      DataType m_data;
  };

  template <class FES>
  GridFunction(const FES& fes)
    -> GridFunction<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>;

  template <class FES, class Data>
  GridFunction(const FES& fes, Data&& data)
    -> GridFunction<FES, Data>;
}

#endif
