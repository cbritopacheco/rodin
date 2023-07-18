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

#include "Rodin/Math.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/MEDIT.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "Component.h"
#include "Restriction.h"
#include "LazyEvaluator.h"
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

  /**
   * @brief Abstract base class for GridFunction objects.
   *
   * @section gridfunction-data-layout Data layout
   *
   * The data of the GridFunctionBase object can be accessed via a call to @ref
   * getData(). The i-th column of the returned Math::Matrix object corresponds
   * to the value of the grid function at the global i-th degree of freedom in
   * the finite element space. Furthermore, the following conditions are
   * satisfied:
   * ```
   *  const Math::Matrix& data = gf.getData();
   *  assert(data.rows() == gf.getFiniteElementSpace().getVectorDimension());
   *  assert(data.cols() == gf.getFiniteElementSpace().getSize());
   * ```
   */
  template <class Derived, class FES>
  class GridFunctionBase : public LazyEvaluator<GridFunctionBase<Derived, FES>>
  {
    public:
      /// Type of finite element
      using Element = typename FES::Element;

      /// Parent class
      using Parent = LazyEvaluator<GridFunctionBase<Derived, FES>>;

      /// Range type of value
      using RangeType = typename FES::RangeType;

      static_assert(std::is_same_v<RangeType, Scalar> || std::is_same_v<RangeType, Math::Vector>);

      GridFunctionBase(const FES& fes)
        : Parent(*this),
          m_fes(fes),
          m_data(fes.getVectorDimension(), fes.getSize())
      {
        m_data.setZero();
      }

      GridFunctionBase(const GridFunctionBase& other)
        : Parent(*this),
          m_fes(other.m_fes),
          m_data(other.m_data)
      {}

      GridFunctionBase(GridFunctionBase&& other)
        : Parent(*this),
          m_fes(std::move(other.m_fes)),
          m_data(std::move(other.m_data))
      {}

      GridFunctionBase& operator=(GridFunctionBase&& other)
      {
        m_fes = std::move(other.m_fes);
        m_data = std::move(other.m_data);
        return *this;
      }

      GridFunctionBase& operator=(const GridFunctionBase&) = delete;

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
      inline
      constexpr
      Scalar max() const
      {
        return m_data.maxCoeff();
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
      inline
      constexpr
      Scalar min() const
      {
        return m_data.minCoeff();
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
      size_t getSize() const
      {
        return getFiniteElementSpace().getSize();
      }


      inline
      Derived& setZero()
      {
        m_data.setZero();
        if (m_weights)
          m_weights->setZero();
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Bulk assigns the value to the whole data array.
       */
      inline
      Derived& operator=(Scalar v)
      {
        static_assert(std::is_same_v<RangeType, Scalar>);
        m_data.setConstant(v);
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator=(const Math::Vector& v)
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        Math::Matrix& data = m_data;
        assert(data.cols() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(data.cols()); i++)
          data.col(i) = v;
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator=(std::function<RangeType(const Geometry::Point&)> fn)
      {
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(getFiniteElementSpace().getVectorDimension() == 1);
          return project(ScalarFunction(fn));
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          return project(VectorFunction(getFiniteElementSpace().getVectorDimension(), fn));
        }
        else
        {
          assert(false);
          return static_cast<Derived&>(*this);
        }
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
        m_data = m_data.array() + rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Substraction of a scalar value.
       */
      inline
      Derived& operator-=(Scalar rhs)
      {
        m_data = m_data.array() - rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Multiplication by a scalar value.
       */
      inline
      Derived& operator*=(Scalar rhs)
      {
        m_data = m_data.array() * rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Division by a scalar value.
       */
      inline
      Derived& operator/=(Scalar rhs)
      {
        m_data = m_data.array() / rhs;
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
          m_data = m_data.array() + rhs.m_data.array();
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
          m_data = m_data.array() - rhs.m_data.array();
        }
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& operator*=(const GridFunctionBase& rhs)
      {
        if (this == &rhs)
        {
          m_data = m_data.array() * m_data.array();
        }
        else
        {
          assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
          m_data = m_data.array() * rhs.m_data.array();
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
          m_data = m_data.array() / rhs.m_data.array();
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

      template <class NestedDerived>
      inline
      Derived& project(const FunctionBase<NestedDerived>& fn)
      {
        return project(fn, std::set<Geometry::Attribute>{});
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
      inline
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
      Derived& project(const FunctionBase<NestedDerived>& fn, const std::set<Geometry::Attribute>& attrs)
      {
        using Function = FunctionBase<NestedDerived>;
        using FunctionRangeType = typename FormLanguage::Traits<Function>::RangeType;
        static_assert(std::is_same_v<RangeType, FunctionRangeType>);
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension();
        for (auto it = mesh.getElement(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
          {
            const auto& i = polytope.getIndex();
            const auto& fe = fes.getFiniteElement(d, i);
            const auto& trans = mesh.getPolytopeTransformation(d, i);
            for (size_t local = 0; local < fe.getCount(); local++)
            {
              const Geometry::Point p(polytope, trans, fe.getNode(local));
              if constexpr (std::is_same_v<RangeType, Scalar>)
              {
                assert(m_data.rows() == 1);
                m_data(fes.getGlobalIndex({ d, i }, local)) = fn.getValue(p);
              }
              else if constexpr (std::is_same_v<RangeType, Math::Vector>)
              {
                m_data.col(fes.getGlobalIndex({ d, i }, local)) = fn.getValue(p);
              }
              else
              {
                assert(false);
              }
            }
          }
        }
        return static_cast<Derived&>(*this);
      }

      inline
      auto& projectOnBoundary(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, std::set<Geometry::Attribute>{attr});
      }

      inline
      auto& projectOnBoundary(std::function<Scalar(const Geometry::Point&)> fn, const std::set<Geometry::Attribute>& attrs = {})
      {
        return projectOnBoundary(ScalarFunction(fn), attrs);
      }

      template <class NestedDerived>
      inline
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn)
      {
        return static_cast<Derived&>(*this).projectOnBoundary(fn, std::set<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      inline
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, std::set<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, const std::set<Geometry::Attribute>& attrs)
      {
        using Function = FunctionBase<NestedDerived>;
        using FunctionRangeType = typename FormLanguage::Traits<Function>::RangeType;
        static_assert(std::is_same_v<RangeType, FunctionRangeType>);
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension() - 1;
        for (auto it = mesh.getBoundary(); !it.end(); ++it)
        {
          const auto& polytope = *it;
          if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
          {
            const auto& i = polytope.getIndex();
            const auto& fe = fes.getFiniteElement(d, i);
            const auto& trans = mesh.getPolytopeTransformation(d, i);
            for (size_t local = 0; local < fe.getCount(); local++)
            {
              const Geometry::Point p(polytope, trans, fe.getNode(local));
              if constexpr (std::is_same_v<RangeType, Scalar>)
              {
                assert(m_data.rows() == 1);
                m_data(fes.getGlobalIndex({ d, i }, local)) = fn.getValue(p);
              }
              else if constexpr (std::is_same_v<RangeType, Math::Vector>)
              {
                m_data.col(fes.getGlobalIndex({ d, i }, local)) = fn.getValue(p);
              }
              else
              {
                assert(false);
              }
            }
          }
        }
        return static_cast<Derived&>(*this);
      }

      Derived& load(
          const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM)
      {
        std::ifstream input(filename.c_str());
        if (!input)
        {
          Alert::Exception()
            << "Failed to open " << filename << " for reading."
            << Alert::Raise;
        }

        switch (fmt)
        {
          case IO::FileFormat::MFEM:
          {
            IO::GridFunctionLoader<IO::FileFormat::MFEM, FES> loader(static_cast<Derived&>(*this));
            loader.load(input);
            break;
          }
          case IO::FileFormat::MEDIT:
          {
            IO::GridFunctionLoader<IO::FileFormat::MEDIT, FES> loader(static_cast<Derived&>(*this));
            loader.load(input);
            break;
          }
          default:
          {
            Alert::Exception()
              << "Loading from \"" << fmt << "\" format unsupported."
              << Alert::Raise;
          }
        }
        return static_cast<Derived&>(*this);
      }

      void save(
          const boost::filesystem::path& filename, IO::FileFormat fmt = IO::FileFormat::MFEM,
          size_t precision = RODIN_DEFAULT_GRIDFUNCTION_SAVE_PRECISION) const
      {
        std::ofstream output(filename.c_str());
        if (!output)
        {
          Alert::Exception()
            << "Failed to open " << filename << " for writing."
            << Alert::Raise;
        }

        output.precision(precision);
        switch (fmt)
        {
          case IO::FileFormat::MFEM:
          {
            IO::GridFunctionPrinter<IO::FileFormat::MFEM, FES> printer(static_cast<const Derived&>(*this));
            printer.print(output);
            break;
          }
          case IO::FileFormat::MEDIT:
          {
            IO::GridFunctionPrinter<IO::FileFormat::MEDIT, FES> printer(static_cast<const Derived&>(*this));
            printer.print(output);
            break;
          }
          default:
          {
            Alert::Exception()
              << "Saving to \"" << fmt << "\" format unsupported."
              << Alert::Raise;
          }
        }
        output.close();
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      template <class Matrix>
      inline
      constexpr
      Derived& setData(Matrix&& data) const
      {
        m_data = std::forward<Matrix>(data);
        return static_cast<Derived&>(*this).setData();
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      inline
      constexpr
      Math::Matrix& getData()
      {
        return m_data;
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      inline
      constexpr
      const Math::Matrix& getData() const
      {
        return m_data;
      }

      inline
      constexpr
      std::optional<Math::Vector>& getWeights()
      {
        return m_weights;
      }

      inline
      constexpr
      const std::optional<Math::Vector>& getWeights() const
      {
        return m_weights;
      }

      inline
      Derived& setWeights()
      {
        return static_cast<Derived&>(*this).setWeights();
      }

      template <class Vector>
      inline
      Derived& setWeights(Vector&& weights)
      {
        return static_cast<Derived&>(*this).setWeights(std::forward<Vector>(weights));
      }

      template <class Vector, class Matrix>
      inline
      Derived& setWeightsAndData(Vector&& weights, Matrix&& data)
      {
        m_weights = std::forward<Vector>(weights);
        m_data = std::forward<Matrix>(data);
        return static_cast<Derived&>(*this);
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getFiniteElementSpace().getVectorDimension(), 1 };
      }

      inline
      auto getValue(Index global) const
      {
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(m_data.size() >= 0);
          assert(global < static_cast<size_t>(m_data.size()));
          return m_data(global);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          assert(m_data.cols() >= 0);
          assert(global < static_cast<size_t>(m_data.cols()));
          return m_data.col(global);
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline
      auto getValue(const std::pair<size_t, Index>& p, size_t local) const
      {
        return getValue(getFiniteElementSpace().getGlobalIndex(p, local));
      }

      /**
       * @brief Gets the interpolated value at the point.
       */
      inline
      auto getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

    private:
      std::reference_wrapper<const FES> m_fes;
      Math::Matrix m_data;
      std::optional<Math::Vector> m_weights;
  };
}

#include "GridFunction.hpp"

#endif
