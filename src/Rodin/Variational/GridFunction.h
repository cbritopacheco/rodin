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
#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "Rodin/Threads/ThreadPool.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "Component.h"
#include "Restriction.h"
#include "LazyEvaluator.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::FormLanguage
{
  template <class Derived, class FESType>
  struct Traits<Variational::GridFunctionBase<Derived, FESType>>
  {
    using FES = FESType;
    using Element = typename Traits<FES>::Element;
    using RangeType = typename Traits<FES>::RangeType;
  };
}

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
   * This class contains the common routines for the behaviour of a
   * GridFunction object. It provides a common interface for the manipulation
   * of its data and weights, as well as projection utilities and convenience
   * functions.
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
  template <class Derived, class FESType = typename FormLanguage::Traits<Derived>::FES>
  class GridFunctionBase : public LazyEvaluator<GridFunctionBase<Derived, FESType>>
  {
    public:
      using FES = FESType;

      /// Type of finite element
      using Element = typename FormLanguage::Traits<FES>::Element;

      /// Parent class
      using Parent = LazyEvaluator<GridFunctionBase<Derived, FES>>;

      /// Range type of value
      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      static_assert(std::is_same_v<RangeType, Scalar> || std::is_same_v<RangeType, Math::Vector>);

      GridFunctionBase(const FES& fes)
        : Parent(std::cref(*this)),
          m_fes(fes),
          m_data(fes.getVectorDimension(), fes.getSize())
      {
        m_data.setZero();
      }

      GridFunctionBase(const GridFunctionBase& other)
        : Parent(std::cref(*this)),
          m_fes(other.m_fes),
          m_data(other.m_data),
          m_weights(other.m_weights)
      {}

      GridFunctionBase(GridFunctionBase&& other)
        : Parent(std::cref(*this)),
          m_fes(std::move(other.m_fes)),
          m_data(std::move(other.m_data)),
          m_weights(std::move(other.m_weights))
      {}

      GridFunctionBase& operator=(GridFunctionBase&& other)
      {
        m_fes = std::move(other.m_fes);
        m_data = std::move(other.m_data);
        m_weights = std::move(other.m_weights);
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
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        return m_data.maxCoeff();
      }

      inline
      constexpr
      Scalar max(Index& idx) const
      {
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        return m_data.maxCoeff(&idx);
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
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        return m_data.minCoeff();
      }

      inline
      constexpr
      Scalar min(Index& idx) const
      {
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        return m_data.minCoeff(&idx);
      }

      inline
      constexpr
      Index argmax() const
      {
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        Index idx = 0;
        m_data.maxCoeff(&idx);
        return idx;
      }

      inline
      constexpr
      Index argmin() const
      {
        static_assert(std::is_same_v<RangeType, Scalar>,
            "GridFunction must be scalar valued.");
        Index idx = 0;
        m_data.minCoeff(&idx);
        return idx;
      }

      inline
      Derived& normalize()
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>,
            "GridFunction must be vector valued.");
        for (size_t i = 0; i < getSize(); i++)
          getData().col(i).normalize();
        return static_cast<Derived&>(*this);
      }

      inline
      Derived& stableNormalize()
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>,
            "GridFunction must be vector valued.");
        for (size_t i = 0; i < getSize(); i++)
          getData().col(i).stableNormalize();
        return static_cast<Derived&>(*this);
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
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        assert(getFiniteElementSpace().getVectorDimension() >= 1);
        return static_cast<const Derived&>(*this).x();
      }

      inline
      constexpr
      auto y() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        assert(getFiniteElementSpace().getVectorDimension() >= 2);
        return static_cast<const Derived&>(*this).y();
      }

      inline
      constexpr
      auto z() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        assert(getFiniteElementSpace().getVectorDimension() >= 3);
        return static_cast<const Derived&>(*this).z();
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
        static_assert(std::is_same_v<RangeType, Scalar>);
        m_data = m_data.array() + rhs;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Substraction of a scalar value.
       */
      inline
      Derived& operator-=(Scalar rhs)
      {
        static_assert(std::is_same_v<RangeType, Scalar>);
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

      /**
       * @brief Projects a scalar valued function on the region of the mesh
       * with the given attribute.
       * @param[in] fn Scalar valued function
       * @param[in] attr Attribute
       */
      inline
      auto& project(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return project(fn, FlatSet<Geometry::Attribute>{attr});
      }

      /**
       * @brief Projects a scalar valued function on the region of the mesh
       * with the given attributes.
       * @param[in] fn Scalar valued function
       * @param[in] attrs Set of attributes
       */
      inline
      auto& project(std::function<Scalar(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return project(ScalarFunction(fn), attrs);
      }

      template <class NestedDerived>
      inline
      Derived& project(const FunctionBase<NestedDerived>& fn)
      {
        return project(fn, FlatSet<Geometry::Attribute>{});
      }

      /**
       * @brief Projects a FunctionBase instance
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attribute.
       *
       * It is a convenience function to call
       * project(const FunctionBase&, const FlatSet<Geometry::Atribute>&) with one
       * attribute.
       */
      template <class NestedDerived>
      inline
      Derived& project(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return project(fn, FlatSet<Geometry::Attribute>{attr});
      }

      /**
       * @brief Projects a FunctionBase instance on the grid function.
       *
       * This function will project a FunctionBase instance on the
       * domain elements with the given attributes. If the attribute set is
       * empty, this function will project over all elements in the mesh.
       */
      template <class NestedDerived>
      Derived& project(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        using Function = FunctionBase<NestedDerived>;
        using FunctionRangeType = typename FormLanguage::Traits<Function>::RangeType;
        static_assert(std::is_same_v<RangeType, FunctionRangeType>);
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension();

        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
#ifdef RODIN_MULTITHREADED
          auto& threadPool = Threads::getGlobalThreadPool();
          auto loop =
            [&](const Index start, const Index end)
            {
              const size_t capacity = fes.getSize() / threadPool.getThreadCount();
              std::vector<Index> is;
              is.reserve(capacity);
              std::vector<Scalar> vs;
              vs.reserve(capacity);
              for (Index i = start; i < end; ++i)
              {
                const auto it = mesh.getCell(i);
                const auto& polytope = *it;
                if (attrs.size() == 0 || attrs.count(polytope.getAttribute()))
                {
                  const auto& i = polytope.getIndex();
                  const auto& fe = fes.getFiniteElement(d, i);
                  const auto& trans = mesh.getPolytopeTransformation(d, i);
                  for (size_t local = 0; local < fe.getCount(); local++)
                  {
                    const Geometry::Point p(polytope, trans, fe.getNode(local));
                    assert(m_data.rows() == 1);
                    is.push_back(fes.getGlobalIndex({ d, i }, local));
                    vs.push_back(fn.getValue(p));
                  }
                }
              }
              assert(is.size() == vs.size());
              m_mutex.lock();
              for (Index i = 0; i < is.size(); i++)
                m_data(is[i]) = vs[i];
              m_mutex.unlock();
            };
          threadPool.pushLoop(0, mesh.getCellCount(), loop);
          threadPool.waitForTasks();
#else
          for (auto it = mesh.getCell(); !it.end(); ++it)
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
                assert(m_data.rows() == 1);
                m_data(fes.getGlobalIndex({ d, i }, local)) = fn.getValue(p);
              }
            }
          }
#endif
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          Math::Vector value;
          for (auto it = mesh.getCell(); !it.end(); ++it)
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
                fn.getValue(value, p);
                m_data.col(fes.getGlobalIndex({ d, i }, local)) = value;
              }
            }
          }
        }
        else
        {
          assert(false);
        }
        return static_cast<Derived&>(*this);
      }

      inline
      auto& projectOnBoundary(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      }

      inline
      auto& projectOnBoundary(std::function<Scalar(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return projectOnBoundary(ScalarFunction(fn), attrs);
      }

      template <class NestedDerived>
      inline
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn)
      {
        return static_cast<Derived&>(*this).projectOnBoundary(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      inline
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnBoundary(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
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

      inline
      auto& projectOnFaces(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      inline
      auto& projectOnFaces(std::function<Scalar(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return projectOnFaces(ScalarFunction(fn), attrs);
      }

      template <class NestedDerived>
      inline
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn)
      {
        return static_cast<Derived&>(*this).projectOnFaces(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      inline
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnFaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnFaces(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        using Function = FunctionBase<NestedDerived>;
        using FunctionRangeType = typename FormLanguage::Traits<Function>::RangeType;
        static_assert(std::is_same_v<RangeType, FunctionRangeType>);
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension() - 1;
        for (auto it = mesh.getFace(); !it.end(); ++it)
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
      auto& projectOnInterfaces(std::function<Scalar(const Geometry::Point&)> fn, Geometry::Attribute attr)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      inline
      auto& projectOnInterfaces(std::function<Scalar(const Geometry::Point&)> fn, const FlatSet<Geometry::Attribute>& attrs = {})
      {
        return projectOnInterfaces(ScalarFunction(fn), attrs);
      }

      template <class NestedDerived>
      inline
      Derived& projectOnInterfaces(const FunctionBase<NestedDerived>& fn)
      {
        return static_cast<Derived&>(*this).projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{});
      }

      template <class NestedDerived>
      inline
      Derived& projectOnInterfaces(const FunctionBase<NestedDerived>& fn, Geometry::Attribute attr)
      {
        return projectOnInterfaces(fn, FlatSet<Geometry::Attribute>{attr});
      }

      template <class NestedDerived>
      Derived& projectOnInterfaces(const FunctionBase<NestedDerived>& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        using Function = FunctionBase<NestedDerived>;
        using FunctionRangeType = typename FormLanguage::Traits<Function>::RangeType;
        static_assert(std::is_same_v<RangeType, FunctionRangeType>);
        const auto& fes = getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension() - 1;
        for (auto it = mesh.getInterface(); !it.end(); ++it)
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

      /**
       * @brief Computes the weights from the data.
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      Derived& setWeights()
      {
        return static_cast<Derived&>(*this).setWeights();
      }

      /**
       * @brief Sets the weights in the GridFunction object and computes the
       * values at all the degrees of freedom.
       * @note CRTP function to be overriden in Derived class.
       */
      template <class Vector>
      inline
      Derived& setWeights(Vector&& weights)
      {
        return static_cast<Derived&>(*this).setWeights(std::forward<Vector>(weights));
      }

      /**
       * @brief Sets the weights and data in the GridFunction object. No
       * computation is performed.
       */
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

      template <class Value>
      inline
      Derived& setValue(const std::pair<size_t, Index>& p, size_t local, Value&& v)
      {
        return setValue(getFiniteElementSpace().getGlobalIndex(p, local), std::forward<Value>(v));
      }

      template <class Value>
      inline
      Derived& setValue(Index global, Value&& v)
      {
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(m_data.size() >= 0);
          assert(global < static_cast<size_t>(m_data.size()));
          m_data.coeffRef(global) = std::forward<Value>(v);
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector>)
        {
          assert(m_data.cols() >= 0);
          assert(global < static_cast<size_t>(m_data.cols()));
          m_data.col(global) = std::forward<Value>(v);
        }
        else
        {
          assert(false);
        }
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the value at the given polytope on the local degree of
       * freedom.
       */
      inline
      auto getValue(const std::pair<size_t, Index>& p, size_t local) const
      {
        return getValue(getFiniteElementSpace().getGlobalIndex(p, local));
      }


      /**
       * @brief Gets the value of the GridFunction at the global degree of
       * freedom index.
       */
      inline
      auto getValue(Index global) const
      {
        if constexpr (std::is_same_v<RangeType, Scalar>)
        {
          assert(m_data.size() >= 0);
          assert(global < static_cast<size_t>(m_data.size()));
          return m_data.coeff(global);
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

      /**
       * @brief Gets the interpolated value at the point.
       */
      inline
      auto getValue(const Geometry::Point& p) const
      {
        RangeType out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& fes = m_fes.get();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          interpolate(out, p);
        }
        else
        {
          if (polytopeMesh.isSubMesh())
          {
            const auto& submesh = polytopeMesh.asSubMesh();
            assert(submesh.getParent() == fes.getMesh());
            interpolate(out, submesh.inclusion(p));
          }
          else if (fesMesh.isSubMesh())
          {
            const auto& submesh = fesMesh.asSubMesh();
            assert(submesh.getParent() == polytopeMesh);
            interpolate(out, submesh.restriction(p));
          }
          else
          {
            assert(false);
            if constexpr (std::is_same_v<RangeType, Scalar>)
              out = NAN;
            else if constexpr (std::is_same_v<RangeType, Math::Vector>)
              out.setConstant(NAN);
          }
        }
        return out;
      }

      inline
      constexpr
      void getValue(Math::Vector& res, const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        interpolate(res, p);
      }

      /**
       * @brief Interpolates the GridFunction at the given point.
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      constexpr
      Scalar interpolate(const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Scalar>);
        return static_cast<const Derived&>(*this).interpolate(p);
      }

      /**
       * @brief Interpolates the GridFunction at the given point.
       * @note CRTP function to be overriden in Derived class.
       */
      inline
      constexpr
      void interpolate(Math::Vector& res, const Geometry::Point& p) const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector>);
        static_cast<const Derived&>(*this).interpolate(res, p);
      }

      /**
       * @brief Deleted function.
       */
      inline
      constexpr
      void getValue(Math::Matrix& res, const Geometry::Point& p) const = delete;

    private:
      void interpolate(Scalar& res, const Geometry::Point& p) const
      {
        res = interpolate(p);
      }

      std::reference_wrapper<const FES> m_fes;
      Math::Matrix m_data;
      std::optional<Math::Vector> m_weights;
      mutable Threads::Mutex m_mutex;
  };
}

#endif
