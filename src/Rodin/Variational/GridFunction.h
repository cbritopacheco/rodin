/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file GridFunction.h
 * @brief Grid function class for representing FEM solutions.
 *
 * This file defines the GridFunction class, which represents functions defined
 * on a finite element mesh by their degrees of freedom. Grid functions are the
 * discrete representation of solutions in finite element analysis.
 *
 * ## Mathematical Foundation
 * A grid function represents a function @f$ u_h \in V_h @f$ by its coefficients:
 * @f[
 *   u_h(x) = \sum_{i=1}^N u_i \phi_i(x)
 * @f]
 * where:
 * - @f$ u_i @f$ are the degrees of freedom (stored in the grid function)
 * - @f$ \phi_i @f$ are the basis functions from the finite element space
 * - @f$ N @f$ is the number of DOFs
 *
 * ## Features
 * - **Storage**: Manages coefficient vector for FEM solutions
 * - **Evaluation**: Point-wise evaluation using basis function interpolation
 * - **I/O**: Export to visualization formats (XDMF, HDF5, MEDIT, etc.)
 * - **Operations**: Arithmetic operations, norms, projections
 * - **Assignment**: Can be set from functions or expressions
 *
 * ## Usage Examples
 * ```cpp
 * P1 Vh(mesh);
 * GridFunction<P1> u(Vh);  // Create grid function
 * 
 * // Set from analytical function
 * u = [](const Point& p) { return sin(p.x()) * cos(p.y()); };
 * 
 * // Evaluate at a point
 * Real value = u(point);
 * 
 * // Export for visualization
 * u.save("solution.vtu");
 * ```
 *
 * @see TrialFunction, FiniteElementSpace
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_GRIDFUNCTION_H

#include <utility>
#include <fstream>
#include <functional>
#include <boost/filesystem.hpp>
#include <type_traits>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeIterator.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/MEDIT.h"
#include "Rodin/IO/MFEM.h"
#include "Rodin/IO/HDF5.h"

#include "ForwardDecls.h"

#include "Function.h"


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

  /**
   * @ingroup RodinVariational
   * @brief Base class for discrete finite element functions.
   *
   * GridFunctionBase represents discrete functions defined on finite element
   * spaces, where the function is represented as a linear combination of basis
   * functions with scalar coefficients stored in a data vector.
   *
   * ## Mathematical Foundation
   * A grid function represents a discrete finite element approximation:
   * @f[
   *   u_h(x) = \sum_{i=1}^N u_i \phi_i(x)
   * @f]
   * where:
   * - @f$ u_h @f$ is the discrete function
   * - @f$ u_i @f$ are the degrees of freedom (coefficients)
   * - @f$ \phi_i @f$ are the finite element basis functions
   * - @f$ N @f$ is the number of degrees of freedom
   *
   * ## Key Features
   * - **DOF Management**: Automatic handling of degrees of freedom storage
   * - **Function Evaluation**: Point-wise evaluation via finite element interpolation
   * - **I/O Support**: Export to various visualization formats (MEDIT, MFEM)
   * - **Space Association**: Strong association with underlying finite element space
   */

  template <class Derived>
  class GridFunctionBaseReference
    : public FunctionBase<GridFunctionBaseReference<Derived>>
  {
    public:
      using Parent = FunctionBase<GridFunctionBaseReference<Derived>>;

      /**
       * @brief R-Values are not allowed.
       */
      GridFunctionBaseReference(Derived&&) = delete;

      /**
       * @brief Prevent implicit copies.
       */
      GridFunctionBaseReference(const Derived& ref) = delete;

      /**
       * @brief Constructs the LazyEvaluator object from a constant reference
       * the data-full object.
       */
      explicit
      constexpr
      GridFunctionBaseReference(std::reference_wrapper<const Derived> ref)
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
      decltype(auto) operator()(const Geometry::Point& p) const
      {
        return m_ref.get().getValue(p);
      }

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return m_ref.get().getValue(p);
      }

      constexpr
      auto x() const
      {
        return m_ref.get().x();
      }

      constexpr
      auto y() const
      {
        return m_ref.get().y();
      }

      constexpr
      auto z() const
      {
        return m_ref.get().z();
      }

      template <class DataType>
      constexpr
      decltype(auto) setData(const DataType& data, size_t offset = 0)
      {
        return m_ref.get().setData(data, offset);
      }

      /**
       * @brief Returns a constant reference to the GridFunction data.
       */
      constexpr
      const auto& getData()
      {
        return m_ref.get().getData();
      }

      constexpr
      const auto& getFiniteElementSpace() const
      {
        return m_ref.get().getFiniteElementSpace();
      }

      constexpr
      size_t getSize() const
      {
        return m_ref.get().getSize();
      }

      Optional<size_t> getOrder(const Geometry::Polytope& geom) const
      {
        return m_ref.get().getOrder(geom);
      }

      GridFunctionBaseReference* copy() const noexcept final override
      {
        return new GridFunctionBaseReference(*this);
      }

    private:
      std::reference_wrapper<const Derived> m_ref;
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
    : public GridFunctionBaseReference<Derived>
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
      using Parent = GridFunctionBaseReference<Derived>;

      static_assert(
          std::is_same_v<RangeType, ScalarType> ||
          std::is_same_v<RangeType, Math::Vector<ScalarType>>);

      /**
       * @brief Constructs a grid function on the given finite element space.
       * @param[in] fes Finite element space
       *
       * Creates a grid function associated with the given FE space. The DOF
       * vector is sized according to the space dimension.
       */
      GridFunctionBase(const FES& fes)
        : Parent(std::cref(static_cast<const Derived&>(*this))),
          m_fes(std::cref(fes))
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Grid function to copy
       */
      GridFunctionBase(const GridFunctionBase& other)
        : Parent(std::cref(static_cast<const Derived&>(*this))),
          m_name(other.m_name),
          m_fes(other.m_fes)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Grid function to move from
       */
      GridFunctionBase(GridFunctionBase&& other)
        : Parent(std::cref(static_cast<const Derived&>(*this))),
          m_name(std::move(other.m_name)),
          m_fes(std::move(other.m_fes))
      {}

      virtual ~GridFunctionBase() = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other Grid function to move from
       * @return Reference to this grid function
       */
      GridFunctionBase& operator=(GridFunctionBase&& other)
      {
        m_name = std::move(other.m_name);
        m_fes = std::move(other.m_fes);
        return *this;
      }

      /**
       * @brief Copy assignment operator.
       * @param[in] other Grid function to copy
       * @return Reference to this grid function
       */
      GridFunctionBase& operator=(const GridFunctionBase& other)
      {
        if (this != &other)
        {
          m_name = other.m_name;
          m_fes = other.m_fes;
        }
        return *this;
      }

      /**
       * @brief Extracts the x-component of a vector-valued grid function.
       * @return Component operator for the first component
       *
       * For vector-valued functions @f$ \mathbf{u} = (u_x, u_y, u_z) @f$,
       * returns @f$ u_x @f$.
       */
      constexpr
      auto x() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 1);
        return Component(static_cast<const Derived&>(*this), 0);
      }

      /**
       * @brief Extracts the y-component of a vector-valued grid function.
       * @return Component operator for the second component
       *
       * For vector-valued functions @f$ \mathbf{u} = (u_x, u_y, u_z) @f$,
       * returns @f$ u_y @f$.
       */
      constexpr
      auto y() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 2);
        return Component(static_cast<const Derived&>(*this), 1);
      }

      /**
       * @brief Extracts the z-component of a vector-valued grid function.
       * @return Component operator for the third component
       *
       * For vector-valued functions @f$ \mathbf{u} = (u_x, u_y, u_z) @f$,
       * returns @f$ u_z @f$.
       */
      constexpr
      auto z() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 3);
        return Component(static_cast<const Derived&>(*this), 2);
      }

      /**
       * @brief Sets the DOF data vector.
       * @param[in] data DOF data to set
       * @param[in] offset Starting offset in the data vector
       * @return Reference to this grid function
       *
       * Copies the provided data into the internal DOF vector starting at
       * the specified offset.
       */
      constexpr
      Derived& setData(const DataType& data, size_t offset = 0)
      {
        return static_cast<Derived&>(*this).setData(data, offset);
      }

      /**
       * @brief Gets the DOF data vector.
       * @return Reference to the DOF data
       *
       * Provides access to the coefficient vector @f$ \mathbf{u} @f$ containing
       * the degrees of freedom.
       */
      constexpr
      auto& getData()
      {
        return static_cast<Derived&>(*this).getData();
      }

      /**
       * @brief Gets the DOF data vector (const version).
       * @return Const reference to the DOF data
       */
      constexpr
      const DataType& getData() const
      {
        return static_cast<const Derived&>(*this).getData();
      }

      /**
       * @brief Gets the associated finite element space.
       * @return Reference to the finite element space
       */
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      /**
       * @brief Gets the number of degrees of freedom.
       * @return Total number of DOFs
       */
      constexpr
      size_t getSize() const
      {
        return m_fes.get().getSize();
      }

      /**
       * @brief Gets the vector dimension of the grid function.
       * @return Number of components (1 for scalar, d for vector-valued)
       *
       * For scalar functions, returns 1. For vector-valued functions,
       * returns the space dimension.
       */
      constexpr
      size_t getDimension() const
      {
        return m_fes.get().getVectorDimension();
      }

      /**
       * @brief Loads grid function data from a file.
       * @param[in] filename Path to the file to load
       * @param[in] fmt File format (default: MFEM)
       * @return Reference to this grid function
       *
       * Reads DOF data from a file in the specified format. Supported formats
       * include MFEM, and MEDIT.
       */
      Derived& load(
          const boost::filesystem::path& filename,
          IO::FileFormat fmt)
      {
        switch (fmt)
        {
          case IO::FileFormat::MFEM:
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
            IO::GridFunctionLoader<IO::FileFormat::MFEM, FESType, DataType>(
              static_cast<Derived&>(*this)).load(input);
            break;
          }
          case IO::FileFormat::MEDIT:
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
            IO::GridFunctionLoader<IO::FileFormat::MEDIT, FES, DataType>(
              static_cast<Derived&>(*this)).load(input);
            break;
          }
          case IO::FileFormat::HDF5:
          {
            IO::GridFunctionLoader<IO::FileFormat::HDF5, FESType, DataType>(
              static_cast<Derived&>(*this)).load(filename);
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
          IO::FileFormat fmt) const
      {
        switch (fmt)
        {
          case IO::FileFormat::MFEM:
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
            IO::GridFunctionPrinter<IO::FileFormat::MFEM, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(output);
            break;
          }
          case IO::FileFormat::MEDIT:
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
            IO::GridFunctionPrinter<IO::FileFormat::MEDIT, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(output);
            break;
          }
          case IO::FileFormat::HDF5:
          {
            IO::GridFunctionPrinter<IO::FileFormat::HDF5, FESType, DataType>(
              static_cast<const Derived&>(*this)).print(filename);
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
      }

      /**
       * @brief Gets the interpolated value at the point.
       */
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local RangeType s_out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& fes = m_fes.get();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          static_cast<const Derived&>(*this).interpolate(s_out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          static_cast<const Derived&>(*this).interpolate(s_out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          if (restriction)
          {
            static_cast<const Derived&>(*this).interpolate(s_out, *restriction);
          }
          else
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Point is not contained in the finite element space mesh."
              << Alert::Raise;
          }
        }
        else
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Point is not contained in the finite element space mesh."
            << Alert::Raise;
        }
        return s_out;
      }

      constexpr
      void interpolate(RangeType& res, const Geometry::Point& p) const
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();

        const auto& fe = fes.getFiniteElement(d, i);
        const size_t count = fe.getCount();
        for (Index local = 0; local < count; ++local)
        {
          const auto mapping = fes.getPushforward({ d, i }, fe.getBasis(local));
          const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * mapping(p);
          if (local == 0)
            res = k;
          else
            res += k;
        }
      }

      template <class Function>
      Derived& project(const std::pair<size_t, Index>& p, const Function& fn)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& [d, i] = p;
        const auto& fe = fes.getFiniteElement(d, i);
        const auto mapping = fes.getPullback({ d, i }, fn);
        for (Index local = 0; local < fe.getCount(); local++)
        {
          const Index global = fes.getGlobalIndex({ d, i }, local);
          this->operator[](global) = fe.getLinearForm(local)(mapping);
        }
        return static_cast<Derived&>(*this);
      }

      Derived& project(const std::pair<size_t, Index>& p, const RangeType& v)
      {
        return static_cast<Derived&>(*this).project(
            p, [&](RangeType& out, const Geometry::Point&){ out = v; });
      }

      template <class T>
      Derived& operator=(const T& v)
      {
        return static_cast<Derived&>(*this).project(v);
      }

      template <class T>
      Derived& project(const T& fn)
      {
        return static_cast<Derived&>(*this).project(
            Geometry::Region::Cells, fn, [](const Geometry::Polytope&) { return true; });
      }

      template <class T>
      Derived& project(const Geometry::Region& region, const T& fn)
      {
        return static_cast<Derived&>(*this).project(region, fn,
            [](const Geometry::Polytope&) { return true; });
      }

      template <class T>
      Derived& project(
          const Geometry::Region& region, const T& fn, const Geometry::Attribute& attr)
      {
        return static_cast<Derived&>(*this).project(region, fn,
            [&](const Geometry::Polytope& polytope)
            { return polytope.getAttribute() == attr; });
      }

      template <class T>
      Derived& project(
          const Geometry::Region& region,
          const T& fn, const FlatSet<Geometry::Attribute>& attrs)
      {
        return static_cast<Derived&>(*this).project(region, fn,
            [&](const Geometry::Polytope& polytope)
            { return attrs.size() == 0 || attrs.count(polytope.getAttribute()); });
      }

      template <class Pred>
      Derived& project(const Geometry::Region& region, const RangeType& fn, const Pred& pred)
      {
        return static_cast<Derived&>(*this).project(
            region,
            [&](const Geometry::Point&) -> decltype(auto)
            { static thread_local RangeType s_out; s_out = fn; return s_out; }, pred);
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      template <class Function, class Pred>
      Derived& project(const Geometry::Region& region, const Function& fn, const Pred& pred)
      {
        return static_cast<Derived&>(*this).project(region, fn, pred);
      }

      /**
       * @brief Searches the minimum value in the grid function data.
       * @returns Minimum value in grid function.
       *
       * This function will compute the minimum value in the grid function
       * data array.
       *
       * @par Complexity
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
       * @par Complexity
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

      Optional<size_t> getOrder(const Geometry::Polytope& geom) const
      {
        return static_cast<const Derived&>(*this).getOrder(geom);
      }

      Optional<StringView> getName() const override
      {
        if (m_name)
          return StringView(m_name->c_str(), m_name->size());
        else
          return std::nullopt;
      }

      GridFunctionBase& setName(const std::string& name)
      {
        m_name = name;
        return *this;
      }

    private:
      Optional<std::string> m_name;
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
      using Parent::project;
      using Parent::min;
      using Parent::max;

      GridFunction(const FESType& fes)
        : Parent(fes)
      {
        auto& data = this->getData();
        data.resize(fes.getSize());
        data.setZero();
      }

      GridFunction(const GridFunction& other)
        : Parent(other),
          m_data(other.m_data)
      {}

      GridFunction(GridFunction&& other)
        : Parent(std::move(other)),
          m_data(std::move(other.m_data))
      {}

      GridFunction& operator=(GridFunction&& other)
      {
        Parent::operator=(std::move(other));
        m_data = std::move(other.m_data);
        return *this;
      }

      GridFunction& operator=(const GridFunction& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_data = other.m_data;
        }
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
        this->getData().array() += rhs;
        return *this;
      }

      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->getData().array() -= rhs;
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
        return *this;
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

      template <class Function, class Pred>
      GridFunction& project(
          const Geometry::Region& region, const Function& v, const Pred& pred)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Geometry::PolytopeIterator it;
        switch (region)
        {
          case Geometry::Region::Cells:
          {
            it = mesh.getCell();
            break;
          }
          case Geometry::Region::Faces:
          {
            it = mesh.getFace();
            break;
          }
          case Geometry::Region::Boundary:
          {
            it = mesh.getBoundary();
            break;
          }
          case Geometry::Region::Interface:
          {
            it = mesh.getInterface();
            break;
          }
        }

        while (it)
        {
          const auto& polytope = *it;
          if (pred(polytope))
            this->project({ polytope.getDimension(), polytope.getIndex() }, v);
          ++it;
        }

        return *this;
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

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const
      {
        return std::nullopt;
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
