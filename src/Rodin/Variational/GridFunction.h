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

#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeIterator.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/EnSight6.h"
#include "Rodin/IO/MEDIT.h"
#include "Rodin/IO/MFEM.h"

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
   * - **I/O Support**: Export to various visualization formats (EnSight, MEDIT, MFEM)
   * - **Space Association**: Strong association with underlying finite element space
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
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return m_ref.get().getValue(p);
      }

      constexpr
      decltype(auto) x() const
      {
        return m_ref.get().x();
      }

      constexpr
      decltype(auto) y() const
      {
        return m_ref.get().y();
      }

      constexpr
      decltype(auto) z() const
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
      auto& getData()
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

      GridFunctionBase& operator=(const GridFunctionBase& other)
      {
        if (this != &other)
        {
          m_fes = other.m_fes;
        }
        return *this;
      }

      constexpr
      auto x() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 1);
        return Component(static_cast<const Derived&>(*this), 0);
      }

      constexpr
      auto y() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 2);
        return Component(static_cast<const Derived&>(*this), 1);
      }

      constexpr
      auto z() const
      {
        static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
        assert(m_fes.get().getVectorDimension() >= 3);
        return Component(static_cast<const Derived&>(*this), 2);
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
        return m_fes.get().getSize();
      }

      constexpr
      size_t getDimension() const
      {
        return m_fes.get().getVectorDimension();
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
            assert(false);
          }
        }
        else
        {
          assert(false);
        }
        return s_out;
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

      void project(const std::pair<size_t, Index>& p, const RangeType& v)
      {
        static_cast<Derived&>(*this).project(p,
            [&](RangeType& out, const Geometry::Point& pt){ out = v; });
      }

      /**
       * @note CRTP function to be overriden in Derived class.
       */
      template <class Function>
      void project(const std::pair<size_t, Index>& p, const Function& fn)
      {
        static_cast<Derived&>(*this).project(fn, p);
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
        : Parent(other)
      {}

      GridFunction(GridFunction&& other)
        : Parent(std::move(other))
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
        for (Index local = 0; local < count; ++local)
        {
          const auto mapping = fes.getInverseMapping({ d, i }, fe.getBasis(local));
          const auto k = this->operator[](fes.getGlobalIndex({ d, i }, local)) * mapping(p);
          if (local == 0)
            res = k; // Initializes the result (resizes)
          else
            res += k; // Accumulates the result (does not resize)
        }
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

      template <class Function>
      void project(const std::pair<size_t, Index>& p, const Function& fn)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& [d, i] = p;
        const auto& fe = fes.getFiniteElement(d, i);
        const auto mapping = fes.getMapping({ d, i }, fn);
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
