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

#include "Rodin/Cast.h"
#include "Rodin/Math.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/IO/ForwardDecls.h"

#include "ForwardDecls.h"
#include "H1.h"
#include "L2.h"
#include "Restriction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "FiniteElementSpace.h"

#include "Function.h"
#include "Exceptions.h"

namespace Rodin::Variational
{

   /**
    * @defgroup GridFunctionSpecializations GridFunction Template Specializations
    * @brief Template specializations of the GridFunction class.
    * @see GridFunction
    */

   class GridFunctionBase : public VectorFunctionBase
   {
      public:
         class GridFunctionValue : public FunctionValue
         {
            public:
               constexpr
               GridFunctionValue(double v)
                  : FunctionValue(v)
               {}

               constexpr
               GridFunctionValue(const mfem::Vector& v)
                  : FunctionValue(v)
               {}

               constexpr
               GridFunctionValue(mfem::Vector&& v)
                  : FunctionValue(std::move(v))
               {}

               GridFunctionValue(const GridFunctionValue&) = default;

               GridFunctionValue(GridFunctionValue&&) = default;
         };

         GridFunctionBase() = default;

         GridFunctionBase(const GridFunctionBase& other)
            : VectorFunctionBase(other)
         {}

         GridFunctionBase(GridFunctionBase&& other)
            : VectorFunctionBase(std::move(other))
         {}

         GridFunctionBase& operator=(GridFunctionBase&& other)
         {
            FunctionBase::operator=(std::move(other));
            return *this;
         }

         GridFunctionBase& operator=(const GridFunctionBase&) = delete;

         GridFunctionValue operator()(const Geometry::Point& v) const
         {
            switch (getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::Vector value;
                  getValue(value, *v.getElementTransformation(), *v.getIntegrationPoint());
                  return GridFunctionValue(value(0));
               }
               case RangeType::Vector:
               {
                  mfem::Vector value;
                  getValue(value, *v.getElementTransformation(), *v.getIntegrationPoint());
                  return GridFunctionValue(std::move(value));
               }
               case RangeType::Matrix:
               {
                  assert(false);
                  return 0.0;
               }
            }
         }

         /**
          * @brief Searches for the maximum value in the grid function data.
          * @returns Maximum value in grid function.
          *
          * This function will compute the maximum value in the grid function
          * data array.
          */
         double max() const;

         /**
          * @brief Searches the minimum value in the grid function data.
          * @returns Minimum value in grid function.
          *
          * This function will compute the minimum value in the grid function
          * data array.
          */
         double min() const;

         Component<FunctionBase> x() const;

         Component<FunctionBase> y() const;

         Component<FunctionBase> z() const;

         /**
          * @brief Updates the state after a refinement in the mesh.
          *
          * This method will update the grid function after a call to the
          * @ref MeshBase::refine() "refine()" method.
          */
         GridFunctionBase& update();

         /**
          * @brief Bulk assigns the value to the whole data array.
          */
         GridFunctionBase& operator=(double v);

         /**
          * @brief Gets the raw data and size of the grid function.
          * @returns `std::pair{data, size}`
          */
         std::pair<const double*, int> getData() const;

         /**
          * @brief Sets the data of the grid function and assumes ownership.
          *
          * @param[in] data Data array
          * @param[in] size Size of the data array
          *
          * @returns Reference to self (for method chaining)
          */
         GridFunctionBase& setData(std::unique_ptr<double[]> data, int size);

         virtual void save(
               const boost::filesystem::path& filename, IO::FileFormat fmt,
               int precision) const = 0;

         virtual GridFunctionBase& load(
               const boost::filesystem::path& filename, IO::FileFormat fmt) = 0;

         /**
          * @brief Addition of a scalar value.
          */
         GridFunctionBase& operator+=(double t);

         GridFunctionBase& operator+=(const GridFunctionBase& rhs);

         /**
          * @brief Substraction of a scalar value.
          */
         GridFunctionBase& operator-=(double t);

         GridFunctionBase& operator-=(const GridFunctionBase& rhs);

         /**
          * @brief Multiplication by a scalar value.
          */
         GridFunctionBase& operator*=(double t);

         GridFunctionBase& operator*=(const GridFunctionBase& rhs);

         /**
          * @brief Division by a scalar value.
          */
         GridFunctionBase& operator/=(double t);

         GridFunctionBase& operator/=(const GridFunctionBase& rhs);

         /**
          * @brief Projection of a function.
          */
         GridFunctionBase& operator=(const FunctionBase& v)
         {
            return project(v);
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
         GridFunctionBase& project(const FunctionBase& v, int attr)
         {
            return project(v, std::set<int>{attr});
         }

         /**
          * @brief Projects a FunctionBase instance on the grid function.
          *
          * This function will project a FunctionBase instance on the
          * domain elements with the given attributes. If the attribute set is
          * empty, this function will project over all elements in the mesh.
          */
         GridFunctionBase& project(
               const FunctionBase& s, const std::set<int>& attrs = {});

         /**
          * @brief Projects the restriction of a scalar coefficient on the given GridFunction.
          * @note The GridFunction must be scalar valued.
          * @param[in] s Scalar coefficient to project
          * @returns Reference to self
          */
         GridFunctionBase& project(const Restriction<FunctionBase>& s);

         /**
          * @brief Transfers the grid function from one finite element space to
          * another.
          * @param[in, out] dst Destination GridFunction for the transfer
          */
         void transfer(GridFunctionBase& dst);

         std::set<Geometry::Point> where(
               const BooleanFunctionBase& p,
               const std::set<int>& attrs = {},
               std::function<int(mfem::ElementTransformation&)> order =
                  [](mfem::ElementTransformation&) { return 1; }) const;


         RangeType getRangeType() const override;

         RangeShape getRangeShape() const override;

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         int getDimension() const override;

         /**
          * @brief Gets the underlying handle to the mfem::GridFunction object.
          * @returns Reference to the underlying object.
          */
         virtual mfem::GridFunction& getHandle() = 0;

         /**
          * @internal
          * @brief Gets the underlying handle to the mfem::GridFunction object.
          * @returns Constant reference to the underlying object.
          */
         virtual const mfem::GridFunction& getHandle() const = 0;

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

      private:
         VectorFunctionBase* copy() const noexcept override;
   };

   /**
    * @ingroup GridFunctionSpecializations
    * @brief Represents a GridFunction which belongs to an L2 finite element
    * space.
    */
   template <class Trait>
   class GridFunction<L2<Trait>> : public GridFunctionBase
   {
      public:
         /**
          * @brief Constructs a grid function on an L2 finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(L2<Trait>& fes)
            :  GridFunctionBase(),
               m_fes(fes),
               m_gf(&fes.getHandle())
         {
            m_gf = 0.0;
         }

         /**
          * @brief Copies the grid function.
          * @param[in] other Other grid function to copy.
          */
         GridFunction(const GridFunction& other)
            :  GridFunctionBase(other),
               m_fes(other.m_fes),
               m_gf(other.m_gf)
         {}

         /**
          * @brief Move constructs the grid function.
          * @param[in] other Other grid function to move.
          */
         GridFunction(GridFunction&& other)
            :  GridFunctionBase(std::move(other)),
               m_fes(std::move(other.m_fes)),
               m_gf(std::move(other.m_gf))
         {}

         /**
          * @brief Move assignment operator.
          */
         GridFunction& operator=(GridFunction&& other) = default;

         GridFunction& operator=(const GridFunction&)  = delete;

         template <class T>
         GridFunction& operator=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator+=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator+=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator-=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator-=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator*=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator*=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator/=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator/=(std::forward<T>(v)));
         }

         template <class ... Args>
         GridFunction& project(Args&&... args)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::project(std::forward<Args>(args)...));
         }

         void save(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM,
               int precision = 16) const override;

         GridFunction& load(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM) override;

         L2<Trait>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const L2<Trait>& getFiniteElementSpace() const override

         {
            return m_fes.get();
         }

         mfem::GridFunction& getHandle() override
         {
            return m_gf;
         }

         const mfem::GridFunction& getHandle() const override
         {
            return m_gf;
         }
      private:
         std::reference_wrapper<L2<Trait>> m_fes;
         mfem::GridFunction m_gf;
   };

   /**
    * @ingroup GridFunctionSpecializations
    * @brief Represents a GridFunction which belongs to an H1 finite element
    * space.
    */
   template <class Trait>
   class GridFunction<H1<Trait>> : public GridFunctionBase
   {
      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(H1<Trait>& fes)
            :  GridFunctionBase(),
               m_fes(fes),
               m_gf(&fes.getHandle())
         {
            m_gf = 0.0;
         }

         /**
          * @brief Copies the grid function.
          * @param[in] other Other grid function to copy.
          */
         GridFunction(const GridFunction& other)
            :  GridFunctionBase(other),
               m_fes(other.m_fes),
               m_gf(other.m_gf)
         {}

         /**
          * @brief Move constructs the grid function.
          * @param[in] other Other grid function to move.
          */
         GridFunction(GridFunction&& other)
            :  GridFunctionBase(std::move(other)),
               m_fes(std::move(other.m_fes)),
               m_gf(std::move(other.m_gf))
         {}

         /**
          * @brief Move assignment operator.
          */
         GridFunction& operator=(GridFunction&& other) = default;

         GridFunction& operator=(const GridFunction&)  = delete;

         template <class T>
         GridFunction& operator=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator+=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator+=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator-=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator-=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator*=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator*=(std::forward<T>(v)));
         }

         template <class T>
         GridFunction& operator/=(T&& v)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::operator/=(std::forward<T>(v)));
         }

         template <class ... Args>
         GridFunction& project(Args&&... args)
         {
            return static_cast<GridFunction&>(
                  GridFunctionBase::project(std::forward<Args>(args)...));
         }

         GridFunction& projectOnBoundary(const FunctionBase& v, int attr)
         {
            return projectOnBoundary(v, std::set<int>{attr});
         }

         GridFunction& projectOnBoundary(
               const FunctionBase& s, const std::set<int>& attrs = {})
         {
            auto va = s.build();
            switch (s.getRangeType())
            {
               case RangeType::Scalar:
               {
                  int maxBdrAttr = getFiniteElementSpace()
                                  .getMesh()
                                  .getHandle().bdr_attributes.Max();
                  mfem::Array<int> marker(maxBdrAttr);
                  if (attrs.size() == 0)
                  {
                     marker = 1;
                     getHandle().ProjectBdrCoefficient(va.get<RangeType::Scalar>(), marker);
                  }
                  else
                  {
                     marker = 0;
                     for (const auto& attr : attrs)
                     {
                        assert(attr - 1 < maxBdrAttr);
                        marker[attr - 1] = 1;
                     }
                     getHandle().ProjectBdrCoefficient(va.get<RangeType::Scalar>(), marker);
                  }
                  break;
               }
               case RangeType::Vector:
               {
                  int maxBdrAttr = getFiniteElementSpace()
                                  .getMesh()
                                  .getHandle().bdr_attributes.Max();
                  mfem::Array<int> marker(maxBdrAttr);
                  if (attrs.size() == 0)
                  {
                     marker = 1;
                     getHandle().ProjectBdrCoefficient(va.get<RangeType::Vector>(), marker);
                  }
                  else
                  {
                     marker = 0;
                     for (const auto& attr : attrs)
                     {
                        assert(attr - 1 < maxBdrAttr);
                        marker[attr - 1] = 1;
                     }
                     getHandle().ProjectBdrCoefficient(va.get<RangeType::Vector>(), marker);
                  }
                  break;
               }
               case RangeType::Matrix:
               {
                  UnexpectedRangeTypeException(
                        {RangeType::Scalar, RangeType::Vector}, RangeType::Matrix).raise();
                  break;
               }
            }
            return *this;
         }

         void save(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM,
               int precision = 16) const override;

         GridFunction& load(
               const boost::filesystem::path& filename,
               IO::FileFormat fmt = IO::FileFormat::MFEM) override;

         H1<Trait>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const H1<Trait>& getFiniteElementSpace() const override

         {
            return m_fes.get();
         }

         mfem::GridFunction& getHandle() override
         {
            return m_gf;
         }

         const mfem::GridFunction& getHandle() const override
         {
            return m_gf;
         }
      private:
         std::reference_wrapper<H1<Trait>> m_fes;
         mfem::GridFunction m_gf;
   };

   template <class FES>
   GridFunction(FES& fes) -> GridFunction<FES>;
}

namespace Rodin::Variational::Internal
{
   class GridFunctionEvaluator : public VectorFunctionBase
   {
      public:
         GridFunctionEvaluator(const GridFunctionBase& gf);

         GridFunctionEvaluator(const GridFunctionEvaluator& other);

         GridFunctionEvaluator(GridFunctionEvaluator&& other);

         RangeShape getRangeShape() const override;

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         int getDimension() const override
         {
            return m_gf.getDimension();
         }

         GridFunctionEvaluator* copy() const noexcept override
         {
            return new GridFunctionEvaluator(*this);
         }

      private:
         const GridFunctionBase& m_gf;
   };
}


#include "GridFunction.hpp"

#endif
