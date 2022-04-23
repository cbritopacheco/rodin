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
#include "Rodin/Core.h"
#include "Rodin/Alert.h"
#include "Rodin/Mesh/SubMesh.h"

#include "ForwardDecls.h"
#include "H1.h"
#include "L2.h"
#include "Restriction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   /**
    * @internal
    * @brief Abstract class for GridFunction objects.
    */
   class GridFunctionBase
   {
      public:
         GridFunctionBase() = default;

         GridFunctionBase(const GridFunctionBase&) = default;

         GridFunctionBase(GridFunctionBase&&) = default;

         GridFunctionBase& operator=(GridFunctionBase&&) = default;

         GridFunctionBase& operator=(const GridFunctionBase&) = delete;

         double max() const;

         double min() const;

         GridFunctionBase& update();

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

         void save(const boost::filesystem::path& filename);

         GridFunctionBase& load(const boost::filesystem::path& filename);

         GridFunctionBase& operator+=(double t);
         GridFunctionBase& operator-=(double t);
         GridFunctionBase& operator*=(double t);
         GridFunctionBase& operator/=(double t);

         GridFunctionBase& operator+=(const ScalarFunctionBase& v);
         GridFunctionBase& operator-=(const ScalarFunctionBase& v);
         GridFunctionBase& operator*=(const ScalarFunctionBase& v);
         GridFunctionBase& operator/=(const ScalarFunctionBase& v);

         GridFunctionBase& operator+=(const VectorFunctionBase& v);
         GridFunctionBase& operator-=(const VectorFunctionBase& v);

         GridFunctionBase& operator=(const ScalarFunctionBase& v)
         {
            return project(v);
         }

         GridFunctionBase& operator=(const VectorFunctionBase& v)
         {
            return project(v);
         }

         GridFunctionBase& project(const ScalarFunctionBase& v, int attr)
         {
            return project(v, std::set<int>{attr});
         }

         GridFunctionBase& project(const VectorFunctionBase& v, int attr)
         {
            return project(v, std::set<int>{attr});
         }

         GridFunctionBase& project(const ScalarFunctionBase& s, const std::set<int>& attrs = {});

         GridFunctionBase& project(const VectorFunctionBase& s, const std::set<int>& attrs = {});

         /**
          * @brief Transfers the grid function from one finite element space to
          * another.
          * @param[in, out] dst Destination GridFunction for the transfer
          */
         void transfer(GridFunctionBase& dst);

         /**
          * @brief Projects the restriction of a scalar coefficient on the given GridFunction.
          * @note The GridFunction must be scalar valued.
          * @param[in] s Scalar coefficient to project
          * @returns Reference to self
          */
         GridFunctionBase& project(const Restriction<ScalarFunctionBase>& s);

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
   };

   /**
    * @brief Represents a grid function which belongs to some finite element space.
    *
    * @tparam FES Finite element collection to which the function belongs.
    *
    * @note Note that the FES template parameter is typically inferred when
    * initializing the grid function, hence it is not necessary to make it
    * explicit.
    */
   template <class Trait>
   class GridFunction<H1, Trait> : public GridFunctionBase
   {
      static_assert(std::is_same_v<Trait, Traits::Serial>,
            "Parallel GridFunction is not supported yet. Please use Trait = Traits::Serial.");
      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(FiniteElementSpace<H1>& fes)
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

         GridFunction& projectOnBoundary(const ScalarFunctionBase& s, int attr)
         {
            return projectOnBoundary(s, std::set<int>{attr});
         }

         GridFunction& projectOnBoundary(const VectorFunctionBase& v, int attr)
         {
            return projectOnBoundary(v, std::set<int>{attr});
         }

         GridFunction& projectOnBoundary(
               const ScalarFunctionBase& s, const std::set<int>& attrs = {})
         {
            assert(getFiniteElementSpace().getVectorDimension() == 1);
            auto iv = s.build();
            int maxBdrAttr = getFiniteElementSpace()
                            .getMesh()
                            .getHandle().bdr_attributes.Max();
            mfem::Array<int> marker(maxBdrAttr);
            if (attrs.size() == 0)
            {
               marker = 1;
               getHandle().ProjectBdrCoefficient(*iv, marker);
            }
            else
            {
               marker = 0;
               for (const auto& attr : attrs)
               {
                  assert(attr - 1 < maxBdrAttr);
                  marker[attr - 1] = 1;
               }
               getHandle().ProjectBdrCoefficient(*iv, marker);
            }
            return *this;
         }

         GridFunction& projectOnBoundary(
               const VectorFunctionBase& v, const std::set<int>& attrs = {})
         {
            assert(getFiniteElementSpace().getVectorDimension() == v.getDimension());
            auto iv = v.build();
            int maxBdrAttr = getFiniteElementSpace()
                            .getMesh()
                            .getHandle().bdr_attributes.Max();
            mfem::Array<int> marker(maxBdrAttr);
            if (attrs.size() == 0)
            {
               marker = 1;
               getHandle().ProjectBdrCoefficient(*iv, marker);
            }
            else
            {
               marker = 0;
               for (const auto& attr : attrs)
               {
                  assert(attr - 1 < maxBdrAttr);
                  marker[attr - 1] = 1;
               }
               getHandle().ProjectBdrCoefficient(*iv, marker);
            }
            return *this;
         }

         FiniteElementSpace<H1>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         /**
          * @brief Gets the finite element space to which the function
          * belongs to.
          * @returns Constant reference to finite element space to which the
          * function belongs to.
          */
         const FiniteElementSpace<H1>& getFiniteElementSpace() const override
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
         std::reference_wrapper<FiniteElementSpace<H1>> m_fes;
         mfem::GridFunction m_gf;
   };

   template <class Trait>
   class GridFunction<L2, Trait> : public GridFunctionBase
   {
      static_assert(std::is_same_v<Trait, Traits::Serial>,
            "Parallel GridFunction is not supported yet. Please use Trait = Traits::Serial.");

      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(FiniteElementSpace<L2>& fes)
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

         FiniteElementSpace<L2>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         /**
          * @brief Gets the finite element space to which the function
          * belongs to.
          * @returns Constant reference to finite element space to which the
          * function belongs to.
          */
         const FiniteElementSpace<L2>& getFiniteElementSpace() const override
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
         std::reference_wrapper<FiniteElementSpace<L2>> m_fes;
         mfem::GridFunction m_gf;
   };

   template <class FEC, class Trait>
   GridFunction(FiniteElementSpace<FEC, Trait>&) -> GridFunction<FEC, Trait>;
}

#endif
