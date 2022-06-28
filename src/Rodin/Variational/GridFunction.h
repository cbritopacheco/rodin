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

         void save(const boost::filesystem::path& filename);

         GridFunctionBase& load(const boost::filesystem::path& filename);

         /**
          * @brief Addition of a scalar value.
          */
         GridFunctionBase& operator+=(double t);

         /**
          * @brief Substraction of a scalar value.
          */
         GridFunctionBase& operator-=(double t);

         /**
          * @brief Multiplication by a scalar value.
          */
         GridFunctionBase& operator*=(double t);

         /**
          * @brief Division by a scalar value.
          */
         GridFunctionBase& operator/=(double t);

         /**
          * @brief Projection of a scalar function.
          */
         GridFunctionBase& operator=(const ScalarFunctionBase& v)
         {
            return project(v);
         }

         /**
          * @brief Projection of a vector function.
          */
         GridFunctionBase& operator=(const VectorFunctionBase& v)
         {
            return project(v);
         }

         /**
          * @brief Projects a ScalarFunctionBase instance
          *
          * This function will project a ScalarFunctionBase instance on the
          * domain elements with the given attribute.
          *
          * It is a convenience function to call
          * project(const ScalarFunctionBase&, const std::set<int>&) with one
          * attribute.
          */
         GridFunctionBase& project(const ScalarFunctionBase& v, int attr)
         {
            return project(v, std::set<int>{attr});
         }

         /**
          * @brief Projects a VectorFunctionBase instance
          *
          * This function will project a VectorFunctionBase instance on the
          * domain elements with the given attribute.
          *
          * It is a convenience function to call
          * project(const VectorFunctionBase&, const std::set<int>&) with one
          * attribute.
          */
         GridFunctionBase& project(const VectorFunctionBase& v, int attr)
         {
            return project(v, std::set<int>{attr});
         }

         /**
          * @brief Projects a ScalarFunctionBase instance on the grid function.
          *
          * This function will project a ScalarFunctionBase instance on the
          * domain elements with the given attributes. If the attribute set is
          * empty, this function will project over all elements in the mesh.
          */
         GridFunctionBase& project(const ScalarFunctionBase& s, const std::set<int>& attrs = {});

         /**
          * @brief Projects a VectorFunctionBase instance on the grid function.
          *
          * This function will project a VectorFunctionBase instance on the
          * domain elements with the given attributes. If the attribute set is
          * empty, this function will project over all elements in the mesh.
          */
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
         GridFunction(FiniteElementSpace<H1, Trait>& fes)
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

         FiniteElementSpace<H1, Trait>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const FiniteElementSpace<H1, Trait>& getFiniteElementSpace() const override

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
         std::reference_wrapper<FiniteElementSpace<H1, Trait>> m_fes;
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
         GridFunction(FiniteElementSpace<L2, Trait>& fes)
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

         FiniteElementSpace<L2, Trait>& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const FiniteElementSpace<L2, Trait>& getFiniteElementSpace() const override
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
         std::reference_wrapper<FiniteElementSpace<L2, Trait>> m_fes;
         mfem::GridFunction m_gf;
   };

   template <class FEC, class Trait>
   GridFunction(FiniteElementSpace<FEC, Trait>&) -> GridFunction<FEC, Trait>;
}

#endif
