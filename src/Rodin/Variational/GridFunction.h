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
#include <filesystem>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @internal
    * @brief Abstract class for GridFunction objects.
    */
   class GridFunctionBase
   {
      public:
         /**
          * @internal
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
         virtual void update() = 0;
         virtual GridFunctionBase& operator*=(double t) = 0;
         virtual GridFunctionBase& operator/=(double t) = 0;
         virtual double max() const = 0;
         virtual double min() const = 0;
   };

   /**
    * @brief Represents a grid function which does not yet have an associated
    * finite element space.
    *
    * To obtain the full functionality of the GridFunction class one must call
    * the @ref setFiniteElementSpace(FiniteElementSpace&) method.
    */
   class IncompleteGridFunction
   {
      public:
         /**
          * @brief Constructs an empty grid function with no associated finite
          * element space.
          */
         IncompleteGridFunction() = default;

         /**
          * @brief Associates a finite element space to the function.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          * @tparam FEC Finite element collection associated to the finite
          * element space.
          */
         template <class FEC>
         GridFunction<FEC> setFiniteElementSpace(FiniteElementSpace<FEC>& fes)
         {
            GridFunction<FEC> res(fes);
            int size = m_gf.Size();
            res.getHandle().SetDataAndSize(m_gf.StealData(), size);
            return res;
         }

         mfem::GridFunction& getHandle()
         {
            return m_gf;
         }

         const mfem::GridFunction& getHandle() const
         {
            return m_gf;
         }

      private:
         mfem::GridFunction m_gf;
   };

   /**
    * @brief Represents a grid function which belongs to some finite element space.
    *
    * @tparam FEC Finite element collection to which the function belongs.
    *
    * @note Note that the FEC template parameter is typically inferred when
    * initializing the grid function, hence it is not necessary to make it
    * explicit.
    */
   template <class FEC>
   class GridFunction : public GridFunctionBase
   {
      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(FiniteElementSpace<FEC>& fes)
            :  m_fes(fes),
               m_gf(&fes.getFES()),
               m_size(m_gf.Size()),
               m_data(m_gf.StealData())
         {
            assert(!m_gf.OwnsData());
            m_gf.SetDataAndSize(m_data.get(), m_size);
            m_gf = 0.0;
         }

         /**
          * @brief Copies the grid function.
          * @param[in] other Other grid function to copy.
          */
         GridFunction(const GridFunction& other)
            :  m_fes(other.m_fes),
               m_gf(other.m_gf),
               m_size(other.m_size)
         {
            assert(!m_gf.OwnsData());
            m_data = std::unique_ptr<double[]>(new double[other.m_size]);
            std::copy(
                  other.m_data.get(), other.m_data.get() + other.m_size,
                  m_data.get()
                  );
            m_gf.SetDataAndSize(m_data.get(), m_size);
         }

         /**
          * @brief Gets the finite element space to which the function
          * belongs to.
          * @returns Reference to finite element space to which the function
          * belongs to.
          */
         FiniteElementSpace<FEC>& getFiniteElementSpace() override
         {
            return m_fes;
         }

         /**
          * @brief Gets the finite element space to which the function
          * belongs to.
          * @returns Constant reference to finite element space to which the
          * function belongs to.
          */
         const FiniteElementSpace<FEC>& getFiniteElementSpace() const override
         {
            return m_fes;
         }

         /**
          * @brief Loads the grid function without assigning a finite element
          * space.
          * @param[in] filename Name of file to which the grid function will be
          * written to.
          */
         static IncompleteGridFunction load(const std::filesystem::path& filename)
         {
            std::ifstream in(filename);
            IncompleteGridFunction res;
            res.getHandle().Load(in);
            return res;
         }

         /**
          * @brief Saves the grid function in MFEMv1.0 format.
          * @param[in] filename Name of file to which the solution file will be
          * written.
          */
         void save(const std::filesystem::path& filename)
         {
            m_gf.Save(filename.string().c_str());
         }

         /**
          * @brief Sets the data of the grid function and assumes ownership.
          *
          * @param[in] data Data array
          * @param[in] size Size of the data array
          *
          * @returns Reference to self (for method chaining)
          */
         GridFunction& setData(std::unique_ptr<double[]> data, int size)
         {
            assert(data.get());
            m_data = std::move(data);
            m_gf.SetDataAndSize(m_data.get(), size);
            return *this;
         }

         /**
          * @brief Gets the raw data and size of the grid function.
          * @returns `std::pair{data, size}`
          */
         std::pair<const double*, int> getData() const
         {
            return {m_data.get(), m_size};
         }

         /**
          * @brief Projects a scalar coefficient on the given GridFunction.
          * @note The GridFunction must be scalar valued.
          * @param[in] s Scalar coefficient to project
          * @returns Reference to self
          */
         GridFunction<FEC>& operator=(const ScalarCoefficientBase& s)
         {
            assert(getFiniteElementSpace().getDimension() == 1);
            std::unique_ptr<ScalarCoefficientBase> sCopy(s.copy());
            sCopy->buildMFEMCoefficient();
            getHandle().ProjectCoefficient(sCopy->getMFEMCoefficient());
            return *this;
         }

         /**
          * @brief Projects a scalar coefficient on the given GridFunction.
          * @note The GridFunction must be vector valued.
          * @param[in] v Scalar coefficient to project
          * @returns Reference to self
          */
         GridFunction<FEC>& operator=(const VectorCoefficientBase& v)
         {
            assert(getFiniteElementSpace().getDimension() == v.getDimension());
            std::unique_ptr<VectorCoefficientBase> vCopy(v.copy());
            vCopy->buildMFEMVectorCoefficient();
            getHandle().ProjectCoefficient(vCopy->getMFEMVectorCoefficient());
            return *this;
         }

         GridFunction<FEC>& operator*=(double t) override
         {
            m_gf *= t;
            return *this;
         }

         GridFunction<FEC>& operator/=(double t) override
         {
            m_gf /= t;
            return *this;
         }

         void update() override
         {
            return m_gf.Update();
         }

         double max() const override
         {
            return m_gf.Max();
         }

         double min() const override
         {
            return m_gf.Min();
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
         std::reference_wrapper<FiniteElementSpace<FEC>> m_fes;
         mfem::GridFunction m_gf;
         int m_size;
         std::unique_ptr<double[]> m_data;
   };
}

#endif
