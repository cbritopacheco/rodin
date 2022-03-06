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
#include <filesystem>
#include <type_traits>

#include <mfem.hpp>

#include "Rodin/Core.h"
#include "Rodin/Alert.h"
#include "Rodin/Mesh/SubMesh.h"

#include "ForwardDecls.h"
#include "Restriction.h"
#include "ScalarCoefficient.h"
#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"
#include "GridFunctionView.h"
#include "GridFunctionIndex.h"

namespace Rodin::Variational
{
   /**
    * @internal
    * @brief Abstract class for GridFunction objects.
    */
   class GridFunctionBase
   {
      public:
         virtual void update() = 0;

         virtual double max() const = 0;
         virtual double min() const = 0;

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

         virtual GridFunctionBase& operator*=(double t) = 0;
         virtual GridFunctionBase& operator/=(double t) = 0;

         virtual GridFunctionBase& project(const ScalarCoefficientBase& s, int attr) = 0;
         virtual GridFunctionBase& project(const VectorCoefficientBase& s, int attr) = 0;
         virtual GridFunctionBase& project(const ScalarCoefficientBase& s, const std::set<int>& attrs) = 0;
         virtual GridFunctionBase& project(const VectorCoefficientBase& s, const std::set<int>& attrs) = 0;

         virtual GridFunctionBase& projectOnBoundary(const ScalarCoefficientBase& s, int attr) = 0;
         virtual GridFunctionBase& projectOnBoundary(const VectorCoefficientBase& s, int attr) = 0;
         virtual GridFunctionBase& projectOnBoundary(const ScalarCoefficientBase& s, const std::set<int>& attrs) = 0;
         virtual GridFunctionBase& projectOnBoundary(const VectorCoefficientBase& s, const std::set<int>& attrs) = 0;

         virtual std::pair<const double*, int> getData() const = 0;
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
          * @tparam FES Finite element space associated
          */
         template <class FES>
         GridFunction<FES> setFiniteElementSpace(FES& fes)
         {
            GridFunction<FES> res(fes);
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
    * @tparam FES Finite element collection to which the function belongs.
    *
    * @note Note that the FES template parameter is typically inferred when
    * initializing the grid function, hence it is not necessary to make it
    * explicit.
    */
   template <class FES>
   class GridFunction : public GridFunctionBase
   {
      static_assert(
            std::is_base_of_v<FiniteElementSpaceBase, FES>,
            "FES must be derived from FiniteElementSpaceBase");
      public:
         /**
          * @brief Constructs a grid function on a finite element space.
          * @param[in] fes Finite element space to which the function belongs
          * to.
          */
         GridFunction(FES& fes)
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

         GridFunction(GridFunction&& other)
            : m_fes(other.m_fes),
              m_gf(std::move(other.m_gf)),
              m_size(other.m_size)
         {}

         FES& getFiniteElementSpace() override
         {
            return m_fes;
         }

         /**
          * @brief Gets the finite element space to which the function
          * belongs to.
          * @returns Constant reference to finite element space to which the
          * function belongs to.
          */
         const FES& getFiniteElementSpace() const override
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
         std::pair<const double*, int> getData() const override
         {
            return {m_data.get(), m_size};
         }

         GridFunction& operator=(double v)
         {
            getHandle() = v;
            return *this;
         }

         template <class T>
         std::enable_if_t<
           std::is_base_of_v<GridFunctionIndexBase, T>, GridFunctionView>
         operator[](T&& idx)
         {
            GridFunctionView res(*this);
            res.setIndex(std::forward<T>(idx));
            return res;
         }

         template <class T>
         GridFunction& operator=(T&& v)
         {
            return project(std::forward<T>(v));
         }

         GridFunction& project(const ScalarCoefficientBase& s, int attr) override
         {
            return project(s, std::set<int>{attr});
         }

         GridFunction& project(const VectorCoefficientBase& v, int attr) override
         {
            return project(v, std::set<int>{attr});
         }

         GridFunction& projectOnBoundary(const ScalarCoefficientBase& s, int attr) override
         {
            return projectOnBoundary(s, std::set<int>{attr});
         }

         GridFunction& projectOnBoundary(const VectorCoefficientBase& v, int attr) override
         {
            return projectOnBoundary(v, std::set<int>{attr});
         }

         GridFunction& project(const ScalarCoefficientBase& s, const std::set<int>& attrs = {}) override
         {
            assert(getFiniteElementSpace().getVectorDimension() == 1);
            auto iv = s.build();

            if (attrs.size() == 0)
               getHandle().ProjectCoefficient(*iv);
            else
            {
               int maxAttr = getFiniteElementSpace()
                            .getMesh()
                            .getHandle().attributes.Max();
               mfem::Array<int> marker(maxAttr);
               marker = 0;
               for (const auto& attr : attrs)
               {
                  assert(attr - 1 < maxAttr);
                  marker[attr - 1] = 1;
               }
               getHandle().ProjectCoefficient(*iv, marker);
            }
            return *this;
         }

         GridFunction& project(const VectorCoefficientBase& s, const std::set<int>& attrs = {}) override
         {
            assert(getFiniteElementSpace().getVectorDimension() == s.getDimension());
            auto iv = s.build();

            if (attrs.size() == 0)
               getHandle().ProjectCoefficient(*iv);
            else
            {
               int maxAttr = getFiniteElementSpace()
                            .getMesh()
                            .getHandle().attributes.Max();
               mfem::Array<int> marker(maxAttr);
               marker = 0;
               for (const auto& attr : attrs)
               {
                  assert(attr - 1 < maxAttr);
                  marker[attr - 1] = 1;
               }
               getHandle().ProjectCoefficient(*iv, marker);
            }
            return *this;
         }

         GridFunction& projectOnBoundary(const ScalarCoefficientBase& s, const std::set<int>& attrs = {}) override
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

         GridFunction& projectOnBoundary(const VectorCoefficientBase& v, const std::set<int>& attrs = {}) override
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

         /**
          * @brief Projects the restriction of a scalar coefficient on the given GridFunction.
          * @note The GridFunction must be scalar valued.
          * @param[in] s Scalar coefficient to project
          * @returns Reference to self
          */
         GridFunction& project(const Restriction<ScalarCoefficientBase>& s)
         {
            assert(getFiniteElementSpace().getVectorDimension() == 1);
            auto iv = s.getScalarCoefficient().build();
            getHandle() = NAN;
            mfem::Array<int> vdofs;
            mfem::Vector vals;
            const auto& fes = getFiniteElementSpace().getFES();
            const auto& attrs = s.getAttributes();
            for (int i = 0; i < fes.GetNE(); i++)
            {
               if (attrs.count(fes.GetAttribute(i)) > 0)
               {
                  fes.GetElementVDofs(i, vdofs);
                  vals.SetSize(vdofs.Size());
                  fes.GetFE(i)->Project(
                        *iv, *fes.GetElementTransformation(i), vals);
                  getHandle().SetSubVector(vdofs, vals);
               }
            }
            return *this;
         }

         /**
          * @brief Transfers the grid function from one finite element space to
          * another.
          */
         template <class OtherFES>
         void transfer(GridFunction<OtherFES>& other)
         {
            assert(getFiniteElementSpace().getVectorDimension() ==
                  other.getFiniteElementSpace().getVectorDimension());
            if (getFiniteElementSpace().getMesh().isSubMesh())
            {
               // If we are here the this means that we are in a submesh of the
               // underlying target finite element space. Hence we should seek
               // out to copy the grid function at the corresponding nodes
               // given by the vertex map given in the Submesh object.
               auto& submesh = static_cast<SubMesh&>(getFiniteElementSpace().getMesh());
               if (&submesh.getParent() == &other.getFiniteElementSpace().getMesh())
               {
                  int vdim = getFiniteElementSpace().getVectorDimension();
                  const auto& s2pv = submesh.getVertexMap();
                  if (vdim == 1)
                  {
                     int size = getHandle().Size();
                     for (int i = 0; i < size; i++)
                        other.getHandle()[i] = getHandle()[s2pv.at(i)];
                  }
                  else
                  {
                     int nv = getFiniteElementSpace().getFES().GetNV();
                     int pnv = other.getFiniteElementSpace().getFES().GetNV();

                     assert(getFiniteElementSpace().getFES().GetOrdering() ==
                              getFiniteElementSpace().getFES().GetOrdering());
                     switch(getFiniteElementSpace().getFES().GetOrdering())
                     {
                        case mfem::Ordering::byNODES:
                        {
                           for (int i = 0; i < vdim; i++)
                              for (int j = 0; j < nv; j++)
                                 other.getHandle()[s2pv.at(j) + i * pnv] = getHandle()[j + i * nv];
                           return;
                        }
                        case mfem::Ordering::byVDIM:
                        {
                           for (int i = 0; i < nv; i++)
                              for (int j = 0; j < vdim; j++)
                                 other.getHandle()[s2pv.at(i) * vdim + j] = getHandle()[i * vdim + j];
                           return;
                        }
                     }
                  }
               }
            }
            else
            {
               // If the meshes are equal or where obtained from refinements
               // one could use the mfem functionality to make a GridTransfer.
               // Alternatively, if the mesh is equal but the finite element
               // spaces are not, mfem also contains the TransferOperator class
               // which can come in useful.
               Alert::Exception("Unimplemented. Sorry.").raise();
            }
         }

         GridFunction& operator*=(double t) override
         {
            m_gf *= t;
            return *this;
         }

         GridFunction& operator/=(double t) override
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
         FES& m_fes;
         mfem::GridFunction m_gf;
         int m_size;
         std::unique_ptr<double[]> m_data;
   };
}

#endif
