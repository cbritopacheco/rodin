/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects of type BilinearForm.
    */
   class BilinearFormBase
   {
      public:
         /**
          * @brief Updates the state after a refinement in the mesh.
          *
          * This method will update the bilinear form after a call to the
          * @ref MeshBase::refine() "refine()" method.
          */
         void update()
         {
            getHandle().Update();
         }

         /**
          * @brief Assembles the bilinear form.
          *
          * This method will assemble the underlying sparse matrix associated
          * the bilinear form.
          *
          * @see getMatrix()
          */
         void assemble()
         {
            getHandle().Assemble();
         }

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Reference to the associated sparse matrix.
          */
         mfem::SparseMatrix& getMatrix()
         {
            return getHandle().SpMat();
         }

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Constant reference to the associated sparse matrix.
          */
         const mfem::SparseMatrix& getMatrix() const
         {
            return getHandle().SpMat();
         }

         /**
          * @internal
          * @returns Underlying internal reference to bilinear form.
          */
         virtual mfem::BilinearForm& getHandle() = 0;

         /**
          * @internal
          * @returns Underlying internal constant reference to bilinear form.
          */
         virtual const mfem::BilinearForm& getHandle() const = 0;

         /**
          * @brief Builds the bilinear form the given bilinear integrator
          * @param[in] bfi Bilinear integrator which will be used to
          * build the bilinear form.
          * @returns Reference to this (for method chaining)
          */
         virtual BilinearFormBase& from(
               const BilinearFormIntegratorBase& bfi) = 0;

         /**
          * @todo
          */
         virtual BilinearFormBase& from(
               const FormLanguage::BilinearFormIntegratorSum& bfi) = 0;

         /**
          * @brief Adds a bilinear integrator to the bilinear form.
          * @returns Reference to this (for method chaining)
          */
         virtual BilinearFormBase& add(
               const BilinearFormIntegratorBase& bfi) = 0;

         /**
          * @todo
          */
         virtual BilinearFormBase& add(
               const FormLanguage::BilinearFormIntegratorSum& lsum) = 0;

         /**
          * @brief Gets the reference to the associated TrialFunction object.
          * @returns Reference to this (for method chaining)
          */
         virtual const ShapeFunctionBase<TrialSpace>& getTrialFunction() const = 0;

         /**
          * @brief Gets the reference to the associated TestFunction object.
          * @returns Reference to this (for method chaining)
          */
         virtual const ShapeFunctionBase<TestSpace>& getTestFunction() const = 0;
   };

   /**
    * @brief Represents a serial bilinear form supported on two finite element
    * spaces originating from two instances of FiniteElementCollection.
    * @tparam TrialFEC Trial FiniteElementCollection
    * @tparam TestFEC Test FiniteElementCollection
    *
    * An object of type BilinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    a : U \times V &\rightarrow \mathbb{R}\\
    *        (u, v) &\mapsto a(u, v)
    * \end{aligned}
    * @f]
    * where @f$ U @f$ and @f$ V @f$ are finite element spaces.
    *
    * A bilinear form may be built using the form language. For example,
    * @code{.cpp}
    * // Define spaces
    * FiniteElementSpace<H1> Vh;
    * TrialFunction u(Vh);
    * TestFunction  v(Vh);
    *
    * // Define mass bilinear form
    * BilinearForm bf(u, v);
    * bf = Integral(u, v);
    * @endcode
    */
   template <class TrialFEC, class TestFEC>
   class BilinearForm<TrialFEC, TestFEC, Traits::Serial>
      : public BilinearFormBase
   {
      static_assert(
            std::is_same_v<TrialFEC, TestFEC>,
            "Different trial and test spaces are currently not supported.");

      public:
         using BFIList = std::vector<std::unique_ptr<BilinearFormIntegratorBase>>;

         /**
          * @brief Constructs a BilinearForm from a TrialFunction and
          * TestFunction.
          *
          * @param[in] u TrialFunction
          * @param[out] v TestFunction
          */
         BilinearForm(
               TrialFunction<TrialFEC, Traits::Serial>& u,
               TestFunction<TestFEC, Traits::Serial>& v);

         /**
          * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
          * v @f$.
          *
          * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
          * compute the action of the bilinear mapping @f$ a(u, v) @f$.
          *
          * @returns The action @f$ a(u, v) @f$ which the bilinear form takes
          * at @f$ ( u, v ) @f$.
          */
         double operator()(
               const GridFunction<TrialFEC, Traits::Serial>& u,
               const GridFunction<TestFEC, Traits::Serial>& v) const;

         /**
          * @brief Builds the bilinear form from a derived instance of
          * BilinearFormIntegratorBase.
          * @param[in] bfi Bilinear form integrator
          */
         BilinearForm& operator=(const BilinearFormIntegratorBase& bfi);

         /**
          * @todo
          */
         BilinearForm& operator=(const FormLanguage::BilinearFormIntegratorSum& lsum);

         const TrialFunction<TrialFEC, Traits::Serial>& getTrialFunction() const override
         {
            return m_u;
         }

         const TestFunction<TestFEC, Traits::Serial>& getTestFunction() const override
         {
            return m_v;
         }

         BilinearForm& add(const BilinearFormIntegratorBase& bfi) override;

         BilinearForm& add(const FormLanguage::BilinearFormIntegratorSum& lsum) override;

         BilinearForm& from(const BilinearFormIntegratorBase& bfi) override;

         BilinearForm& from(const FormLanguage::BilinearFormIntegratorSum& bfi) override;

         mfem::BilinearForm& getHandle() override
         {
            return *m_bf;
         }

         const mfem::BilinearForm& getHandle() const override
         {
            return *m_bf;
         }

      private:
         TrialFunction<TrialFEC, Traits::Serial>& m_u;
         TestFunction<TestFEC, Traits::Serial>&   m_v;
         std::unique_ptr<mfem::BilinearForm> m_bf;
         BFIList m_bfiDomainList;
         BFIList m_bfiBoundaryList;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_domAttrMarkers;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_bdrAttrMarkers;
   };
   template <class TrialFEC, class TestFEC, class Trait>
   BilinearForm(TrialFunction<TrialFEC, Trait>&, TestFunction<TestFEC, Trait>&)
      -> BilinearForm<TrialFEC, TestFEC, Trait>;
}

#include "BilinearForm.hpp"

#endif

