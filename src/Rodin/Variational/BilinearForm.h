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
    * @internal
    * @brief Abstract base class for objects representing bilinear forms.
    */
   class BilinearFormBase
   {
      public:
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
          * @brief Builds the bilinear form the given BilinearFormDomainIntegrator.
          * @param[in] bfi BilinearFormDomainIntegrator which will be used to
          * build the bilinear form.
          */
         virtual BilinearFormBase& from(
               const BilinearFormIntegratorBase& bfi) = 0;

         virtual BilinearFormBase& from(
               const FormLanguage::BilinearFormIntegratorSum& bfi) = 0;

         /**
          * @brief Adds a BilinearFormDomainIntegrator to the bilinear form.
          */
         virtual BilinearFormBase& add(
               const BilinearFormIntegratorBase& bfi) = 0;

         virtual BilinearFormBase& add(
               const FormLanguage::BilinearFormIntegratorSum& lsum) = 0;

         /**
          * @brief Reflects the changes in the mesh.
          */
         virtual void update() = 0;

         /**
          * @brief Assembles the bilinear form.
          */
         virtual void assemble() = 0;

         virtual const ShapeFunctionBase<Trial>& getTrialFunction() const = 0;

         virtual const ShapeFunctionBase<Test>& getTestFunction() const = 0;
   };

   /**
    * @brief Represents a bilinear form on some finite element space
    *
    * An object of type BilinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    a : U \times V &\rightarrow \mathbb{R}\\
    *        (u, v) &\mapsto a(u, v)
    * \end{aligned}
    * @f]
    * where @f$ U @f$ and @f$ V @f$ are finite element spaces. A bilinear form
    * can be specified by from one or more bilinear integrators (e.g.
    * DiffusionIntegrator).
    */
   template <class TrialFES, class TestFES>
   class BilinearForm : public BilinearFormBase
   {
      static_assert(
            std::is_same_v<TrialFES, TestFES>,
            "Different trial and test spaces are currently not supported.");

      public:
         using BFIList = std::vector<std::unique_ptr<BilinearFormIntegratorBase>>;

         /**
          * @brief Constructs a bilinear form defined on some finite element
          * space.
          *
          * @param[in] fes Reference to the finite element space
          */
         BilinearForm(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v);

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
               const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const;

         /**
          * @brief Builds the bilinear form from a derived instance of
          * BilinearFormIntegratorBase.
          * @param[in] bfi Bilinear form integrator
          */
         BilinearForm& operator=(const BilinearFormIntegratorBase& bfi);

         BilinearForm& operator=(const FormLanguage::BilinearFormIntegratorSum& lsum);

         void assemble() override;

         void update() override;

         const TrialFunction<TrialFES>& getTrialFunction() const override
         {
            return m_u;
         }

         const TestFunction<TestFES>& getTestFunction() const override
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
         TrialFunction<TrialFES>& m_u;
         TestFunction<TestFES>&   m_v;
         std::unique_ptr<mfem::BilinearForm> m_bf;
         BFIList m_bfiDomainList;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_domAttrMarkers;
   };
}

#include "BilinearForm.hpp"

#endif

