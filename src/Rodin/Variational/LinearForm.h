/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_H
#define RODIN_VARIATIONAL_LINEARFORM_H

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "TestFunction.h"
#include "LinearFormIntegrator.h"
#include "FormLanguage/LinearFormIntegratorSum.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects of type LinearForm.
    */
   class LinearFormBase
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
            return getHandle().Update();
         }

         /**
          * @brief Assembles the linear form.
          *
          * This method will assemble the underlying vector associated
          * the linear form.
          *
          * @see getVector()
          */
         void assemble()
         {
            return getHandle().Assemble();
         }

         /**
          * @brief Gets the reference to the (local) associated vector
          * to the LinearForm.
          */
         mfem::Vector& getVector()
         {
            return static_cast<mfem::Vector&>(getHandle());
         }

         /**
          * @brief Gets the reference to the (local) associated vector
          * to the LinearForm.
          */
         const mfem::Vector& getVector() const
         {
            return static_cast<const mfem::Vector&>(getHandle());
         }

         /**
          * @brief Builds the linear form the given LinearFormIntegratorBase
          * instance
          * @param[in] lfi Integrator which will be used to build the linear form.
          * @returns Reference to this (for method chaining)
          */
         virtual LinearFormBase& from(const LinearFormIntegratorBase& lfi) = 0;

         virtual LinearFormBase& from(const FormLanguage::LinearFormIntegratorSum& lsum) = 0;

         /**
          * @brief Builds the linear form the given LinearFormIntegratorBase
          * instance
          * @param[in] lfi Integrator which will be used to build the linear form.
          * @returns Reference to this (for method chaining)
          */
         virtual LinearFormBase& add(const LinearFormIntegratorBase& lfi) = 0;

         virtual LinearFormBase& add(const FormLanguage::LinearFormIntegratorSum& lsum) = 0;

         virtual const ShapeFunctionBase<Test>& getTestFunction() const = 0;

         virtual mfem::LinearForm& getHandle() = 0;

         virtual const mfem::LinearForm& getHandle() const = 0;
   };

   /**
    * @brief Represents a linear form defined over some finite element space
    *
    * An object of type LinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    L : V &\rightarrow \mathbb{R}\\
    *        v &\mapsto L(v)
    * \end{aligned}
    * @f]
    * where @f$ V @f$ is a finite element space.
    *
    * A linear form can be specified by from one or more
    * LinearFormIntegratorBase instances.
    */
   template <class FEC>
   class LinearForm<FEC, Traits::Serial> : public LinearFormBase
   {
      public:
         using LFIList = std::vector<std::unique_ptr<LinearFormIntegratorBase>>;

         /**
          * @brief Constructs a linear form defined on some finite element
          * space
          * @param[in] fes Reference to the finite element space
          */
         LinearForm(TestFunction<FEC, Traits::Serial>& v);

         LinearForm& operator=(const LinearFormIntegratorBase& lfi);

         LinearForm& operator=(const FormLanguage::LinearFormIntegratorSum& lsum);

         /**
          * @brief Evaluates the linear form at the function @f$ u @f$.
          *
          * Given a grid function @f$ u @f$, this function will compute the
          * action of the linear mapping @f$ L(u) @f$.
          *
          * @returns The value which the linear form takes at @f$ u @f$.
          */
         double operator()(const GridFunction<FEC, Traits::Serial>& u) const;

         LinearForm& add(const LinearFormIntegratorBase& lfi) override;

         LinearForm& add(const FormLanguage::LinearFormIntegratorSum& lsum) override;

         LinearForm& from(const LinearFormIntegratorBase& lfi) override;

         LinearForm& from(const FormLanguage::LinearFormIntegratorSum& lsum) override;

         const TestFunction<FEC, Traits::Serial>& getTestFunction() const override
         {
            return m_v;
         }

         mfem::LinearForm& getHandle() override
         {
            return *m_lf;
         }

         const mfem::LinearForm& getHandle() const override
         {
            return *m_lf;
         }

      private:
         TestFunction<FEC, Traits::Serial>& m_v;
         std::unique_ptr<mfem::LinearForm> m_lf;
         LFIList m_lfiDomainList;
         LFIList m_lfiBoundaryList;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_bdrAttrMarkers;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_domAttrMarkers;
   };
   template <class FEC, class Trait>
   LinearForm(TestFunction<FEC, Trait>&)
      -> LinearForm<FEC, Trait>;
}

#include "LinearForm.hpp"

#endif
