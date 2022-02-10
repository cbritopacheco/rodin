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
#include "LinearFormIntegrator.h"
#include "FormLanguage/ProblemBody.h"
#include "FormLanguage/LinearFormIntegratorSum.h"

namespace Rodin::Variational
{
   class LinearFormBase
   {
      public:
         virtual mfem::LinearForm& getHandle() = 0;
         virtual const mfem::LinearForm& getHandle() const = 0;

         virtual LinearFormBase& from(const LinearFormIntegratorBase& lfi) = 0;
         virtual LinearFormBase& from(const FormLanguage::LinearFormIntegratorSum& lsum) = 0;

         virtual LinearFormBase& add(const LinearFormIntegratorBase& lfi) = 0;
         virtual LinearFormBase& add(const FormLanguage::LinearFormIntegratorSum& lsum) = 0;

         virtual void assemble() = 0;
         virtual void update() = 0;
   };

   /**
    * @brief Represents a linear form on some finite element space
    * @tparam FEC Finite element collection on which the linear form operates
    *
    * An object of type LinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    L : V &\rightarrow \mathbb{R}\\
    *        v &\mapsto L(v)
    * \end{aligned}
    * @f]
    * where @f$ V @f$ is a finite element space. A linear form can be specified
    * by from one or more linear integrators (e.g. DomainLFIntegrator).
    */
   template <class FEC>
   class LinearForm : public LinearFormBase
   {
      public:
         using LFIList = std::vector<std::unique_ptr<LinearFormIntegratorBase>>;

         /**
          * @brief Constructs a linear form defined on some finite element
          * space
          * @param[in] fes Reference to the finite element space
          */
         LinearForm(FiniteElementSpace<FEC>& fes);

         LinearForm<FEC>& operator=(const LinearFormIntegratorBase& lfi);
         LinearForm<FEC>& operator=(const FormLanguage::LinearFormIntegratorSum& lsum);

         /**
          * @brief Evaluates the linear form at the function @f$ u @f$.
          *
          * Given a grid function @f$ u @f$, this function will compute the
          * action of the linear mapping @f$ L(u) @f$.
          *
          * @returns The value which the linear form takes at @f$ u @f$.
          */
         double operator()(const GridFunction<FEC>& u) const;

         void update() override;

         void assemble() override;

         LinearForm<FEC>& add(const LinearFormIntegratorBase& lfi) override;
         LinearForm<FEC>& add(const FormLanguage::LinearFormIntegratorSum& lsum) override;

         LinearForm<FEC>& from(const LinearFormIntegratorBase& lfi) override;
         LinearForm<FEC>& from(const FormLanguage::LinearFormIntegratorSum& lsum) override;

         mfem::LinearForm& getHandle() override
         {
            return *m_lf;
         }

         const mfem::LinearForm& getHandle() const override
         {
            return *m_lf;
         }

      private:
         FiniteElementSpace<FEC>& m_fes;
         std::unique_ptr<mfem::LinearForm> m_lf;
         LFIList m_lfiDomainList;
         LFIList m_lfiBoundaryList;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_bdrAttrMarkers;
         std::vector<std::unique_ptr<mfem::Array<int>>> m_domAttrMarkers;
   };
}

#include "LinearForm.hpp"

#endif
