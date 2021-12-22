/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include <mfem.hpp>


#include "FormLanguage/List.h"

#include "ForwardDecls.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   class BilinearFormBase
   {
      public:
         virtual mfem::BilinearForm& getHandle() = 0;
         virtual const mfem::BilinearForm& getHandle() const = 0;

         virtual BilinearFormBase& from(
               const BilinearFormDomainIntegrator& bfi) = 0;

         virtual BilinearFormBase& add(
               const BilinearFormDomainIntegrator& bfi) = 0;

         virtual void assemble() = 0;
   };

   /**
    * @brief Represents a bilinear form on some finite element space
    * @tparam FEC Finite element collection on which the bilinear form operates
    *
    * An object of type BilinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    a : U \times V &\rightarrow \mathbb{R}\\
    *        (u, v) &\mapsto a(u, v)
    * \end{aligned}
    * @f]
    * where $U$ and $V$ are finite element spaces. A bilinear form can be
    * specified by from one or more bilinear integrators (e.g.
    * DiffusionIntegrator).
    */
   template <class FEC>
   class BilinearForm : public BilinearFormBase
   {
      public:
         /**
          * @brief Constructs a bilinear form defined on some finite element
          * space.
          *
          * @param[in] fes Reference to the finite element space
          */
         BilinearForm(FiniteElementSpace<FEC>& fes);

         /**
          * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
          * v @f$.
          *
          * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
          * compute the action of the bilinear mapping @f$ a(u, v) @f$.
          *
          * @returns The value which the bilinear form takes at
          * @f$ ( u, v ) @f$.
          */
         double operator()(
               const GridFunction<FEC>& u, const GridFunction<FEC>& v) const;

         template <class T>
         std::enable_if_t<
            std::is_base_of_v<BilinearFormIntegratorBase, T>, BilinearForm<FEC>&>
         operator=(const T& bfi)
         {
            from(bfi).assemble();
            return *this;
         }

         void assemble() override;

         BilinearForm<FEC>& add(const BilinearFormDomainIntegrator& bfi) override;

         BilinearForm<FEC>& from(const BilinearFormDomainIntegrator& bfi) override;

         mfem::BilinearForm& getHandle() override
         {
            return *m_bf;
         }

         const mfem::BilinearForm& getHandle() const override
         {
            return *m_bf;
         }

      private:
         FiniteElementSpace<FEC>& m_fes;
         std::unique_ptr<mfem::BilinearForm> m_bf;
         FormLanguage::List<BilinearFormDomainIntegrator> m_bfiDomainList;
   };
}

#include "BilinearForm.hpp"

#endif

