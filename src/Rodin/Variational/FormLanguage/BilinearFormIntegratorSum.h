/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATORSUM_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATORSUM_H

#include <memory>
#include <utility>

#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class BilinearFormIntegratorSum : public FormLanguage::Base
   {
      public:
         using BFIList = std::vector<std::unique_ptr<BilinearFormIntegratorBase>>;

         BilinearFormIntegratorSum(
               const BilinearFormIntegratorBase& lhs,
               const BilinearFormIntegratorBase& rhs);

         BilinearFormIntegratorSum(
               const BilinearFormIntegratorSum& lhs,
               const BilinearFormIntegratorBase& rhs);

         BilinearFormIntegratorSum(
               const BilinearFormIntegratorSum& lhs,
               const BilinearFormIntegratorSum& rhs);

         BilinearFormIntegratorSum(BilinearFormIntegratorSum&& other)
            :  FormLanguage::Base(std::move(other)),
               m_bfiDomainList(std::move(other.m_bfiDomainList)),
               m_bfiBoundaryList(std::move(other.m_bfiBoundaryList))
         {}


         BilinearFormIntegratorSum(const BilinearFormIntegratorSum& other);

         BFIList& getBilinearFormDomainIntegratorList()
         {
            return m_bfiDomainList;
         }

         const BFIList& getBilinearFormDomainIntegratorList() const
         {
            return m_bfiDomainList;
         }

         BFIList& getBilinearFormBoundaryIntegratorList()
         {
            return m_bfiBoundaryList;
         }

         const BFIList& getBilinearFormBoundaryIntegratorList() const
         {
            return m_bfiBoundaryList;
         }

         BilinearFormIntegratorSum* copy() const noexcept override
         {
            return new BilinearFormIntegratorSum(*this);
         }

      private:
         BFIList m_bfiDomainList;
         BFIList m_bfiBoundaryList;
   };

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorBase& rhs);

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorBase& rhs);

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorSum& rhs);

   BilinearFormIntegratorSum operator+(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorSum& rhs);

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorBase& rhs);

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorBase& rhs);

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorBase& lhs, const BilinearFormIntegratorSum& rhs);

   BilinearFormIntegratorSum operator-(
         const BilinearFormIntegratorSum& lhs, const BilinearFormIntegratorSum& rhs);

}

#endif

