/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATORSUM_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATORSUM_H

#include <memory>
#include <utility>

#include "Rodin/Variational/LinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class LinearFormIntegratorSum : public FormLanguage::Base
   {
      public:
         using LFIList = std::vector<std::unique_ptr<LinearFormIntegratorBase>>;

         LinearFormIntegratorSum(
               const LinearFormIntegratorBase& lhs,
               const LinearFormIntegratorBase& rhs);

         LinearFormIntegratorSum(
               const LinearFormIntegratorSum& lhs,
               const LinearFormIntegratorBase& rhs);

         LinearFormIntegratorSum(
               const LinearFormIntegratorSum& lhs,
               const LinearFormIntegratorSum& rhs);

         LinearFormIntegratorSum(LinearFormIntegratorSum&& other)
            :  FormLanguage::Base(std::move(other)),
               m_lfiDomainList(std::move(other.m_lfiDomainList)),
               m_lfiBoundaryList(std::move(other.m_lfiBoundaryList))
         {}

         LinearFormIntegratorSum(const LinearFormIntegratorSum& other);

         LFIList& getLinearFormDomainIntegratorList()
         {
            return m_lfiDomainList;
         }

         LFIList& getLinearFormBoundaryIntegratorList()
         {
            return m_lfiBoundaryList;
         }

         const LFIList& getLinearFormDomainIntegratorList() const
         {
            return m_lfiDomainList;
         }

         const LFIList& getLinearFormBoundaryIntegratorList() const
         {
            return m_lfiBoundaryList;
         }

         LinearFormIntegratorSum* copy() const noexcept override
         {
            return new LinearFormIntegratorSum(*this);
         }

      private:
         std::vector<std::unique_ptr<LinearFormIntegratorBase>> m_lfiDomainList;
         std::vector<std::unique_ptr<LinearFormIntegratorBase>> m_lfiBoundaryList;
   };

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs);

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorBase& rhs);

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorSum& rhs);

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorSum& rhs);

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs);

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorBase& rhs);

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorSum& rhs);

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorSum& lhs, const LinearFormIntegratorSum& rhs);

}

#endif

