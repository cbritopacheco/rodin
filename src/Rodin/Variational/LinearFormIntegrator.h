#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <vector>
#include <optional>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class LinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         virtual ~LinearFormIntegratorBase() = default;

         virtual LinearFormIntegratorBase& over(int attr)
         {
            return over(std::vector<int>{attr});
         }

         virtual LinearFormIntegratorBase& over(const std::vector<int>& attrs)
         {
            m_attrs = attrs;
            return *this;
         }

         const std::vector<int>& getAttributes() const
         {
            return m_attrs;
         }

         virtual void buildMFEMLinearFormIntegrator() = 0;

         virtual mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() = 0;

         /**
          * @internal
          * @brief Releases ownership of the mfem::LinearFormIntegrator.
          *
          * @note After this call, calling getMFEMLinearFormIntegrator() will
          * result in undefined behaviour.
          *
          * @warning The LinearFormIntegratorBase instance must still be kept
          * in memory since it might contain objects which the
          * mfem::LinearFormIntegrator instance refers to.
          */
         virtual mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() = 0;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;

      private:
         std::vector<int> m_attrs;
   };

   class LinearFormDomainIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormDomainIntegrator& over(int attr) override
         {
            return LinearFormDomainIntegrator::over(std::vector<int>{attr});
         }

         LinearFormDomainIntegrator& over(const std::vector<int>& attrs) override
         {
            return static_cast<LinearFormDomainIntegrator&>(
                  LinearFormIntegratorBase::over(attrs));
         }

         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
   };

   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormBoundaryIntegrator() = default;

         LinearFormBoundaryIntegrator(
               const LinearFormBoundaryIntegrator&) = default;

         LinearFormBoundaryIntegrator(LinearFormBoundaryIntegrator&&) = default;

         LinearFormBoundaryIntegrator& over(int attr) override
         {
            return LinearFormBoundaryIntegrator::over(std::vector<int>{attr});
         }

         LinearFormBoundaryIntegrator& over(const std::vector<int>& attrs) override
         {
            return static_cast<LinearFormBoundaryIntegrator&>(
                  LinearFormIntegratorBase::over(attrs));
         }

         virtual ~LinearFormBoundaryIntegrator() = default;

         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;
   };

}

#endif
