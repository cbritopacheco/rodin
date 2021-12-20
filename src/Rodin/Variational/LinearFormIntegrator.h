#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class LinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         virtual void buildMFEMLinearFormIntegrator() = 0;

         virtual mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() = 0;

         /**
          * @brief Releases ownership of the mfem::LinearFormIntegrator.
          * @note After this call, calling @see getMFEMLinearFormIntegrator()
          * will result in undefined behaviour.
          * @warning The LinearFormIntegratorBase instance must still be kept
          * in memory since it might contain objects which the
          * mfem::LinearFormIntegrator instance refers to.
          */
         virtual mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() = 0;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;
   };
}

#endif
