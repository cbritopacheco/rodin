#ifndef RODIN_VARIATIONAL_VECTORDOMAINLFDIVINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORDOMAINLFDIVINTEGRATOR_H

#include "ForwardDecls.h"

#include "ScalarCoefficient.h"
#include "LinearFormIntegrator.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class VectorDomainLFDivIntegrator : public mfem::LinearFormIntegrator
      {
         public:
            VectorDomainLFDivIntegrator(mfem::Coefficient& q)
               : m_q(q)
            {}

            void AssembleRHSElementVect(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& Tr,
                  mfem::Vector& elvect) override
            {
               int dim  = el.GetDim();
               int dof = el.GetDof();

               m_dshape.SetSize(dof, dim);
               m_divshape.SetSize(dof * dim);

               elvect.SetSize(dof * dim);
               elvect = 0.0;

               const mfem::IntegrationRule *ir = IntRule ? IntRule : GetIntRule();
               if (ir == nullptr)
               {
                  int intorder = 2 * el.GetOrder();
                  ir = &mfem::IntRules.Get(el.GetGeomType(), intorder);
               }

               for (int i = 0; i < ir->GetNPoints(); i++)
               {
                  const mfem::IntegrationPoint &ip = ir->IntPoint(i);
                  Tr.SetIntPoint(&ip);
                  el.CalcPhysDShape(Tr, m_dshape);
                  double val = ip.weight * Tr.Weight() * m_q.Eval(Tr, ip);
                  m_dshape.GradToDiv(m_divshape);
                  elvect.Add(val, m_divshape);
               }
            }

         private:
            mfem::Vector m_divshape;
            mfem::DenseMatrix m_dshape;
            mfem::Coefficient& m_q;
      };
   }

   class VectorDomainLFDivIntegrator
      : public LinearFormDomainIntegrator
   {
      public:
         VectorDomainLFDivIntegrator(const ScalarCoefficientBase& f)
            : m_f(f.copy()),
              m_mfemScalar(m_f->build()),
              m_mfemLFI(*m_mfemScalar)
         {}

         VectorDomainLFDivIntegrator(const VectorDomainLFDivIntegrator& other)
            :  LinearFormDomainIntegrator(other),
               m_f(other.m_f->copy()),
               m_mfemScalar(m_f->build()),
               m_mfemLFI(*m_mfemScalar)
         {}

         void getElementVector(
                  const mfem::FiniteElement& el,
                  mfem::ElementTransformation& trans,
                  mfem::Vector& vec) override
         {
            m_mfemLFI.AssembleRHSElementVect(el, trans, vec);
         }

         VectorDomainLFDivIntegrator* copy() const noexcept override
         {
            return new VectorDomainLFDivIntegrator(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_f;
         std::unique_ptr<mfem::Coefficient> m_mfemScalar;
         Internal::VectorDomainLFDivIntegrator m_mfemLFI;
   };
}

#endif

