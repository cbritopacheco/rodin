#ifndef RODIN_VARIATIONAL_VECTORDOMAINLFDIVINTEGRATOR_H
#define RODIN_VARIATIONAL_VECTORDOMAINLFDIVINTEGRATOR_H

#include "ForwardDecls.h"

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
               m_gshape.SetSize(dof, dim);
               m_Jadj.SetSize(dim);
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
                  el.CalcDShape(ip, m_dshape);
                  Tr.SetIntPoint(&ip);
                  double val = ip.weight * m_q.Eval(Tr, ip);
                  CalcAdjugate(Tr.Jacobian(), m_Jadj);
                  Mult(m_dshape, m_Jadj, m_gshape);
                  m_gshape.GradToDiv(m_divshape);
                  elvect.Add(val, m_divshape);
               }
            }

         private:
            mfem::Vector m_divshape;
            mfem::DenseMatrix m_dshape;
            mfem::DenseMatrix m_gshape;
            mfem::DenseMatrix m_Jadj;
            mfem::Coefficient& m_q;
      };
   }

   class VectorDomainLFDivIntegrator
      : public LinearFormDomainIntegrator
   {
      public:
         VectorDomainLFDivIntegrator(const ScalarCoefficientBase& f);

         VectorDomainLFDivIntegrator(const VectorDomainLFDivIntegrator& other);

         const std::set<int>& getAttributes() const override
         {
            return m_attr;
         }

         VectorDomainLFDivIntegrator& over(int attr) override
         {
            return over(std::set{attr});
         }

         VectorDomainLFDivIntegrator& over(const std::set<int>& attrs) override
         {
            m_attr = attrs;
            return *this;
         }

         void buildMFEMLinearFormIntegrator() override;

         mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() override;
         mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() override;

         VectorDomainLFDivIntegrator* copy() const noexcept override
         {
            return new VectorDomainLFDivIntegrator(*this);
         }
      private:
         std::set<int> m_attr;
         std::unique_ptr<ScalarCoefficientBase> m_f;
         std::unique_ptr<Internal::VectorDomainLFDivIntegrator> m_mfemLFI;
   };
}

#endif

