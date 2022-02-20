#include "ScalarCoefficient.h"

#include "VectorBoundaryFluxLFIntegrator.h"

namespace Rodin::Variational
{
   VectorBoundaryFluxLFIntegrator::VectorBoundaryFluxLFIntegrator(
         const ScalarCoefficientBase& f)
      :  m_f(f.copy()),
         m_mfemScalar(m_f->build()),
         m_mfemLFI(*m_mfemScalar)
   {}

   VectorBoundaryFluxLFIntegrator
   ::VectorBoundaryFluxLFIntegrator(const VectorBoundaryFluxLFIntegrator& other)
      : LinearFormBoundaryIntegrator(other),
        m_f(other.m_f->copy()),
        m_mfemScalar(m_f->build()),
        m_mfemLFI(*m_mfemScalar)
   {}
}
