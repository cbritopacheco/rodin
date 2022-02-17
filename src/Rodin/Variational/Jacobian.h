/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class VectorGradientCoefficient : public mfem::MatrixCoefficient
      {
         public:
            VectorGradientCoefficient(mfem::GridFunction& u)
               :  mfem::MatrixCoefficient(
                     u.FESpace()->GetVDim(), u.FESpace()->GetMesh()->Dimension()),
                  m_u(u)
            {}

            VectorGradientCoefficient(mfem::GridFunction& u, int traceDomain)
               : mfem::MatrixCoefficient(
                     u.FESpace()->GetVDim(), u.FESpace()->GetMesh()->Dimension()),
                  m_u(u),
                  m_traceDomain(traceDomain)
            {}

            virtual void Eval(
                  mfem::DenseMatrix& grad,
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT
                     && T.mesh->FaceIsInterior((T.mesh->GetBdrFace(T.ElementNo))))
               {
                  if (!m_traceDomain)
                  {
                     Alert::Exception()
                        << "Integration over an interior boundary requires "
                        << "specifying a trace domain via setTraceDomain(int)"
                        << Alert::Raise;
                  }

                  int faceId = m_u.FESpace()->GetMesh()->GetBdrFace(T.ElementNo);
                  mfem::FaceElementTransformations* ft =
                     m_u.FESpace()->GetMesh()->GetFaceElementTransformations(faceId);
                  ft->SetAllIntPoints(&ip);
                  if (ft->GetElement1Transformation().Attribute == *m_traceDomain)
                     m_u.GetVectorGradient(ft->GetElement1Transformation(), grad);
                  else if (ft->GetElement2Transformation().Attribute == *m_traceDomain)
                     m_u.GetVectorGradient(ft->GetElement2Transformation(), grad);
                  else
                  {
                     // The boundary over which we are evaluating must be the
                     // interface between the trace domain and some other
                     // domain!
                     Alert::Exception()
                        << "Invalid boundary for trace domain " << *m_traceDomain
                        << Alert::Raise;
                  }
               }
               else
                  m_u.GetVectorGradient(T, grad);
            }

         private:
            mfem::GridFunction& m_u;
            std::optional<int> m_traceDomain;
      };
   }

   /**
    * @brief Represents the Jacobian matrix @f$ \mathbf{J}_u @f$ of the
    * function @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^n \rightarrow \mathbb{R}^m @f$, the Jacobian matrix
    * @f$ \mathbf{J}_u(x) @f$ at any point @f$ x = (x_1, \ldots, x_n) @f$ is
    * defined by the @f$ m \times n @f$ matrix
    * @f[
    * \mathbf{J}_u = \begin{bmatrix}
    * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_1}{\partial x_n}\\
    * \vdots & \ddots & \vdots\\
    * \dfrac{\partial u_m}{\partial x_1} & \ldots & \dfrac{\partial u_m}{\partial x_n}
    * \end{bmatrix}
    * @f]
    * 
    */
   class Jacobian : public MatrixCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         constexpr
         Jacobian(GridFunction<H1>& u)
            :  m_u(u)
         {}

         constexpr
         Jacobian(const Jacobian& other)
            : m_u(other.m_u),
              m_traceDomain(other.m_traceDomain)
         {}

         Jacobian& setTraceDomain(int domain);

         int getRows() const override;

         int getColumns() const override;

         void build() override;

         mfem::MatrixCoefficient& get() override;

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1>& m_u;
         std::optional<Internal::VectorGradientCoefficient> m_mfemMatrixCoefficient;
         std::optional<int> m_traceDomain;
   };
}

#endif
