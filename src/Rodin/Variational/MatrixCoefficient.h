#ifndef RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H
#define RODIN_VARIATIONAL_MATRIXCOEFFICIENT_H

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/Base.h"

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class MatrixCoefficientBase : public FormLanguage::Base
   {
      public:
         virtual ~MatrixCoefficientBase() = default;

         virtual int getDimension() const = 0;
         virtual void buildMFEMMatrixCoefficient() = 0;
         virtual mfem::MatrixCoefficient& getMFEMMatrixCoefficient() = 0;
         virtual MatrixCoefficientBase* copy() const noexcept override = 0;
         virtual ScalarCoefficientBase& operator()(int i, int j) = 0;
   };

   class MatrixCoefficient : public MatrixCoefficientBase
   {
      public:
         MatrixCoefficient(int dimension)
            :  m_dimension(dimension),
               m_mfemMatrixCoefficient(dimension)
         {}

         MatrixCoefficient(const MatrixCoefficient& other)
            :  m_dimension(other.m_dimension),
               m_mfemMatrixCoefficient(other.m_mfemMatrixCoefficient)
         {}

         virtual ~MatrixCoefficient() = default;

         int getDimension() const override
         {
            return m_dimension;
         }

         virtual void buildMFEMMatrixCoefficient() override
         {
            for (int i = 0; i < m_dimension; i++)
            {
               for (int j = 0; j < m_dimension; j++)
               {
                  auto& s = m_scalarCoefficients[i * m_dimension + j];
                  s->buildMFEMCoefficient();
                  m_mfemMatrixCoefficient.Set(
                        i, j, &s->getMFEMCoefficient(),
                        false // Do not own the data
                        );
               }
            }
         }

         virtual mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override
         {
            return m_mfemMatrixCoefficient;
         }

         virtual MatrixCoefficient* copy() const noexcept override
         {
            return new MatrixCoefficient(*this);
         }

         virtual ScalarCoefficientBase& operator()(int i, int j) override
         {
            return *m_scalarCoefficients[i * m_dimension + j];
         }

      private:
         int m_dimension;
         std::vector<std::unique_ptr<ScalarCoefficientBase>> m_scalarCoefficients;
         mfem::MatrixArrayCoefficient m_mfemMatrixCoefficient;
   };
}

#endif
