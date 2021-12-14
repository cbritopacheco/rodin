/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
         virtual Transpose T() const;

         virtual int getRows() const = 0;
         virtual int getColumns() const = 0;
         virtual void buildMFEMMatrixCoefficient() = 0;
         virtual mfem::MatrixCoefficient& getMFEMMatrixCoefficient() = 0;
         virtual MatrixCoefficientBase* copy() const noexcept override = 0;
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

         int getRows() const override
         {
            return m_dimension;
         }

         int getColumns() const override
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

         virtual ScalarCoefficientBase& operator()(int i, int j)
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
