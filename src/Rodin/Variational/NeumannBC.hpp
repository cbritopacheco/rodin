/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_NEUMANNBC_HPP
#define RODIN_VARIATIONAL_NEUMANNBC_HPP

#include "NeumannBC.h"

namespace Rodin::Variational
{
   template <class T>
   NeumannBC<T>::NeumannBC(int bdrAttr, const T& value)
      : m_bdrAttr(bdrAttr), m_value(ScalarCoefficient<T>(value))
   {}

   template <class T>
   NeumannBC<T>::NeumannBC(const NeumannBC& other)
      :  m_bdrAttr(other.m_bdrAttr),
         m_value(other.m_value),
         m_problem(other.m_problem)
   {}

   template <class T>
   NeumannBC<T>& NeumannBC<T>::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   template <class T>
   void NeumannBC<T>::eval()
   {
      assert(m_problem.has_value());
      assert(m_problem->get().getSolution().getHandle().VectorDim() == 1);

      int maxBdrAttr = m_problem->get()
                                 .getSolution()
                                 .getHandle()
                                 .FESpace()
                                 ->GetMesh()
                                 ->bdr_attributes.Max();

      if (m_bdrAttr > maxBdrAttr)
         Rodin::Alert::Exception("NeumannBC boundary attribute is out of range.")();

      m_nbcBdr = mfem::Array<int>(maxBdrAttr);
      m_nbcBdr = 0;
      m_nbcBdr[m_bdrAttr - 1] = 1;

      m_value.eval();
      m_problem->get().getLinearForm().getHandle().AddBoundaryIntegrator(
            new mfem::BoundaryLFIntegrator(m_value.coeff()), m_nbcBdr);
   }

   template <class T>
   template <class ... Args>
   NeumannBC<T>* NeumannBC<T>::create(Args&&... args) noexcept
   {
      return new NeumannBC(std::forward<Args>(args)...);
   }

   template <class T>
   NeumannBC<T>* NeumannBC<T>::copy() const noexcept
   {
      return new NeumannBC(*this);
   }

   // ---- VectorCoeff<T> ----------------------------------------------------
   template <class ... Values>
   NeumannBC<VectorCoefficient<Values...>>::NeumannBC(
         int bdrAttr, const VectorCoefficient<Values...>& value)
      : m_bdrAttr(bdrAttr), m_value(value)
   {}

   template <class ... Values>
   NeumannBC<VectorCoefficient<Values...>>::NeumannBC(const NeumannBC& other)
      :  m_bdrAttr(other.m_bdrAttr),
         m_value(other.m_value),
         m_problem(other.m_problem)
   {}

   template <class ... Values>
   NeumannBC<VectorCoefficient<Values...>>&
   NeumannBC<VectorCoefficient<Values...>>::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   template <class ... Values>
   void NeumannBC<VectorCoefficient<Values...>>::eval()
   {
      assert(m_problem.has_value());
      assert(m_problem->get().getSolution().getHandle().VectorDim() == m_value.dimension());

      int maxBdrAttr = m_problem->get()
                                 .getSolution()
                                 .getHandle()
                                 .FESpace()
                                 ->GetMesh()
                                 ->bdr_attributes.Max();

      if (m_bdrAttr > maxBdrAttr)
         Rodin::Alert::Exception("NeumannBC boundary attribute is out of range.")();

      m_nbcBdr = mfem::Array<int>(maxBdrAttr);
      m_nbcBdr = 0;
      m_nbcBdr[m_bdrAttr - 1] = 1;

      m_value.eval();
      m_problem->get()
               .getLinearForm()
               .getHandle()
               .AddBoundaryIntegrator(
            new mfem::VectorBoundaryLFIntegrator(m_value.coeff()), m_nbcBdr);
   }

   template <class ... Values>
   template <class ... Args>
   NeumannBC<VectorCoefficient<Values...>>*
   NeumannBC<VectorCoefficient<Values...>>::create(Args&&... args) noexcept
   {
      return new NeumannBC(std::forward<Args>(args)...);
   }

   template <class ... Values>
   NeumannBC<VectorCoefficient<Values...>>*
   NeumannBC<VectorCoefficient<Values...>>::copy() const noexcept
   {
      return new NeumannBC(*this);
   }
}

#endif

