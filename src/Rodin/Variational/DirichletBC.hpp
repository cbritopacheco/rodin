/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIRICHLETBC_HPP
#define RODIN_VARIATIONAL_DIRICHLETBC_HPP

#include "Rodin/Alert.h"

#include "DirichletBC.h"

namespace Rodin::Variational
{
   template <class T>
   DirichletBC<T>::DirichletBC(int bdrAttr, const T& value)
      : m_bdrAttr(bdrAttr), m_value(ScalarCoefficient<T>(value))
   {}

   template <class T>
   DirichletBC<T>::DirichletBC(const DirichletBC& other)
      :  m_bdrAttr(other.m_bdrAttr),
         m_value(other.m_value),
         m_problem(other.m_problem)
   {}

   template <class T>
   DirichletBC<T>& DirichletBC<T>::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   template <class T>
   void DirichletBC<T>::eval()
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
         Rodin::Alert::Exception("DirichletBC boundary attribute is out of range.")();

      // Project the coefficient onto the boundary
      m_essBdr = mfem::Array<int>(maxBdrAttr);
      m_essBdr = 0;
      m_essBdr[m_bdrAttr - 1] = 1;

      m_value.eval();
      m_problem->get().getSolution()
                      .getHandle()
                      .ProjectBdrCoefficient(m_value.coeff(), m_essBdr);

      // Keep track of the boundary attributes that have been projected
      m_problem->get().getEssentialBoundary()[m_bdrAttr - 1] = 1;
   }

   template <class T>
   template <class ... Args>
   DirichletBC<T>* DirichletBC<T>::create(Args&&... args) noexcept
   {
      return new DirichletBC(std::forward<Args>(args)...);
   }

   template <class T>
   DirichletBC<T>* DirichletBC<T>::copy() const noexcept
   {
      return new DirichletBC(*this);
   }

   // ---- VectorCoeff<T> ----------------------------------------------------
   template <class ... Values>
   DirichletBC<VectorCoefficient<Values...>>::DirichletBC(
         int bdrAttr, const VectorCoefficient<Values...>& value)
      : m_bdrAttr(bdrAttr), m_value(value)
   {}

   template <class ... Values>
   DirichletBC<VectorCoefficient<Values...>>::DirichletBC(const DirichletBC& other)
      :  m_bdrAttr(other.m_bdrAttr),
         m_value(other.m_value),
         m_problem(other.m_problem)
   {}

   template <class ... Values>
   DirichletBC<VectorCoefficient<Values...>>&
   DirichletBC<VectorCoefficient<Values...>>::setProblem(ProblemBase& problem)
   {
      m_problem.emplace(problem);
      return *this;
   }

   template <class ... Values>
   void DirichletBC<VectorCoefficient<Values...>>::eval()
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
         Rodin::Alert::Exception("DirichletBC boundary attribute is out of range.")();

      // Project the coefficient onto the boundary
      m_essBdr = mfem::Array<int>(maxBdrAttr);
      m_essBdr = 0;
      m_essBdr[m_bdrAttr - 1] = 1;

      m_value.eval();
      m_problem->get().getSolution().getHandle().ProjectBdrCoefficient(m_value.coeff(), m_essBdr);

      // Keep track of the boundary attributes that have been projected
      m_problem->get().getEssentialBoundary()[m_bdrAttr - 1] = 1;
   }

   template <class ... Values>
   template <class ... Args>
   DirichletBC<VectorCoefficient<Values...>>*
   DirichletBC<VectorCoefficient<Values...>>::create(Args&&... args) noexcept
   {
      return new DirichletBC(std::forward<Args>(args)...);
   }

   template <class ... Values>
   DirichletBC<VectorCoefficient<Values...>>*
   DirichletBC<VectorCoefficient<Values...>>::copy() const noexcept
   {
      return new DirichletBC(*this);
   }
}
#endif
