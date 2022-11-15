/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include <mfem.hpp>

#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   /**
    * @brief Abstract base class for objects of type BilinearForm.
    */
   template <class OperatorType>
   class BilinearFormBase : public FormLanguage::Base
   {
      public:
         constexpr
         BilinearFormBase() = default;

         constexpr
         BilinearFormBase(const BilinearFormBase& other)
            :  FormLanguage::Base(other),
               m_bfis(other.m_bfis)
         {}

         constexpr
         BilinearFormBase(BilinearFormBase&& other)
            :  FormLanguage::Base(std::move(other)),
               m_bfis(std::move(other.m_bfis))
         {}

         constexpr
         const FormLanguage::List<BilinearFormIntegratorBase>& getIntegrators() const
         {
            return m_bfis;
         }

         /**
          * @brief Updates the state after a refinement in the mesh.
          *
          * This method will update the bilinear form after a call to the
          * @ref MeshBase::refine() "refine()" method.
          */
         virtual BilinearFormBase& update() = 0;

         /**
          * @brief Assembles the bilinear form.
          *
          * This method will assemble the underlying sparse matrix associated
          * the bilinear form.
          *
          * @see getMatrix()
          */
         virtual void assemble() = 0;

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Reference to the associated sparse matrix.
          */
         virtual OperatorType& getOperator() = 0;

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Constant reference to the associated sparse matrix.
          */
         virtual const OperatorType& getOperator() const = 0;

         /**
          * @brief Builds the bilinear form the given bilinear integrator
          * @param[in] bfi Bilinear integrator which will be used to
          * build the bilinear form.
          * @returns Reference to this (for method chaining)
          */
         virtual BilinearFormBase& from(const BilinearFormIntegratorBase& bfi)
         {
            m_bfis.clear();
            add(bfi).assemble();
            return *this;
         }

         virtual BilinearFormBase& from(
               const FormLanguage::List<BilinearFormIntegratorBase>& bfi)
         {
            m_bfis.clear();
            add(bfi).assemble();
            return *this;
         }

         virtual BilinearFormBase& operator=(const BilinearFormIntegratorBase& bfi)
         {
            from(bfi).assemble();
            return *this;
         }

         /**
          * @todo
          */
         virtual BilinearFormBase& operator=(
               const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
         {
            from(bfis).assemble();
            return *this;
         }

         /**
          * @brief Adds a bilinear integrator to the bilinear form.
          * @returns Reference to this (for method chaining)
          */
         virtual BilinearFormBase& add(const BilinearFormIntegratorBase& bfi)
         {
            m_bfis.add(bfi);
            return *this;
         }

         virtual BilinearFormBase& add(const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
         {
            m_bfis.add(bfis);
            return *this;
         }

         /**
          * @brief Gets the reference to the associated TrialFunction object.
          * @returns Reference to this (for method chaining)
          */
         virtual const ShapeFunctionBase<TrialSpace>& getTrialFunction() const = 0;

         /**
          * @brief Gets the reference to the associated TestFunction object.
          * @returns Reference to this (for method chaining)
          */
         virtual const ShapeFunctionBase<TestSpace>& getTestFunction() const = 0;

         virtual BilinearFormBase* copy() const noexcept override = 0;

      private:
         FormLanguage::List<BilinearFormIntegratorBase> m_bfis;
   };

   /**
    * @brief Represents a serial bilinear form supported on two finite element
    * spaces originating from two instances of FiniteElementCollection.
    * @tparam TrialFES Trial FiniteElementCollection
    * @tparam TestFES Test FiniteElementCollection
    *
    * An object of type BilinearForm represents a linear map
    * @f[
    * \begin{aligned}
    *    a : U \times V &\rightarrow \mathbb{R}\\
    *        (u, v) &\mapsto a(u, v)
    * \end{aligned}
    * @f]
    * where @f$ U @f$ and @f$ V @f$ are finite element spaces.
    *
    * A bilinear form may be built using the form language. For example,
    * @code{.cpp}
    * // Define spaces
    * FiniteElementSpace<H1> Vh;
    * TrialFunction u(Vh);
    * TestFunction  v(Vh);
    *
    * // Define mass bilinear form
    * BilinearForm bf(u, v);
    * bf = Integral(u, v);
    * @endcode
    */
   template <class TrialFES, class TestFES>
   class BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>
      : public BilinearFormBase<mfem::Operator>
   {
      static_assert(
            std::is_same_v<TrialFES, TestFES>,
            "Different trial and test spaces are currently not supported.");

      static_assert(std::is_same_v<typename TrialFES::Context, Context::Serial>);

      public:
         using Context = typename TrialFES::Context;
         using OperatorType = mfem::Operator;
         using Parent = BilinearFormBase<mfem::Operator>;

         /**
          * @brief Constructs a BilinearForm from a TrialFunction and
          * TestFunction.
          *
          * @param[in] u Trial function argument
          * @param[in] v Test function argument
          */
         constexpr
         BilinearForm(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
            :  m_u(u), m_v(v)
         {}

         constexpr
         BilinearForm(const BilinearForm& other)
            :  BilinearFormBase(other),
               m_u(other.m_u), m_v(other.m_v)
         {}

         constexpr
         BilinearForm(BilinearForm&& other)
            :  BilinearFormBase(std::move(other)),
               m_u(other.m_u), m_v(other.m_v),
               m_bf(std::move(other.m_bf))
         {}

         /**
          * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
          * v @f$.
          *
          * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
          * compute the action of the bilinear mapping @f$ a(u, v) @f$.
          *
          * @returns The action @f$ a(u, v) @f$ which the bilinear form takes
          * at @f$ ( u, v ) @f$.
          */
         constexpr
         double operator()(
               const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const;

         constexpr
         void formLinearSystem(
               const mfem::Array<int>& essTrueDofs,
               mfem::Vector& x,
               mfem::Vector& b,
               mfem::SparseMatrix& stiffness,
               mfem::Vector& guess,
               mfem::Vector& mass)
         {
            m_bf->FormLinearSystem(essTrueDofs, x, b, stiffness, guess, mass);
         }

         void assemble() override;

         BilinearForm& update() override;

         const TrialFunction<TrialFES>& getTrialFunction() const override
         {
            return m_u;
         }

         const TestFunction<TestFES>& getTestFunction() const override
         {
            return m_v;
         }

         BilinearForm& operator=(const BilinearFormIntegratorBase& bfi) override
         {
            from(bfi).assemble();
            return *this;
         }

         /**
          * @todo
          */
         BilinearForm& operator=(
               const FormLanguage::List<BilinearFormIntegratorBase>& bfis) override
         {
            from(bfis).assemble();
            return *this;
         }

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Reference to the associated sparse matrix.
          */
         virtual OperatorType& getOperator() override
         {
            assert(m_bf->HasSpMat());
            return m_bf->SpMat();
         }

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Constant reference to the associated sparse matrix.
          */
         virtual const OperatorType& getOperator() const override
         {
            assert(m_bf->HasSpMat());
            return m_bf->SpMat();
         }

         virtual BilinearForm* copy() const noexcept override
         {
            return new BilinearForm(*this);
         }

      private:
         TrialFunction<TrialFES>& m_u;
         TestFunction<TestFES>&   m_v;

         std::unique_ptr<mfem::BilinearForm> m_bf;
   };
   template <class TrialFES, class TestFES>
   BilinearForm(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
      -> BilinearForm<TrialFES, TestFES, typename TrialFES::Context, mfem::Operator>;

   template <class TrialFES, class TestFES>
   class BilinearForm<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix>
      : public BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>
   {
      public:
         using Context = typename TrialFES::Context;
         using OperatorType = mfem::SparseMatrix;
         using Parent = BilinearForm<TrialFES, TestFES, Context, mfem::Operator>;

         constexpr
         BilinearForm(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
            : Parent(u, v)
         {}

         constexpr
         BilinearForm(const BilinearForm& other)
            : Parent(other)
         {}

         constexpr
         BilinearForm(BilinearForm&& other)
            : Parent(std::move(other))
         {}

         virtual OperatorType& getOperator() override
         {
            assert(dynamic_cast<mfem::SparseMatrix*>(&Parent::getOperator()));
            return static_cast<mfem::SparseMatrix&>(Parent::getOperator());
         }

         /**
          * @brief Gets the reference to the (local) associated sparse matrix
          * to the bilinear form.
          * @returns Constant reference to the associated sparse matrix.
          */
         virtual const OperatorType& getOperator() const override
         {
            assert(dynamic_cast<const mfem::SparseMatrix*>(&Parent::getOperator()));
            return static_cast<const mfem::SparseMatrix&>(Parent::getOperator());
         }
   };
}

#include "BilinearForm.hpp"

#endif

