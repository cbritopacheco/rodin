/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <set>
#include <variant>
#include <functional>
#include <boost/mp11.hpp>

#include "Rodin/Pair.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/Solver/Solver.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/System.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/BlockSparseMatrix.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Tuple.h"
#include "Rodin/Utility/Extract.h"
#include "Rodin/Utility/Product.h"
#include "Rodin/Utility/Wrap.h"

#include "ForwardDecls.h"

#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"
#include "TrialFunction.h"
#include "TestFunction.h"


namespace Rodin::Variational
{
  /**
   * @defgroup ProblemSpecializations Problem Template Specializations
   * @brief Template specializations of the Problem class.
   * @see Problem
   */

  /**
   * @brief Abstract base class for variational problems.
   */
  template <class Operator, class Vector, class Scalar>
  class ProblemBase : public FormLanguage::Base
  {
    public:
      using VectorType = Vector;

      using OperatorType = Operator;

      using VectorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<Vector>>::ScalarType;

      using OperatorScalarType =
        typename FormLanguage::Traits<
          std::remove_reference_t<Operator>>::ScalarType;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ScalarType = Scalar;

      ProblemBase() = default;

      ProblemBase(ProblemBase&& other) = default;

      ProblemBase(const ProblemBase& other) = default;

      virtual ProblemBase& operator=(
          const ProblemBody<OperatorType, VectorType, ScalarType>& rhs) = 0;

      virtual void solve(
          Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) = 0;

      /**
       * @brief Assembles the underlying linear system to solve.
       */
      virtual ProblemBase& assemble() = 0;

      virtual LinearSystemType& getLinearSystem() = 0;

      virtual const LinearSystemType& getLinearSystem() const = 0;

      virtual ProblemBase* copy() const noexcept override = 0;
  };

  template <class TrialFES, class TestFES, class Operator, class Vector>
  class Problem<TrialFES, TestFES, Operator, Vector>;

  /**
   * @ingroup ProblemSpecializations
   * @brief General class to assemble linear systems with `Operator`
   * and `Vector` generic types in a sequential context.
   */
  template <class TrialFES, class TestFES, class Operator, class Vector>
  class Problem<TrialFES, TestFES, Operator, Vector>
    : public ProblemBase<Operator, Vector,
        typename FormLanguage::Mult<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>
  {
    public:
      using TrialFESType = TrialFES;

      using TestFESType = TestFES;

      using OperatorType = Operator;

      using VectorType = Vector;

      using ScalarType =
        typename FormLanguage::Mult<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using TrialFESScalarType =
        typename FormLanguage::Traits<TrialFES>::ScalarType;

      using TestFESScalarType =
        typename FormLanguage::Traits<TestFES>::ScalarType;

      using OperatorScalarType =
        typename FormLanguage::Mult<TrialFESScalarType, TestFESScalarType>::Type;

      using VectorScalarType = TestFESScalarType;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<TestFESScalarType>;

      using ContextType = Context::Local;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using Parent = ProblemBase<OperatorType, VectorType, ScalarType>;

      constexpr
      Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
        : Problem(u, v, LinearSystemType())
      {}

      /**
       * @brief Constructs an empty problem involving the trial function @f$ u @f$
       * and the test function @f$ v @f$.
       *
       * @param[in,out] u Trial function
       * @param[in,out] v %Test function
       */
      constexpr
      Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v, const LinearSystemType& axb)
         :  m_trialFunction(u), m_testFunction(v),
            m_axb(axb),
            m_assembled(false)
      {}

      /**
       * @brief Deleted copy constructor.
       */
      Problem(const Problem& other) = delete;

      /**
       * @brief Deleted copy assignment operator.
       */
      void operator=(const Problem& other) = delete;

      constexpr
      TrialFunction<TrialFES>& getTrialFunction()
      {
        return m_trialFunction;
      }

      constexpr
      TestFunction<TestFES>& getTestFunction()
      {
        return m_testFunction;
      }

      constexpr
      const TrialFunction<TrialFES>& getTrialFunction() const
      {
        return m_trialFunction.get();
      }

      constexpr
      const TestFunction<TestFES>& getTestFunction() const
      {
        return m_testFunction.get();
      }

      Problem& assemble() override
      {
        auto& pb = m_pb;
        auto& u = getTrialFunction();
        auto& v = getTestFunction();
        auto& axb = getLinearSystem();
        auto& mass = axb.getVector();
        auto& stiffness = axb.getOperator();
        auto& bfs = pb.getBFs();
        auto& dbcs = pb.getDBCs();
        auto& pbcs = pb.getPBCs();

        LinearForm lf(v, mass);
        for (auto& lfi : pb.getLFIs())
          lf.add(UnaryMinus(lfi)); // Negate every linear form integrator
        lf.assemble();

        BilinearForm bf(u, v, stiffness);
        for (auto& bfi : pb.getLocalBFIs())
          bf.add(bfi);
        for (auto& bfi : pb.getGlobalBFIs())
          bf.add(bfi);
        for (auto& bf : bfs)
        {
          bf.assemble();
          stiffness += bf.getOperator();
        }
        bf.assemble();

        // Impose Dirichlet boundary conditions
        auto& trial = getTrialFunction();
        const auto& trialFES = trial.getFiniteElementSpace();
        const auto& test = getTestFunction();
        const auto& testFES = test.getFiniteElementSpace();
        for (auto& dbc : dbcs)
        {
          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          if (dbc.isComponent())
          {
            assert(false);
          }
          else
          {
            m_axb.eliminate(dofs);
          }
        }

        // Impose periodic boundary conditions
        if (trialFES == testFES)
        {
          for (auto& pbc : pbcs)
          {
            pbc.assemble();
            const auto& dofs = pbc.getDOFs();

            if (pbc.isComponent())
            {
              assert(false);
            }
            else
            {
              m_axb.merge(dofs);
            }
          }
        }
        else
        {
          assert(false); // Not handled yet
        }

        m_assembled = true;

        return *this;
      }

      void solve(
          Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) override
      {
         // Assemble the system
         if (!m_assembled)
            assemble();

         // Solve the system AX = B
         auto& axb = getLinearSystem();
         auto& a = axb.getOperator();
         auto& x = axb.getGuess();
         auto& b = axb.getVector();
         solver.solve(a, x, b);

         // Recover solution
         getTrialFunction().emplace().getSolution().setWeights(x);
      }

      Problem& operator=(const ProblemBody<OperatorType, VectorType, ScalarType>& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      Problem* copy() const noexcept override
      {
        assert(false);
        return nullptr;
      }

    private:
      std::reference_wrapper<TrialFunction<TrialFES>> m_trialFunction;
      std::reference_wrapper<TestFunction<TestFES>>   m_testFunction;

      LinearSystemType  m_axb;
      Boolean           m_assembled;

      ProblemBody<OperatorType, VectorType, ScalarType> m_pb;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class TrialFES, class TestFES>
  Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
    -> Problem<
        TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Math::Vector<
          typename FormLanguage::Traits<TestFES>::ScalarType>>;

  /**
   * @ingroup RodinCTAD
   */
  template <class TrialFES, class TestFES, class Operator, class Vector>
  Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v, Math::LinearSystem<Operator, Vector>& axb)
    -> Problem<TrialFES, TestFES, Operator, Vector>;

  template <class Operator, class Vector, class U1, class U2, class ... Us>
  class Problem<
      Tuple<U1, U2, Us...>, Operator, Vector>
    : public ProblemBase<Operator, Vector, Real>
  {

    template <class T>
    struct IsTrialOrTestFunction
    {
      static constexpr Boolean Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
    };

    static_assert(Utility::ParameterPack<U1, U2, Us...>::template All<IsTrialOrTestFunction>::Value);

    public:
      using ScalarType = Real;

      using ContextType = Context::Local;

      using OperatorType = Operator;

      using VectorType = Vector;

      using Parent = ProblemBase<OperatorType, VectorType, Real>;

    private:
      template <class T>
      struct GetFES;

      template <class T>
      struct GetFES<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::FES;
      };

      template <class T>
      struct IsTrialFunctionReferenceWrapper
      {
        static constexpr Boolean Value = false;
      };

      template <class T>
      struct IsTrialFunctionReferenceWrapper<std::reference_wrapper<T>>
      {
        static constexpr Boolean Value = IsTrialFunction<T>::Value;
      };

      template <class T>
      struct IsTestFunctionReferenceWrapper
      {
        static constexpr Boolean Value = false;
      };

      template <class T>
      struct IsTestFunctionReferenceWrapper<std::reference_wrapper<T>>
      {
        static constexpr Boolean Value = IsTestFunction<T>::Value;
      };

      using TrialFunctionTuple =
        decltype(std::declval<
          Tuple<
            std::reference_wrapper<U1>,
            std::reference_wrapper<U2>,
            std::reference_wrapper<Us>...>>()
            .template filter<IsTrialFunctionReferenceWrapper>());

      using TestFunctionTuple =
        decltype(std::declval<
          Tuple<
            std::reference_wrapper<U1>,
            std::reference_wrapper<U2>,
            std::reference_wrapper<Us>...>>()
            .template filter<IsTestFunctionReferenceWrapper>());

      using TrialFESTuple = typename Utility::Extract<TrialFunctionTuple>::template Type<GetFES>;

      using TestFESTuple = typename Utility::Extract<TestFunctionTuple>::template Type<GetFES>;

      template <class TrialFES, class TestFES>
      using BilinearFormType = BilinearForm<TrialFES, TestFES, OperatorType>;

      template <class TestFES>
      using LinearFormType = LinearForm<TestFES, VectorType>;

      using BilinearFormTuple =
        typename Utility::Product<TrialFESTuple, TestFESTuple>::template Type<BilinearFormType>;

      using LinearFormTuple =
        typename Utility::Wrap<TestFESTuple>::template Type<LinearFormType>;

    public:
      using BilinearFormTupleSequentialAssembly =
        Assembly::Sequential<OperatorType, BilinearFormTuple>;

      using LinearFormTupleSequentialAssembly =
        Assembly::Sequential<VectorType, LinearFormTuple>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      Problem(U1& u1, U2& u2, Us&... us)
        : m_us(
            Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
            .template filter<IsTrialFunctionReferenceWrapper>()),
          m_vs(
            Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
            .template filter<IsTestFunctionReferenceWrapper>()),
          m_lft(m_vs.map(
                [](const auto& v)
                { return LinearFormType<
                    typename std::decay_t<
                    typename Utility::UnwrapRefDecay<decltype(v)>::Type>::FES>(v.get());
                })),
          m_bft(m_us.product(
                [](const auto& u, const auto& v) { return Pair(u, v); }, m_vs)
                    .map(
                      [](const auto& uv)
                      { return BilinearFormType<
                          typename std::decay_t<
                          typename Utility::UnwrapRefDecay<decltype(uv.first())>::Type>::FES,
                          typename std::decay_t<
                          typename Utility::UnwrapRefDecay<decltype(uv.second())>::Type>::FES>(
                              uv.first().get(), uv.second().get());
                      })),
          m_assembled(false),
          m_axb(m_stiffness, m_guess, m_mass)
      {
        m_bfa.reset(new BilinearFormTupleSequentialAssembly);
        m_lfa.reset(new LinearFormTupleSequentialAssembly);
        m_us.iapply([&](size_t i, const auto& u)
            { m_trialUUIDMap.right.insert({ i, u.get().getUUID() }); });
        m_vs.iapply([&](size_t i, const auto& v)
            { m_testUUIDMap.right.insert({ i, v.get().getUUID() }); });
      }

      Problem& assemble() override
      {
        auto& axb = m_axb;
        auto bt =
          m_bft.map(
              [](auto& bf)
              {
                auto& u = bf.getTrialFunction();
                auto& v = bf.getTestFunction();
                return Assembly::BilinearFormAssemblyInput(
                    u.getFiniteElementSpace(), v.getFiniteElementSpace(),
                    bf.getLocalIntegrators(), bf.getGlobalIntegrators());
              });

        auto lt =
          m_lft.map(
              [](auto& lf)
              {
                auto& v = lf.getTestFunction();
                return Assembly::LinearFormAssemblyInput(
                    v.getFiniteElementSpace(), lf.getIntegrators());
              });

        // Compute trial offsets
        {
          std::array<size_t, TrialFunctionTuple::Size> sz;
          m_us.map(
                [](const auto& u) { return u.get().getFiniteElementSpace().getSize(); })
              .iapply(
                [&](const Index i, size_t s) { sz[i] = s; });
          m_trialOffsets[0] = 0;
          for (size_t i = 0; i < TrialFunctionTuple::Size - 1; i++)
            m_trialOffsets[i + 1] = sz[i] + m_trialOffsets[i];
        }

        // Compute test offsets
        {
          std::array<size_t, TestFunctionTuple::Size> sz;
          m_vs.map(
                [](const auto& u) { return u.get().getFiniteElementSpace().getSize(); })
              .iapply(
                [&](const Index i, size_t s) { sz[i] = s; });
          m_testOffsets[0] = 0;
          for (size_t i = 0; i < TrialFunctionTuple::Size - 1; i++)
            m_testOffsets[i + 1] = sz[i] + m_testOffsets[i];
        }

        size_t rows =
          m_vs.reduce(
            [](const auto& a, const auto& b)
            {
              return a.get().getFiniteElementSpace().getSize() + b.get().getFiniteElementSpace().getSize();
            });

        size_t cols =
          m_us.reduce(
            [](const auto& a, const auto& b)
            {
              return a.get().getFiniteElementSpace().getSize() + b.get().getFiniteElementSpace().getSize();
            });

        // Compute block offsets to build the triplets
        std::array<Pair<size_t, size_t>, decltype(bt)::Size> boffsets;
        std::array<size_t, decltype(lt)::Size> loffsets;

        m_bft.iapply(
            [&](const Index i, const auto& bf)
            {
              auto ui = m_trialUUIDMap.left.find(bf.getTrialFunction().getUUID());
              auto vi = m_testUUIDMap.left.find(bf.getTestFunction().getUUID());
              if (ui != m_trialUUIDMap.left.end() && vi != m_testUUIDMap.left.end())
                boffsets[i] = Pair{ m_trialOffsets[ui->second], m_testOffsets[vi->second] };
            });

        m_lft.iapply(
            [&](const Index i, const auto& lf)
            {
              auto vi = m_testUUIDMap.left.find(lf.getTestFunction().getUUID());
              if (vi != m_testUUIDMap.left.end())
                loffsets[i] = m_testOffsets[vi->second];
            });

        // Assemble stiffness operator
        m_bfa->execute(m_stiffness,
          Assembly::BilinearFormTupleAssemblyInput(rows, cols, boffsets, bt));

        // Assemble mass vector
        m_lfa->execute(m_mass,
          Assembly::LinearFormTupleAssemblyInput(rows, loffsets, lt));

        // Impose Dirichlet boundary conditions
        m_us.apply(
            [&](const auto& u)
            {
              const auto ui = m_trialUUIDMap.left.find(u.get().getUUID());
              size_t offset = m_trialOffsets[ui->second];
              for (auto& dbc : m_dbcs)
              {
                if (dbc.getOperand().getUUID() == u.get().getUUID())
                {
                  dbc.assemble();
                  const auto& dofs = dbc.getDOFs();
                  axb.eliminate(dofs, offset);
                }
              }
            });

        m_assembled = true;

        return *this;
      }

      void solve(Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) override
      {
         // Assemble the system
         if (!m_assembled)
            assemble();

         // Solve the system AX = B
         solver.solve(m_stiffness, m_guess, m_mass);

         // Recover solutions
         m_us.iapply(
             [&](size_t i, auto& u)
             {
              const size_t n = u.get().getFiniteElementSpace().getSize();
              u.get().emplace().getSolution().setWeights(m_guess.segment(m_trialOffsets[i], n));
             });
      }

      Problem& operator=(const ProblemBody<OperatorType, VectorType, Real>& rhs) override
      {
        m_bft.apply([&](auto& bf) { bf.clear(); });
        m_lft.apply([&](auto& lf) { lf.clear(); });

        for (auto& bfi : rhs.getLocalBFIs())
        {
          m_bft.apply(
              [&](auto& bf)
              {
                if (bfi.getTrialFunction().getUUID() == bf.getTrialFunction().getUUID() &&
                    bfi.getTestFunction().getUUID() == bf.getTestFunction().getUUID())
                {
                  bf.add(bfi);
                }
              });
        }

        for (auto& bfi : rhs.getGlobalBFIs())
        {
          m_bft.apply(
              [&](auto& bf)
              {
                if (bfi.getTrialFunction().getUUID() == bf.getTrialFunction().getUUID() &&
                    bfi.getTestFunction().getUUID() == bf.getTestFunction().getUUID())
                {
                  bf.add(bfi);
                }
              });
        }

        for (auto& lfi : rhs.getLFIs())
        {
          m_lft.apply(
              [&](auto& lf)
              {
                if (lfi.getTestFunction().getUUID() == lf.getTestFunction().getUUID())
                {
                  lf.add(UnaryMinus(lfi));
                }
              });
        }

        m_dbcs = rhs.getDBCs();

        m_assembled = false;

        return *this;
      }

      const auto& getTrialOffsets() const
      {
        return m_trialOffsets;
      }

      const auto& getTestOffsets() const
      {
        return m_testOffsets;
      }

      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      Problem* copy() const noexcept override
      {
        assert(false);
        return nullptr;
      }

    private:
      TrialFunctionTuple m_us;
      TestFunctionTuple  m_vs;

      LinearFormTuple   m_lft;
      BilinearFormTuple m_bft;

      EssentialBoundary<ScalarType> m_dbcs;

      Boolean             m_assembled;
      VectorType          m_mass;
      VectorType          m_guess;
      OperatorType        m_stiffness;
      LinearSystemType    m_axb;

      std::array<size_t, TrialFunctionTuple::Size> m_trialOffsets;
      std::array<size_t, TestFunctionTuple::Size> m_testOffsets;

      boost::bimap<FormLanguage::Base::UUID, size_t> m_trialUUIDMap;
      boost::bimap<FormLanguage::Base::UUID, size_t> m_testUUIDMap;

      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearFormTuple>> m_bfa;
      std::unique_ptr<Assembly::AssemblyBase<VectorType, LinearFormTuple>> m_lfa;
  };

  template <class U1, class U2, class ... Us>
  Problem(U1& u1, U2& u2, Us&... us)
    -> Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Real>, Math::Vector<Real>>;
}

#include "Problem.hpp"

#endif
