#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

/**
 * @file Problem.h
 * @brief PETSc specialization of variational problems.
 *
 * Provides two partial specializations of @ref Rodin::Variational::Problem
 * that assemble into a @ref Rodin::PETSc::Math::LinearSystem:
 *
 * 1. **Two-field** (`Problem<LinearSystem, U, V>`): A single trial / test
 *    function pair producing a scalar linear system @f$ A\mathbf{x} = \mathbf{b} @f$.
 * 2. **Multi-field** (`Problem<LinearSystem, U1, U2, U3, Us...>`): An
 *    arbitrary number of coupled trial / test functions producing a
 *    block-structured linear system suitable for field-split
 *    preconditioning.
 *
 * Both specializations support `Context::Local` (sequential) and
 * `Context::MPI` (distributed) mesh contexts.
 *
 * @see Rodin::PETSc::Variational::TrialFunction,
 *      Rodin::PETSc::Variational::TestFunction,
 *      Rodin::PETSc::Math::LinearSystem,
 *      Rodin::Solver::KSP
 */

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <type_traits>

#include <petscmat.h>
#include <petscvec.h>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/PETSc/Variational/TrialFunction.h"
#include "Rodin/Variational/Problem.h"

#include "Rodin/Assembly/Default.h"
#include "Rodin/PETSc/Math/LinearSystem.h"

#include "Rodin/PETSc/Assembly/Sequential.h"
#include "Rodin/Variational/TestFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief PETSc variational problem for a single trial / test function pair.
   *
   * Assembles a variational formulation into a
   * @ref Rodin::PETSc::Math::LinearSystem containing the system matrix
   * @f$ A @f$, the right-hand side @f$ \mathbf{b} @f$, and the solution
   * @f$ \mathbf{x} @f$.  After calling `solve()`, the solution is
   * automatically scattered back into the trial function's grid function.
   *
   * @tparam U Trial function type (e.g.
   *           `TrialFunction<GridFunction<P1<…>>, P1<…>>`).
   * @tparam V Test function type (e.g. `TestFunction<P1<…>>`).
   *
   * @see Rodin::PETSc::Variational::Problem (convenience alias),
   *      Rodin::Solver::KSP
   */
  template <class U, class V>
  class Problem<PETSc::Math::LinearSystem, U, V>
    : public Variational::ProblemUVBase<PETSc::Math::LinearSystem, U, V>
  {
    public:
      /// @brief Type of the PETSc linear system.
      using LinearSystemType =
        PETSc::Math::LinearSystem;

      /// @brief Base solver type for this problem.
      using SolverBaseType =
        Solver::LinearSolverBase<LinearSystemType>;

      /// @brief PETSc matrix operator type.
      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      /// @brief PETSc vector type.
      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      /// @brief PETSc scalar type.
      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      /// @brief Problem body type encapsulating bilinear/linear forms.
      using ProblemBodyType = Variational::ProblemBody<OperatorType, VectorType, ScalarType>;

      /// @brief Finite element space type for the trial function.
      using TrialFESType =
        typename FormLanguage::Traits<U>::FESType;

      /// @brief Finite element space type for the test function.
      using TestFESType =
        typename FormLanguage::Traits<V>::FESType;

      /// @brief Mesh type for the trial finite element space.
      using TrialFESMeshType =
        typename FormLanguage::Traits<TrialFESType>::MeshType;

      /// @brief Context type (Local or MPI) for the trial mesh.
      using TrialFESMeshContextType =
        typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      /// @brief Mesh type for the test finite element space.
      using TestFESMeshType =
        typename FormLanguage::Traits<TestFESType>::MeshType;

      /// @brief Context type (Local or MPI) for the test mesh.
      using TestFESMeshContextType =
        typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      /// @brief Parent class type.
      using Parent =
        Variational::ProblemUVBase<LinearSystemType, U, V>;

      /// @brief Assembly type for this problem.
      using AssemblyType =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<LinearSystemType, Problem>;

      static_assert(
          std::is_same_v<TrialFESMeshContextType, TestFESMeshContextType>);

      /**
       * @brief Constructs the problem from a trial and test function pair.
       * @param u Trial function.
       * @param v Test function.
       */
      Problem(U& u, V& v)
        : Parent(u, v),
          m_axb(
              [&]() -> MPI_Comm
              {
                if constexpr (std::is_same_v<TrialFESMeshContextType, Context::Local>)
                {
                  return PETSC_COMM_SELF;
                }
                else if constexpr (std::is_same_v<TrialFESMeshContextType, Context::MPI>)
                {
                  const auto& fes = u.getFiniteElementSpace();
                  const auto& mesh = fes.getMesh();
                  const auto& ctx = mesh.getContext();
                  const MPI_Comm comm = ctx.getCommunicator();
                  return comm;
                }
                else
                {
                  assert(false);
                }
              }())
      {}

      /// @brief Copy constructor.
      constexpr
      Problem(const Problem& other)
        : Parent(other),
          m_axb(other.m_axb)
      {}

      /// @brief Move constructor.
      constexpr
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_axb(std::move(other.m_axb))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Problem to copy.
       * @return Reference to this problem.
       */
      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_axb = other.m_axb;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Problem to move from.
       * @return Reference to this problem.
       */
      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_axb = std::move(other.m_axb);
        }
        return *this;
      }

      /**
       * @brief Assigns a problem body (bilinear and linear forms).
       * @param[in] rhs Problem body to assign.
       * @return Reference to this problem.
       */
      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      /// @brief Assembles the variational formulation into the linear system.
      Problem& assemble() override
      {
        m_assembly.execute(m_axb, { m_pb, this->getTrialFunction(), this->getTestFunction() });
        m_assembled = true;
        return *this;
      }

      /**
       * @brief Solves the assembled linear system and scatters the solution.
       * @param[in] solver Linear solver to use.
       */
      void solve(SolverBaseType& solver) override
      {
        auto& axb = this->getLinearSystem();
        if (!m_assembled)
           this->assemble();
        solver.solve(axb);
        this->getTrialFunction().getSolution().setData(axb.getSolution());
      }

      /// @brief Returns a mutable reference to the linear system.
      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      /// @brief Returns a read-only reference to the linear system.
      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      /// @brief Creates a heap-allocated copy of this problem.
      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

    private:
      Boolean m_assembled;       ///< Whether the problem has been assembled.
      ProblemBodyType m_pb;      ///< The problem body (bilinear/linear forms).
      LinearSystemType m_axb;    ///< The assembled linear system.
      AssemblyType m_assembly;   ///< The assembly strategy.
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for two-field PETSc Problem.
   */
  template <class Solution, class TrialFES, class TestFES>
  Problem(
      PETSc::Variational::TrialFunction<Solution, TrialFES>&,
      PETSc::Variational::TestFunction<TestFES>&)
    -> Problem<
          PETSc::Math::LinearSystem,
          TrialFunction<Solution, TrialFES>,
          TestFunction<TestFES>>;

  /**
   * @brief PETSc variational problem for multiple coupled trial / test
   *        functions.
   *
   * Supports coupled multi-physics systems (e.g. Stokes, fluid-structure
   * interaction) by accepting an arbitrary number of PETSc trial and test
   * function arguments.  Assembles a block-structured
   * @ref Rodin::PETSc::Math::LinearSystem whose global matrix and vectors
   * are partitioned by field DOF offsets.
   *
   * After calling `solve()`, the global solution vector is sliced and
   * scattered back into each trial function's grid function using the
   * precomputed offset arrays.
   *
   * ## Block Structure
   *
   * For @f$ K @f$ trial fields with sizes @f$ N_1, \ldots, N_K @f$, the
   * global system has size @f$ N = \sum_{k=1}^K N_k @f$.  The DOF
   * offset of the @f$ k @f$-th field is @f$ O_k = \sum_{j=1}^{k-1} N_j @f$
   * (with @f$ O_1 = 0 @f$).
   *
   * ## Field-Split Preconditioning
   *
   * Call `setFieldSplits()` after `assemble()` to create PETSc `IS`
   * objects for each trial field.  These are passed to `PCFIELDSPLIT` so
   * that the preconditioner can be configured per-block from the command
   * line or programmatically.
   *
   * @tparam U1  First function type (trial or test).
   * @tparam U2  Second function type.
   * @tparam U3  Third function type.
   * @tparam Us  Additional function types.
   *
   * @see Rodin::PETSc::Variational::Problem (convenience alias),
   *      Rodin::PETSc::Math::LinearSystem::FieldSplits
   */
  template <class U1, class U2, class U3, class... Us>
  class Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>
    : public ProblemBase<PETSc::Math::LinearSystem>
  {
    private:
      template <class T>
      struct IsTrialOrTestFunction
      {
        static constexpr Boolean Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
      };

      static_assert(
        Utility::ParameterPack<U1, U2, U3, Us...>::template All<IsTrialOrTestFunction>::Value);

    public:
      /// @brief Type of the PETSc linear system.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief PETSc matrix operator type.
      using OperatorType = typename FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      /// @brief PETSc vector type.
      using VectorType   = typename FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec
      /// @brief PETSc scalar type.
      using ScalarType   = typename FormLanguage::Traits<LinearSystemType>::ScalarType;   // PetscScalar

      /// @brief Problem body type encapsulating forms.
      using ProblemBodyType = ProblemBody<OperatorType, VectorType, ScalarType>;
      /// @brief Parent class type.
      using Parent          = ProblemBase<LinearSystemType>;

    private:
      // --------------------------
      // Helpers to build tuples
      // --------------------------
      template <class T>
      struct GetFES;

      template <class T>
      struct GetFES<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::FESType;
      };

      template <class T>
      struct GetSolution;

      template <class T>
      struct GetSolution<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::SolutionType;
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

      using AllTuple =
        Tuple<
          std::reference_wrapper<U1>,
          std::reference_wrapper<U2>,
          std::reference_wrapper<U3>,
          std::reference_wrapper<Us>...>;

      using TrialFunctionTuple =
        decltype(std::declval<AllTuple>()
                 .template filter<IsTrialFunctionReferenceWrapper>());

      using TestFunctionTuple =
        decltype(std::declval<AllTuple>()
                 .template filter<IsTestFunctionReferenceWrapper>());

      using TrialFESTuple =
        typename Utility::Extract<TrialFunctionTuple>::template Type<GetFES>;

      using TestFESTuple =
        typename Utility::Extract<TestFunctionTuple>::template Type<GetFES>;

      using SolutionTuple =
        typename Utility::Extract<TrialFunctionTuple>::template Type<GetSolution>;

      using U1FESType = typename GetFES<std::reference_wrapper<U1>>::Type;

      using U2FESType = typename GetFES<std::reference_wrapper<U2>>::Type;

      using U3FESType = typename GetFES<std::reference_wrapper<U3>>::Type;

      using U1FESMeshType = typename FormLanguage::Traits<U1FESType>::MeshType;

      using U1FESMeshContextType = typename FormLanguage::Traits<U1FESMeshType>::ContextType;

      using U2FESMeshType = typename FormLanguage::Traits<U2FESType>::MeshType;

      using U2FESMeshContextType = typename FormLanguage::Traits<U2FESMeshType>::ContextType;

    public:
      /// @brief Assembly type for this problem.
      using AssemblyType =
        typename Assembly::Default<U1FESMeshContextType, U2FESMeshContextType>
          ::template Type<LinearSystemType, Problem>;

      /// @brief Base solver type for this problem.
      using SolverBaseType =
        Solver::LinearSolverBase<LinearSystemType>;

      /// @brief Input data structure for the assembly pipeline.
      using AssemblyInput =
        Assembly::ProblemAssemblyInput<ProblemBodyType, U1, U2, U3, Us...>;

      // --------------------------
      // Ctors / assignment
      // --------------------------
      /**
       * @brief Constructs a multi-field problem from trial and test functions.
       * @param u1 First function (trial or test).
       * @param u2 Second function.
       * @param u3 Third function.
       * @param us Additional functions.
       */
      Problem(U1& u1, U2& u2, U3& u3, Us&... us)
        : m_assembled(false),
          m_us(AllTuple{std::ref(u1), std::ref(u2), std::ref(u3), std::ref(us)...}
                .template filter<IsTrialFunctionReferenceWrapper>()),
          m_vs(AllTuple{std::ref(u1), std::ref(u2), std::ref(u3), std::ref(us)...}
                .template filter<IsTestFunctionReferenceWrapper>()),
          m_axb(deduceCommunicator())
      {
        buildUUIDMaps();
      }

      /// @brief Copy constructor.
      Problem(const Problem& other)
        : Parent(other),
          m_assembled(other.m_assembled),
          m_us(other.m_us),
          m_vs(other.m_vs),
          m_pb(other.m_pb),
          m_trialOffsets(other.m_trialOffsets),
          m_testOffsets(other.m_testOffsets),
          m_trialUUIDMap(other.m_trialUUIDMap),
          m_testUUIDMap(other.m_testUUIDMap),
          m_totalTrial(other.m_totalTrial),
          m_totalTest(other.m_totalTest),
          m_axb(other.m_axb),
          m_assembly(other.m_assembly)
      {}

      /// @brief Move constructor.
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_assembled(std::exchange(other.m_assembled, false)),
          m_us(std::move(other.m_us)),
          m_vs(std::move(other.m_vs)),
          m_pb(std::move(other.m_pb)),
          m_trialOffsets(std::move(other.m_trialOffsets)),
          m_testOffsets(std::move(other.m_testOffsets)),
          m_trialUUIDMap(std::move(other.m_trialUUIDMap)),
          m_testUUIDMap(std::move(other.m_testUUIDMap)),
          m_totalTrial(std::exchange(other.m_totalTrial, 0)),
          m_totalTest(std::exchange(other.m_totalTest, 0)),
          m_axb(std::move(other.m_axb)),
          m_assembly(std::move(other.m_assembly))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Problem to copy.
       * @return Reference to this problem.
       */
      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_assembled    = other.m_assembled;
          m_us           = other.m_us;
          m_vs           = other.m_vs;
          m_pb           = other.m_pb;
          m_trialOffsets = other.m_trialOffsets;
          m_testOffsets  = other.m_testOffsets;
          m_trialUUIDMap = other.m_trialUUIDMap;
          m_testUUIDMap  = other.m_testUUIDMap;
          m_totalTrial   = other.m_totalTrial;
          m_totalTest    = other.m_totalTest;
          m_axb          = other.m_axb;
          m_assembly     = other.m_assembly;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Problem to move from.
       * @return Reference to this problem.
       */
      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_assembled    = std::exchange(other.m_assembled, false);
          m_us           = std::move(other.m_us);
          m_vs           = std::move(other.m_vs);
          m_pb           = std::move(other.m_pb);
          m_trialOffsets = std::move(other.m_trialOffsets);
          m_testOffsets  = std::move(other.m_testOffsets);
          m_trialUUIDMap = std::move(other.m_trialUUIDMap);
          m_testUUIDMap  = std::move(other.m_testUUIDMap);
          m_totalTrial   = std::exchange(other.m_totalTrial, 0);
          m_totalTest    = std::exchange(other.m_totalTest, 0);
          m_axb          = std::move(other.m_axb);
          m_assembly     = std::move(other.m_assembly);
        }
        return *this;
      }

      // --------------------------
      // ProblemBody binding
      // --------------------------
      /**
       * @brief Assigns a problem body (bilinear and linear forms).
       * @param[in] rhs Problem body to assign.
       * @return Reference to this problem.
       */
      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      // --------------------------
      // Assembly / solve
      // --------------------------
      /// @brief Assembles the block-structured variational formulation.
      Problem& assemble() override
      {
        computeOffsets();

        AssemblyInput in{
          m_pb, m_us, m_vs,
          m_trialOffsets, m_testOffsets,
          m_trialUUIDMap, m_testUUIDMap,
          m_totalTrial, m_totalTest
        };

        m_assembly.execute(m_axb, in);

        m_assembled = true;
        return *this;
      }

      /**
       * @brief Solves the assembled block system and scatters sub-solutions.
       * @param[in] solver Linear solver to use.
       */
      void solve(SolverBaseType& solver) override
      {
        auto& axb = getLinearSystem();
        if (!m_assembled)
          assemble();

        solver.solve(axb);

        // Scatter global solution back into each trial function solution
        m_us.iapply(
          [&](size_t i, auto& u)
          {
            u.get().getSolution().setData(axb.getSolution(), m_trialOffsets[i]);
          });
      }

      // --------------------------
      // Accessors (useful for solvers / debugging)
      // --------------------------
      /// @brief Returns a mutable reference to the linear system.
      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      /// @brief Returns a read-only reference to the linear system.
      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      /// @brief Returns the DOF offset array for trial fields.
      const auto& getTrialOffsets() const { return m_trialOffsets; }
      /// @brief Returns the DOF offset array for test fields.
      const auto& getTestOffsets()  const { return m_testOffsets;  }

      /// @brief Returns the total number of trial DOFs across all fields.
      size_t getTotalTrialSize() const { return m_totalTrial; }
      /// @brief Returns the total number of test DOFs across all fields.
      size_t getTotalTestSize()  const { return m_totalTest;  }

      /// @brief Returns the UUID-to-index map for trial functions.
      const auto& getTrialUUIDMap() const { return m_trialUUIDMap; }
      /// @brief Returns the UUID-to-index map for test functions.
      const auto& getTestUUIDMap()  const { return m_testUUIDMap;  }

      /// @brief Creates a heap-allocated copy of this problem.
      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

      /**
       * @brief Configures PETSc field-split index sets for block preconditioning.
       *
       * Must be called after assemble(). Creates one `IS` per trial field
       * using the computed DOF offsets.
       */
      void setFieldSplits()
      {
        if (!m_assembled)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "setFieldSplits() can only be called after assemble()."
            << Alert::Raise;
        }

        assert(m_assembled);

        PetscErrorCode ierr;

        using Split = typename LinearSystemType::FieldSplits::Split;

        std::vector<Split> splits;
        splits.reserve(TrialFunctionTuple::Size);

        // Sizes in the same order as m_us (== offsets order)
        std::array<PetscInt, TrialFunctionTuple::Size> sz{};
        m_us
          .map([](const auto& u)
          {
            return static_cast<PetscInt>(u.get().getFiniteElementSpace().getSize());
          })
          .iapply([&](const Index i, const PetscInt s) { sz[i] = s; });

        // Build one IS per trial field (block layout: [u][p][lambda]...)
        for (size_t k = 0; k < TrialFunctionTuple::Size; ++k)
        {
          const PetscInt n     = sz[k];
          const PetscInt start = static_cast<PetscInt>(m_trialOffsets[k]);

          ::IS is = PETSC_NULLPTR;
          ierr = ISCreateStride(
            m_axb.getCommunicator(), // PETSC_COMM_SELF in sequential
            n,
            start,
            1,
            &is);
          assert(ierr == PETSC_SUCCESS);

          std::string name;
          // Fetch the name from the k-th trial function (same order as splits)
          m_us.iapply([&](size_t i, const auto& uref)
          {
            if (i != k)
              return;

            const auto opt = uref.get().getName(); // Optional<StringView>
            if (opt.has_value())
              name = std::string(opt->data(), opt->size());
          });

          if (name.empty())
          {
            name = std::to_string(k); // deterministic fallback
          }

          splits.push_back(Split{ std::move(name), is });
        }

        // Store (takes ownership of IS handles)
        m_axb.setFieldSplits(typename LinearSystemType::FieldSplits{ std::move(splits) });

        (void) ierr;
      }

    private:

      /// @brief Builds UUID-to-index maps for all trial and test functions.
      void buildUUIDMaps()
      {
        m_trialUUIDMap.clear();
        m_testUUIDMap.clear();

        m_us.iapply(
          [&](size_t i, const auto& u)
          {
            m_trialUUIDMap.right.insert({ i, u.get().getUUID() });
          });

        m_vs.iapply(
          [&](size_t i, const auto& v)
          {
            m_testUUIDMap.right.insert({ i, v.get().getUUID() });
          });
      }

      /// @brief Computes DOF offset arrays and totals for trial/test fields.
      void computeOffsets()
      {
        // Trial offsets + total
        {
          std::array<size_t, TrialFunctionTuple::Size> sz{};
          m_us
            .map([](const auto& u)
            {
              return static_cast<size_t>(u.get().getFiniteElementSpace().getSize());
            })
            .iapply([&](const Index i, const size_t s) { sz[i] = s; });

          m_trialOffsets[0] = 0;
          for (size_t i = 0; i + 1 < TrialFunctionTuple::Size; ++i)
            m_trialOffsets[i + 1] = m_trialOffsets[i] + sz[i];

          m_totalTrial = 0;
          for (size_t i = 0; i < TrialFunctionTuple::Size; ++i)
            m_totalTrial += sz[i];
        }

        // Test offsets + total
        {
          std::array<size_t, TestFunctionTuple::Size> sz{};
          m_vs
            .map([](const auto& v)
            {
              return static_cast<size_t>(v.get().getFiniteElementSpace().getSize());
            })
            .iapply([&](const Index i, const size_t s) { sz[i] = s; });

          m_testOffsets[0] = 0;
          for (size_t i = 0; i + 1 < TestFunctionTuple::Size; ++i)
            m_testOffsets[i + 1] = m_testOffsets[i] + sz[i];

          m_totalTest = 0;
          for (size_t i = 0; i < TestFunctionTuple::Size; ++i)
            m_totalTest += sz[i];
        }
      }

      /// @brief Deduces the MPI communicator from the trial mesh contexts.
      MPI_Comm deduceCommunicator() const
      {
        // Take mesh context from the first trial function in the tuple.
        MPI_Comm comm = PETSC_COMM_SELF;

        m_us.apply(
          [&](const auto& uref)
          {
            if (comm != PETSC_COMM_SELF)
              return;

            const auto& fes  = uref.get().getFiniteElementSpace();
            const auto& mesh = fes.getMesh();

            using MeshType = std::decay_t<decltype(mesh)>;
            using Ctx      = typename FormLanguage::Traits<MeshType>::ContextType;

            if constexpr (std::is_same_v<Ctx, Context::Local>)
            {
              comm = PETSC_COMM_SELF;
            }
            else if constexpr (std::is_same_v<Ctx, Context::MPI>)
            {
              comm = mesh.getContext().getCommunicator();
            }
            else
            {
              static_assert(!sizeof(Ctx), "Unsupported mesh context for PETSc Problem.");
            }
          });

        return comm;
      }

    private:
      Boolean m_assembled = false;  ///< Whether the problem has been assembled.

      TrialFunctionTuple m_us;  ///< Tuple of trial function references.
      TestFunctionTuple  m_vs;  ///< Tuple of test function references.

      ProblemBodyType m_pb;  ///< The problem body (bilinear/linear forms).

      std::array<size_t, TrialFunctionTuple::Size> m_trialOffsets{};  ///< DOF offsets for trial fields.
      std::array<size_t, TestFunctionTuple::Size>  m_testOffsets{};   ///< DOF offsets for test fields.

      boost::bimap<FormLanguage::Base::UUID, size_t> m_trialUUIDMap;  ///< Trial UUID-to-index map.
      boost::bimap<FormLanguage::Base::UUID, size_t> m_testUUIDMap;   ///< Test UUID-to-index map.

      size_t m_totalTrial = 0;  ///< Total trial DOF count.
      size_t m_totalTest  = 0;  ///< Total test DOF count.

      LinearSystemType m_axb;    ///< The assembled linear system.
      AssemblyType     m_assembly;  ///< The assembly strategy.
  };

  template <class T>
  struct IsPETScTrialFunction : std::false_type {};

  template <class Sol, class FES>
  struct IsPETScTrialFunction<Rodin::PETSc::Variational::TrialFunction<Sol, FES>>
    : std::true_type {};

  template <class T>
  struct IsPETScTestFunction : std::false_type {};

  template <class FES>
  struct IsPETScTestFunction<Rodin::PETSc::Variational::TestFunction<FES>>
    : std::true_type {};

  template <class... Ts>
  struct AllPETScTrialOrTest;

  template <>
  struct AllPETScTrialOrTest<> : std::true_type {};

  template <class T, class... Ts>
  struct AllPETScTrialOrTest<T, Ts...>
    : std::bool_constant<
        (IsPETScTrialFunction<std::decay_t<T>>::value ||
         IsPETScTestFunction<std::decay_t<T>>::value) &&
        AllPETScTrialOrTest<Ts...>::value> {};

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for multi-field PETSc Problem.
   */
  // PETSc-only CTAD guide (enabled only if ALL args are PETSc trial/test wrappers)
  template <class U1, class U2, class U3, class... Us>
    requires AllPETScTrialOrTest<U1, U2, U3, Us...>::value
  Problem(U1&, U2&, U3&, Us&...)
    -> Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>;
}

namespace Rodin::PETSc::Variational
{
  /**
   * @brief Convenient PETSc alias for Rodin::Variational::Problem.
   */
  template <class ... Us>
  using Problem =
    Rodin::Variational::Problem<PETSc::Math::LinearSystem, Us...>;
}

#endif
