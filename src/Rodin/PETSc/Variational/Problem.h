#ifndef RODIN_PETSC_VARIATIONAL_PROBLEM_H
#define RODIN_PETSC_VARIATIONAL_PROBLEM_H

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <type_traits>

#include "Rodin/Context/Local.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/PETSc/Variational/TestFunction.h"
#include "Rodin/PETSc/Variational/TrialFunction.h"
#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

#include "Rodin/PETSc/Assembly/Generic.h"
#include "Rodin/PETSc/Assembly/Sequential.h"
#include "Rodin/Variational/TestFunction.h"

namespace Rodin::Variational
{
  template <class U, class V>
  class Problem<PETSc::Math::LinearSystem, U, V>
    : public Variational::ProblemUVBase<PETSc::Math::LinearSystem, U, V>
  {
    public:
      using LinearSystemType =
        PETSc::Math::LinearSystem;

      using AssemblyType =
        PETSc::Assembly::Generic<LinearSystemType, Problem>;

      using SolverBaseType =
        Solver::SolverBase<LinearSystemType>;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using ProblemBodyType = Variational::ProblemBody<OperatorType, VectorType, ScalarType>;

      using TrialFESType =
        typename FormLanguage::Traits<U>::FESType;

      using TestFESType =
        typename FormLanguage::Traits<V>::FESType;

      using TrialFESMeshType =
        typename FormLanguage::Traits<TrialFESType>::MeshType;

      using TrialFESMeshContextType =
        typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      using TestFESMeshType =
        typename FormLanguage::Traits<TestFESType>::MeshType;

      using TestFESMeshContextType =
        typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      using Parent =
        Variational::ProblemUVBase<LinearSystemType, U, V>;

      static_assert(
          std::is_same_v<TrialFESMeshContextType, Context::Local> ||
          std::is_same_v<TrialFESMeshContextType, Context::MPI>);

      static_assert(
          std::is_same_v<TrialFESMeshContextType, TestFESMeshContextType>);

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

      constexpr
      Problem(const Problem& other)
        : Parent(other),
          m_axb(other.m_axb)
      {}

      constexpr
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_axb(std::move(other.m_axb))
      {}

      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_axb = other.m_axb;
        }
        return *this;
      }

      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_axb = std::move(other.m_axb);
        }
        return *this;
      }

      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      Problem& assemble() override
      {
        m_assembly.execute(m_axb, { m_pb, this->getTrialFunction(), this->getTestFunction() });
        m_assembled = true;
        return *this;
      }

      void solve(SolverBaseType& solver) override
      {
        auto& axb = this->getLinearSystem();
        if (!m_assembled)
           this->assemble();
        solver.solve(axb);
        this->getTrialFunction().getSolution().setData(axb.getSolution());
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
        return new Problem(*this);
      }

    private:
      Boolean m_assembled;
      ProblemBodyType m_pb;
      LinearSystemType m_axb;
      AssemblyType m_assembly;
  };

  template <class Solution, class TrialFES, class TestFES>
  Problem(
      PETSc::Variational::TrialFunction<Solution, TrialFES>&,
      PETSc::Variational::TestFunction<TestFES>&)
    -> Problem<
          PETSc::Math::LinearSystem,
          TrialFunction<Solution, TrialFES>,
          TestFunction<TestFES>>;

  // -----------------------------------------------------------------------------
  // PETSc multi-variable Problem specialization
  // - Owns: pb, trial/test tuples, UUID→block index maps, offsets, PETSc LinearSystem
  // - Delegates assembly to PETSc::Assembly::Generic<LinearSystemType, Problem>
  //   via an AssemblyInput that exposes everything Generic needs.
  // -----------------------------------------------------------------------------
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
      using LinearSystemType = PETSc::Math::LinearSystem;

      using OperatorType = typename FormLanguage::Traits<LinearSystemType>::OperatorType; // ::Mat
      using VectorType   = typename FormLanguage::Traits<LinearSystemType>::VectorType;   // ::Vec
      using ScalarType   = typename FormLanguage::Traits<LinearSystemType>::ScalarType;   // PetscScalar

      using ProblemBodyType = ProblemBody<OperatorType, VectorType, ScalarType>;
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

    public:
      // Assembly backend
      using AssemblyType =
        PETSc::Assembly::Sequential<LinearSystemType, Problem>;

      using SolverBaseType =
        Solver::SolverBase<LinearSystemType>;

      using AssemblyInput =
        Assembly::ProblemAssemblyInput<ProblemBodyType, U1, U2, U3, Us...>;

      // --------------------------
      // Ctors / assignment
      // --------------------------
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
      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      // --------------------------
      // Assembly / solve
      // --------------------------
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
      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      const auto& getTrialOffsets() const { return m_trialOffsets; }
      const auto& getTestOffsets()  const { return m_testOffsets;  }

      size_t getTotalTrialSize() const { return m_totalTrial; }
      size_t getTotalTestSize()  const { return m_totalTest;  }

      const auto& getTrialUUIDMap() const { return m_trialUUIDMap; }
      const auto& getTestUUIDMap()  const { return m_testUUIDMap;  }

      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

    private:
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
      Boolean m_assembled = false;

      TrialFunctionTuple m_us;
      TestFunctionTuple  m_vs;

      ProblemBodyType m_pb;

      std::array<size_t, TrialFunctionTuple::Size> m_trialOffsets{};
      std::array<size_t, TestFunctionTuple::Size>  m_testOffsets{};

      boost::bimap<FormLanguage::Base::UUID, size_t> m_trialUUIDMap;
      boost::bimap<FormLanguage::Base::UUID, size_t> m_testUUIDMap;

      size_t m_totalTrial = 0;
      size_t m_totalTest  = 0;

      LinearSystemType m_axb;
      AssemblyType     m_assembly;
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

  // PETSc-only CTAD guide (enabled only if ALL args are PETSc trial/test wrappers)
  template <class U1, class U2, class U3, class... Us>
    requires AllPETScTrialOrTest<U1, U2, U3, Us...>::value
  Problem(U1&, U2&, U3&, Us&...)
    -> Problem<Rodin::PETSc::Math::LinearSystem, U1, U2, U3, Us...>;
}

namespace Rodin::PETSc::Variational
{
  template <class ... Us>
  using Problem =
    Rodin::Variational::Problem<PETSc::Math::LinearSystem, Us...>;
}

#endif
