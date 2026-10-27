/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Shared PETSc scalar diffusion p and hp local/MPI verification entry point. */

#include "PETScDiffusionRefinement.h"
#include <type_traits>
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::PETScDiffusionRefinementTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeStudy(Polytope::Type geometry)
  {
    return PETScDiffusionRefinement(
      [geometry](size_t n) {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#ifdef RODIN_USE_MPI
        else
        {
          Context::MPI context(*environment, *world);
          return DistributedUniformGrid(context, geometry).makeMesh(n);
        }
#endif
      },
#ifdef RODIN_DIFFUSION_REFINEMENT_POISSON
      true
#else
      false
#endif
    );
  }

  template <class ContextType>
  class Fixture : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
      void checkPath() const
      {
#ifdef RODIN_DIFFUSION_REFINEMENT_HP
        makeStudy<ContextType>(this->GetParam()).checkCombinedPath();
#else
        makeStudy<ContextType>(this->GetParam()).checkDegreePath();
#endif
      }
      void checkThreshold() const
      {
        makeStudy<ContextType>(this->GetParam()).checkDegreeThreshold();
      }
      void checkPatch() const
      {
        makeStudy<ContextType>(this->GetParam()).template checkPatch<HighestDegree>();
      }
      void checkControl() const
      {
        makeStudy<ContextType>(this->GetParam())
          .template checkOperatorControl<HighestDegree>();
      }
      void checkSensitivity() const
      {
        makeStudy<ContextType>(this->GetParam())
          .template checkSensitivity<HighestDegree>();
      }

    private:
#ifdef RODIN_DIFFUSION_REFINEMENT_HP
      static constexpr size_t HighestDegree = 3;
#else
      static constexpr size_t HighestDegree = 4;
#endif
  };

#define RODIN_DIFFUSION_REFINEMENT_TESTS(Name)                                           \
  TEST_P(Name, ErrorsImproveAtEveryInterval)                                             \
  {                                                                                      \
    checkPath();                                                                         \
  }                                                                                      \
  TEST_P(Name, QuadraticReproductionStartsAtDegreeTwo)                                   \
  {                                                                                      \
    checkThreshold();                                                                    \
  }                                                                                      \
  TEST_P(Name, HighestDegreeReproducesQuadraticPatch)                                    \
  {                                                                                      \
    checkPatch();                                                                        \
  }                                                                                      \
  TEST_P(Name, RejectsWrongDiffusionOperator)                                            \
  {                                                                                      \
    checkControl();                                                                      \
  }                                                                                      \
  TEST_P(Name, IndependentQuadratureAndSolverSensitivity)                                \
  {                                                                                      \
    checkSensitivity();                                                                  \
  }                                                                                      \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_DIFFUSION_REFINEMENT_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_DIFFUSION_REFINEMENT_TESTS(MPITest)
#endif
#undef RODIN_DIFFUSION_REFINEMENT_TESTS
}

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::PETScDiffusionRefinementTests::environment = &env;
    Rodin::Tests::Convergence::PETScDiffusionRefinementTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
