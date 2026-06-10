/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <functional>

#include <gtest/gtest.h>

#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static boost::mpi::environment*  g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  Mesh<Context::Local> makeShardableMesh()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 5, 5 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    return mesh;
  }

  Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);

    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh();
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }

    return sharder.gather(0);
  }

  TEST(PETSc_MPI_Form, DistributedLinearFormUsesMeshCommunicatorAndAssembles)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    PETSc::Variational::TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    auto& b = lf.getVector();

    MPI_Comm comm = MPI_COMM_NULL;
    PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(b), &comm);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    int commSize = 0;
    MPI_Comm_size(comm, &commSize);
    EXPECT_EQ(commSize, world.size());

    size_t begin = 0;
    size_t end = 0;
    fes.getOwnershipRange(begin, end);

    PetscInt localSize = 0;
    PetscInt globalSize = 0;
    ierr = VecGetLocalSize(b, &localSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(b, &globalSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localSize, static_cast<PetscInt>(end - begin));
    EXPECT_EQ(globalSize, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_MPI_Form, DistributedBilinearFormUsesMeshCommunicatorAndAssembles)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    auto& A = bf.getOperator();

    MPI_Comm comm = MPI_COMM_NULL;
    PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(A), &comm);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    int commSize = 0;
    MPI_Comm_size(comm, &commSize);
    EXPECT_EQ(commSize, world.size());

    size_t begin = 0;
    size_t end = 0;
    fes.getOwnershipRange(begin, end);

    PetscInt localRows = 0;
    PetscInt localCols = 0;
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    ierr = MatGetLocalSize(A, &localRows, &localCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(A, &globalRows, &globalCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    const auto ownedSize = static_cast<PetscInt>(end - begin);
    EXPECT_EQ(localRows, ownedSize);
    EXPECT_EQ(localCols, ownedSize);
    EXPECT_EQ(globalRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalCols, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_MPI_Form, DistributedIdentificationProjectsOwnedSlaveRow)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    PETSc::Variational::TrialFunction eta(fes);
    PETSc::Variational::TestFunction  zeta(fes);

    constexpr PetscScalar gamma = 2.0;
    constexpr PetscScalar defect = 3.0;
    const PetscInt n = static_cast<PetscInt>(fes.getSize());
    size_t begin = 0;
    size_t end = 0;
    fes.getOwnershipRange(begin, end);
    const PetscInt localSize = static_cast<PetscInt>(end - begin);

    auto dbc =
      DirichletBC(u, RealFunction(gamma) * eta, RealFunction(defect));
    dbc.assemble();
    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());

    bool hasOwnedConstraint = false;
    PetscInt slave = -1;
    PetscInt master = -1;
    PetscScalar coefficient = 0.0;
    for (const auto& [s, pair] : dofs)
    {
      if (begin <= s && s < end && pair.first.size() == 1)
      {
        hasOwnedConstraint = true;
        slave = static_cast<PetscInt>(s);
        master = n + static_cast<PetscInt>(pair.first.coeff(0));
        coefficient = static_cast<PetscScalar>(pair.second.coeff(0));
        break;
      }
    }

    BilinearForm uu(u, v);
    auto& op = uu.getOperator();
    PetscErrorCode ierr = MatSetSizes(op, localSize, localSize, n, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    if (hasOwnedConstraint)
    {
      ierr = MatSetValue(op, slave, slave, 2.0, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    LinearForm loadU(v);
    auto& load = loadU.getVector();
    ierr = VecSetSizes(load, localSize, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    if (hasOwnedConstraint)
    {
      const PetscScalar val = 7.0;
      ierr = VecSetValue(load, slave, val, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = VecAssemblyBegin(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    Problem problem(u, v, eta, zeta);
    problem = uu + dbc - loadU;
    problem.assemble();

    if (hasOwnedConstraint)
    {
      auto& A = problem.getLinearSystem().getOperator();
      auto& b = problem.getLinearSystem().getVector();
      PetscInt r;
      PetscInt c;
      PetscScalar value;

      r = master; c = master;
      ierr = MatGetValues(A, 1, &r, 1, &c, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, coefficient * 2.0 * coefficient, 1e-14);

      r = slave; c = slave;
      ierr = MatGetValues(A, 1, &r, 1, &c, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, 1.0, 1e-14);

      r = slave; c = master;
      ierr = MatGetValues(A, 1, &r, 1, &c, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, -coefficient, 1e-14);

      r = master;
      ierr = VecGetValues(b, 1, &r, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, coefficient * 7.0 - coefficient * 2.0 * defect, 1e-14);

      r = slave;
      ierr = VecGetValues(b, 1, &r, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, defect, 1e-14);
    }
  }

  TEST(PETSc_MPI_Form, DistributedSelfIdentificationMatchesZeroValueConstraint)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1 refFES(mesh);
    PETSc::Variational::TrialFunction uRef(refFES);
    PETSc::Variational::TestFunction  vRef(refFES);

    Problem refProblem(uRef, vRef);
    refProblem = Integral(Grad(uRef), Grad(vRef))
               - Integral(RealFunction(1.0), vRef)
               + DirichletBC(uRef, Zero());
    refProblem.assemble();

    P1 idFES(mesh);
    PETSc::Variational::TrialFunction uId(idFES);
    PETSc::Variational::TestFunction  vId(idFES);

    Problem idProblem(uId, vId);
    idProblem = Integral(Grad(uId), Grad(vId))
              - Integral(RealFunction(1.0), vId)
              + DirichletBC(uId, -uId);
    idProblem.assemble();

    auto& ARef = refProblem.getLinearSystem().getOperator();
    auto& AId  = idProblem.getLinearSystem().getOperator();
    auto& bRef = refProblem.getLinearSystem().getVector();
    auto& bId  = idProblem.getLinearSystem().getVector();

    size_t begin = 0;
    size_t end   = 0;
    refFES.getOwnershipRange(begin, end);

    const PetscInt n = static_cast<PetscInt>(refFES.getSize());
    for (PetscInt i = static_cast<PetscInt>(begin);
         i < static_cast<PetscInt>(end); i++)
    {
      PetscErrorCode ierr;
      PetscScalar refValue = 0;
      PetscScalar idValue  = 0;

      ierr = VecGetValues(bRef, 1, &i, &refValue);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecGetValues(bId, 1, &i, &idValue);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(refValue, idValue, 1e-12) << "row " << i;

      for (PetscInt j = 0; j < n; j++)
      {
        ierr = MatGetValues(ARef, 1, &i, 1, &j, &refValue);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        ierr = MatGetValues(AId, 1, &i, 1, &j, &idValue);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_NEAR(refValue, idValue, 1e-12)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  TEST(PETSc_MPI_Form, DistributedVectorMasterIdentificationProjectsOwnedSlaveRow)
  {
    const auto& world = *g_world;
    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx);

    P1 slaveFES(mesh);
    P1 masterFES(mesh, mesh.getSpaceDimension());

    PETSc::Variational::TrialFunction u(slaveFES);
    PETSc::Variational::TestFunction  v(slaveFES);
    PETSc::Variational::TrialFunction eta(masterFES);
    PETSc::Variational::TestFunction  zeta(masterFES);

    const PetscInt nSlave  = static_cast<PetscInt>(slaveFES.getSize());
    const PetscInt nMaster = static_cast<PetscInt>(masterFES.getSize());

    size_t slaveBegin = 0;
    size_t slaveEnd = 0;
    slaveFES.getOwnershipRange(slaveBegin, slaveEnd);
    const PetscInt localSlaveSize =
      static_cast<PetscInt>(slaveEnd - slaveBegin);

    size_t masterBegin = 0;
    size_t masterEnd = 0;
    masterFES.getOwnershipRange(masterBegin, masterEnd);

    auto dbc =
      DirichletBC(
          u,
          RealFunction(2.0) * eta.x() + RealFunction(-0.5) * eta.y());
    dbc.assemble();
    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(dbc.getDOFs()));
    const auto& dofs = std::get<IdentifiedDOFs>(dbc.getDOFs());

    bool hasOwnedConstraint = false;
    PetscInt slave = -1;
    std::vector<PetscInt> masters;
    std::vector<PetscScalar> coefficients;
    for (const auto& [s, pair] : dofs)
    {
      if (slaveBegin <= s && s < slaveEnd && pair.first.size() >= 2)
      {
        hasOwnedConstraint = true;
        slave = static_cast<PetscInt>(s);
        for (Index k = 0; k < static_cast<Index>(pair.first.size()); k++)
        {
          masters.push_back(nSlave + static_cast<PetscInt>(pair.first[k]));
          coefficients.push_back(static_cast<PetscScalar>(pair.second[k]));
        }
        break;
      }
    }

    BilinearForm uu(u, v);
    auto& op = uu.getOperator();
    PetscErrorCode ierr =
      MatSetSizes(op, localSlaveSize, localSlaveSize, nSlave, nSlave);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    if (hasOwnedConstraint)
    {
      ierr = MatSetValue(op, slave, slave, 2.0, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    LinearForm loadU(v);
    auto& load = loadU.getVector();
    ierr = VecSetSizes(load, localSlaveSize, nSlave);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    if (hasOwnedConstraint)
    {
      const PetscScalar val = 7.0;
      ierr = VecSetValue(load, slave, val, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = VecAssemblyBegin(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    Problem problem(u, v, eta, zeta);
    problem = uu + dbc - loadU;
    problem.assemble();

    if (hasOwnedConstraint)
    {
      ASSERT_GE(masters.size(), static_cast<size_t>(2));

      auto& A = problem.getLinearSystem().getOperator();
      auto& b = problem.getLinearSystem().getVector();

      PetscScalar value = 0;
      for (size_t a = 0; a < masters.size(); a++)
      {
        PetscInt r = masters[a];
        ierr = VecGetValues(b, 1, &r, &value);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_NEAR(value, coefficients[a] * 7.0, 1e-14)
          << "projected vector row " << r;

        for (size_t cidx = 0; cidx < masters.size(); cidx++)
        {
          r = masters[a];
          PetscInt c = masters[cidx];
          ierr = MatGetValues(A, 1, &r, 1, &c, &value);
          ASSERT_EQ(ierr, PETSC_SUCCESS);
          EXPECT_NEAR(
              value,
              coefficients[a] * 2.0 * coefficients[cidx],
              1e-14)
            << "projected matrix entry (" << r << ", " << c << ")";
        }
      }

      PetscInt r = slave;
      PetscInt c = slave;
      ierr = MatGetValues(A, 1, &r, 1, &c, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, 1.0, 1e-14);

      for (size_t a = 0; a < masters.size(); a++)
      {
        r = slave;
        c = masters[a];
        ierr = MatGetValues(A, 1, &r, 1, &c, &value);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_NEAR(value, -coefficients[a], 1e-14)
          << "reconstruction entry (" << r << ", " << c << ")";
      }

      r = slave;
      ierr = VecGetValues(b, 1, &r, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, 0.0, 1e-14);
    }

    (void) nMaster;
    (void) masterBegin;
    (void) masterEnd;
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
