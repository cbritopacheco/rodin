/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <mpi.h>
#include <petsc.h>

#include <boost/bimap.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <template <class, class> class Assembler>
  void checkPETScVectorMasterIdentification()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 slaveFES(mesh);
    P1 masterFES(mesh, mesh.getSpaceDimension());

    PETSc::Variational::TrialFunction u(slaveFES);
    PETSc::Variational::TestFunction  v(slaveFES);
    PETSc::Variational::TrialFunction eta(masterFES);
    PETSc::Variational::TestFunction  zeta(masterFES);

    auto bc =
      DirichletBC(
          u,
          RealFunction(2.0) * eta.x() + RealFunction(-0.5) * eta.y());
    bc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(bc.getDOFs()));
    const auto& ident = std::get<IdentifiedDOFs>(bc.getDOFs());
    ASSERT_FALSE(ident.empty());

    bool sawMultiMaster = false;
    for (const auto& [slave, row] : ident)
    {
      (void) slave;
      if (row.first.size() >= 2)
        sawMultiMaster = true;
    }
    ASSERT_TRUE(sawMultiMaster);

    const PetscInt nSlave  = static_cast<PetscInt>(slaveFES.getSize());
    const PetscInt nMaster = static_cast<PetscInt>(masterFES.getSize());
    const PetscInt nTotal  = nSlave + nMaster;

    struct MatrixEntry
    {
      PetscInt row;
      PetscInt col;
      PetscScalar value;
    };
    const std::vector<MatrixEntry> matrixEntries{
      { 0, 0, 2.0 },
      { 0, 1, 3.0 },
      { 1, 0, 5.0 },
      { 1, 1, 7.0 }
    };
    const std::vector<std::pair<PetscInt, PetscScalar>> vectorEntries{
      { 0, 11.0 },
      { 1, -13.0 }
    };

    BilinearForm uu(u, v);
    auto& op = uu.getOperator();
    PetscErrorCode ierr = MatSetSizes(op, nSlave, nSlave, nSlave, nSlave);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    for (const auto& e : matrixEntries)
    {
      ierr = MatSetValue(op, e.row, e.col, e.value, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    LinearForm loadU(v);
    auto& load = loadU.getVector();
    ierr = VecSetSizes(load, nSlave, nSlave);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    for (const auto& [row, value] : vectorEntries)
    {
      ierr = VecSetValue(load, row, value, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = VecAssemblyBegin(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    auto body = uu + bc - loadU;

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType =
      Problem<LinearSystemType,
              decltype(u), decltype(v), decltype(eta), decltype(zeta)>;

    auto trialFunctions = Tuple{ std::ref(u), std::ref(eta) };
    auto testFunctions  = Tuple{ std::ref(v), std::ref(zeta) };

    std::array<size_t, 2> trialOffsets{
      0, static_cast<size_t>(nSlave)
    };
    std::array<size_t, 2> testOffsets{
      0, static_cast<size_t>(nSlave)
    };

    boost::bimap<FormLanguage::Base::UUID, size_t> trialUUIDMap;
    boost::bimap<FormLanguage::Base::UUID, size_t> testUUIDMap;
    trialUUIDMap.right.insert({ 0, u.getUUID() });
    trialUUIDMap.right.insert({ 1, eta.getUUID() });
    testUUIDMap.right.insert({ 0, v.getUUID() });
    testUUIDMap.right.insert({ 1, zeta.getUUID() });

    Assembly::ProblemAssemblyInput<
      std::decay_t<decltype(body)>,
      decltype(u), decltype(v), decltype(eta), decltype(zeta)> input(
          body,
          trialFunctions,
          testFunctions,
          trialOffsets,
          testOffsets,
          trialUUIDMap,
          testUUIDMap,
          static_cast<size_t>(nTotal),
          static_cast<size_t>(nTotal));

    LinearSystemType ls(PETSC_COMM_SELF);
    Assembler<LinearSystemType, ProblemType> assembler;
    assembler.execute(ls, input);

    Eigen::MatrixXd expectedA =
      Eigen::MatrixXd::Zero(nTotal, nTotal);
    Eigen::VectorXd expectedB =
      Eigen::VectorXd::Zero(nTotal);

    auto expand = [&](PetscInt local)
    {
      std::vector<std::pair<PetscInt, PetscScalar>> res;
      const auto it = ident.find(static_cast<Index>(local));
      if (it == ident.end())
      {
        res.emplace_back(local, 1.0);
        return res;
      }
      const auto& masters = it->second.first;
      const auto& coeffs  = it->second.second;
      for (Index k = 0; k < static_cast<Index>(masters.size()); k++)
        res.emplace_back(nSlave + static_cast<PetscInt>(masters[k]), coeffs[k]);
      return res;
    };

    for (const auto& e : matrixEntries)
      for (const auto& r : expand(e.row))
        for (const auto& c : expand(e.col))
          expectedA(r.first, c.first) += r.second * e.value * c.second;

    for (const auto& [row, value] : vectorEntries)
      for (const auto& r : expand(row))
        expectedB(r.first) += r.second * value;

    for (const auto& [slave, row] : ident)
    {
      expectedA.row(static_cast<Eigen::Index>(slave)).setZero();
      expectedA(slave, slave) = 1.0;
      const auto& masters = row.first;
      const auto& coeffs  = row.second;
      for (Index k = 0; k < static_cast<Index>(masters.size()); k++)
        expectedA(slave, nSlave + static_cast<PetscInt>(masters[k])) -= coeffs[k];
      expectedB(slave) = 0.0;
    }

    auto& A = ls.getOperator();
    auto& b = ls.getVector();
    for (PetscInt i = 0; i < nTotal; i++)
    {
      PetscScalar bValue = 0;
      ierr = VecGetValues(b, 1, &i, &bValue);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(bValue, expectedB(i), 1e-14) << "row " << i;

      for (PetscInt j = 0; j < nTotal; j++)
      {
        PetscScalar aValue = 0;
        ierr = MatGetValues(A, 1, &i, 1, &j, &aValue);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_NEAR(aValue, expectedA(i, j), 1e-14)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  TEST(PETSc_Form, SequentialLinearFormUsesSelfCommunicatorAndAssembles)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
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
    EXPECT_EQ(commSize, 1);

    PetscInt localSize = 0;
    PetscInt globalSize = 0;
    ierr = VecGetLocalSize(b, &localSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(b, &globalSize);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localSize, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalSize, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_Form, SequentialBilinearFormUsesSelfCommunicatorAndAssembles)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
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
    EXPECT_EQ(commSize, 1);

    PetscInt localRows = 0;
    PetscInt localCols = 0;
    PetscInt globalRows = 0;
    PetscInt globalCols = 0;
    ierr = MatGetLocalSize(A, &localRows, &localCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(A, &globalRows, &globalCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    EXPECT_EQ(localRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(localCols, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalRows, static_cast<PetscInt>(fes.getSize()));
    EXPECT_EQ(globalCols, static_cast<PetscInt>(fes.getSize()));
  }

  TEST(PETSc_Form, SequentialIdentificationProjectsMatrixAndVector)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    PETSc::Variational::TrialFunction eta(fes);
    PETSc::Variational::TestFunction  zeta(fes);

    constexpr PetscScalar gamma = 2.0;
    const PetscInt n = static_cast<PetscInt>(fes.getSize());

    BilinearForm uu(u, v);
    auto& op = uu.getOperator();
    PetscErrorCode ierr = MatSetSizes(op, n, n, n, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    PetscInt rows[3] = { 0, 0, 1 };
    PetscInt cols[3] = { 0, 1, 0 };
    PetscScalar vals[3] = { 2.0, 3.0, 5.0 };
    for (PetscInt i = 0; i < 3; i++)
    {
      ierr = MatSetValue(op, rows[i], cols[i], vals[i], INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    LinearForm loadU(v);
    auto& load = loadU.getVector();
    ierr = VecSetSizes(load, n, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    PetscInt rhsRows[2] = { 0, 1 };
    PetscScalar rhsVals[2] = { 7.0, 11.0 };
    ierr = VecSetValues(load, 2, rhsRows, rhsVals, INSERT_VALUES);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyBegin(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(load);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    Problem problem(u, v, eta, zeta);
    problem = uu + DirichletBC(u, RealFunction(gamma) * eta) - loadU;
    problem.assemble();

    auto& A = problem.getLinearSystem().getOperator();
    auto& b = problem.getLinearSystem().getVector();

    PetscInt r, c;
    PetscScalar value;

    r = n + 0; c = n + 0;
    ierr = MatGetValues(A, 1, &r, 1, &c, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, gamma * 2.0 * gamma, 1e-14);

    r = n + 0; c = n + 1;
    ierr = MatGetValues(A, 1, &r, 1, &c, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, gamma * 3.0 * gamma, 1e-14);

    r = n + 1; c = n + 0;
    ierr = MatGetValues(A, 1, &r, 1, &c, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, gamma * 5.0 * gamma, 1e-14);

    r = 0; c = 0;
    ierr = MatGetValues(A, 1, &r, 1, &c, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, 1.0, 1e-14);

    r = 0; c = n + 0;
    ierr = MatGetValues(A, 1, &r, 1, &c, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, -gamma, 1e-14);

    r = 0;
    ierr = VecGetValues(b, 1, &r, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, 0.0, 1e-14);

    r = n + 0;
    ierr = VecGetValues(b, 1, &r, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, gamma * 7.0, 1e-14);

    r = n + 1;
    ierr = VecGetValues(b, 1, &r, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, gamma * 11.0, 1e-14);
  }

  TEST(PETSc_Form, SequentialVectorMasterIdentificationProjectsMatrixAndVector)
  {
    checkPETScVectorMasterIdentification<PETSc::Assembly::Sequential>();
  }

  TEST(PETSc_Form, OpenMPVectorMasterIdentificationProjectsMatrixAndVector)
  {
    checkPETScVectorMasterIdentification<PETSc::Assembly::OpenMP>();
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
