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
  class ToggleLocalBilinearIntegrator final
    : public LocalBilinearFormIntegratorBase<PetscScalar>
  {
    public:
      using Parent = LocalBilinearFormIntegratorBase<PetscScalar>;

      template <class TrialFunctionType, class TestFunctionType>
      ToggleLocalBilinearIntegrator(
          const TrialFunctionType& u,
          const TestFunctionType& v,
          const PetscScalar* offdiag)
        : Parent(u, v),
          m_offdiag(offdiag)
      {}

      ToggleLocalBilinearIntegrator(const ToggleLocalBilinearIntegrator& other)
        : Parent(other),
          m_polytope(other.m_polytope),
          m_offdiag(other.m_offdiag)
      {}

      const Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      ToggleLocalBilinearIntegrator& setPolytope(const Polytope& polytope) override
      {
        m_polytope = &polytope;
        return *this;
      }

      PetscScalar integrate(size_t tr, size_t te) override
      {
        return tr == te ? PetscScalar(1) : *m_offdiag;
      }

      Region getRegion() const override
      {
        return Region::Cells;
      }

      ToggleLocalBilinearIntegrator* copy() const noexcept override
      {
        return new ToggleLocalBilinearIntegrator(*this);
      }

    private:
      const Polytope* m_polytope = nullptr;
      const PetscScalar* m_offdiag = nullptr;
  };

  template <template <class, class> class Assembler>
  void checkPETScReassemblyKeepsExplicitZeroStructuralEntries()
  {
    auto mesh = Mesh<Context::Local>::Builder()
      .initialize(1)
      .nodes(2)
      .vertex({0.0})
      .vertex({1.0})
      .polytope(Polytope::Type::Segment, {0, 1})
      .finalize();

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    PetscScalar offdiag = 0;
    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    using ProblemBodyType = typename ProblemType::ProblemBodyType;

    ToggleLocalBilinearIntegrator integrator(u, v, &offdiag);
    ProblemBodyType body(integrator);

    Assembly::ProblemAssemblyInput<
      ProblemBodyType,
      decltype(u), decltype(v)> input(body, u, v);

    LinearSystemType ls(PETSC_COMM_SELF);
    Assembler<LinearSystemType, ProblemType> assembler;
    assembler.execute(ls, input);

    auto& A = ls.getOperator();
    PetscErrorCode ierr =
      MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    offdiag = 2;
    assembler.execute(ls, input);

    PetscInt row = 0;
    PetscInt col = 1;
    PetscScalar value = 0;
    ierr = MatGetValues(A, 1, &row, 1, &col, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, offdiag, 1e-14);

    row = 1;
    col = 0;
    value = 0;
    ierr = MatGetValues(A, 1, &row, 1, &col, &value);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    EXPECT_NEAR(value, offdiag, 1e-14);
  }

  template <template <class, class> class Assembler>
  void checkPETScAffineIdentificationDefect()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);
    PETSc::Variational::TrialFunction eta(fes);
    PETSc::Variational::TestFunction  zeta(fes);

    constexpr PetscScalar gamma = 2.0;
    constexpr PetscScalar defect = 3.0;
    const PetscInt n = static_cast<PetscInt>(fes.getSize());
    const PetscInt nTotal = 2 * n;

    auto bc =
      DirichletBC(u, RealFunction(gamma) * eta, RealFunction(defect));

    BilinearForm uu(u, v);
    auto& op = uu.getOperator();
    PetscErrorCode ierr = MatSetSizes(op, n, n, n, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetFromOptions(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatSetUp(op);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    for (PetscInt i = 0; i < n; i++)
    {
      ierr = MatSetValue(op, i, i, 1.0, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    }
    ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    LinearForm zero(v);
    auto& zeroVec = zero.getVector();
    ierr = VecSetSizes(zeroVec, n, n);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetFromOptions(zeroVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecSetUp(zeroVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecZeroEntries(zeroVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyBegin(zeroVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecAssemblyEnd(zeroVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    auto body = uu + bc - zero;

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType =
      Problem<LinearSystemType,
              decltype(u), decltype(v), decltype(eta), decltype(zeta)>;

    auto trialFunctions = Tuple{ std::ref(u), std::ref(eta) };
    auto testFunctions  = Tuple{ std::ref(v), std::ref(zeta) };

    std::array<size_t, 2> offsets{ 0, static_cast<size_t>(n) };

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
          offsets,
          offsets,
          trialUUIDMap,
          testUUIDMap,
          static_cast<size_t>(nTotal),
          static_cast<size_t>(nTotal));

    LinearSystemType ls(PETSC_COMM_SELF);
    Assembler<LinearSystemType, ProblemType> assembler;
    assembler.execute(ls, input);

    auto& A = ls.getOperator();
    auto& b = ls.getVector();

    for (PetscInt i = 0; i < n; i++)
    {
      PetscScalar value = 0;
      ierr = VecGetValues(b, 1, &i, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, defect, 1e-14) << "row " << i;

      PetscInt col = i;
      ierr = MatGetValues(A, 1, &i, 1, &col, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, 1.0, 1e-14) << "entry (" << i << ", " << col << ")";

      col = n + i;
      ierr = MatGetValues(A, 1, &i, 1, &col, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, -gamma, 1e-14) << "entry (" << i << ", " << col << ")";

      PetscInt row = n + i;
      ierr = VecGetValues(b, 1, &row, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, -gamma * defect, 1e-14) << "projected row " << row;

      col = n + i;
      ierr = MatGetValues(A, 1, &row, 1, &col, &value);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(value, gamma * gamma, 1e-14)
        << "projected entry (" << row << ", " << col << ")";
    }
  }

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

  template <template <class, class> class Assembler>
  void checkPETScSelfIdentificationMatchesZeroValueConstraint()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 refFES(mesh);
    PETSc::Variational::TrialFunction uRef(refFES);
    PETSc::Variational::TestFunction  vRef(refFES);

    const PetscInt n = static_cast<PetscInt>(refFES.getSize());

    auto setOperator = [&](auto& form)
    {
      auto& op = form.getOperator();
      PetscErrorCode ierr = MatSetSizes(op, n, n, n, n);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = MatSetFromOptions(op);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = MatSetUp(op);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      const PetscInt rows[4] = { 0, 0, 1, 1 };
      const PetscInt cols[4] = { 0, 1, 0, 1 };
      const PetscScalar vals[4] = { 2.0, 3.0, 5.0, 7.0 };
      for (PetscInt k = 0; k < 4; k++)
      {
        ierr = MatSetValue(op, rows[k], cols[k], vals[k], INSERT_VALUES);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
      }
      ierr = MatAssemblyBegin(op, MAT_FINAL_ASSEMBLY);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = MatAssemblyEnd(op, MAT_FINAL_ASSEMBLY);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    };

    auto setVector = [&](auto& form)
    {
      auto& vec = form.getVector();
      PetscErrorCode ierr = VecSetSizes(vec, n, n);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecSetFromOptions(vec);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      const PetscInt rows[2] = { 0, 1 };
      const PetscScalar vals[2] = { 11.0, -13.0 };
      ierr = VecSetValues(vec, 2, rows, vals, INSERT_VALUES);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecAssemblyBegin(vec);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecAssemblyEnd(vec);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
    };

    BilinearForm uuRef(uRef, vRef);
    setOperator(uuRef);
    LinearForm loadRef(vRef);
    setVector(loadRef);

    auto refBody = uuRef + DirichletBC(uRef, Zero()) - loadRef;

    P1 idFES(mesh);
    PETSc::Variational::TrialFunction uId(idFES);
    PETSc::Variational::TestFunction  vId(idFES);

    BilinearForm uuId(uId, vId);
    setOperator(uuId);
    LinearForm loadId(vId);
    setVector(loadId);

    auto idBody = uuId + DirichletBC(uId, -uId) - loadId;

    using LinearSystemType = PETSc::Math::LinearSystem;
    using ProblemType =
      Problem<LinearSystemType, decltype(uRef), decltype(vRef)>;

    LinearSystemType refLS(PETSC_COMM_SELF);
    LinearSystemType idLS(PETSC_COMM_SELF);
    Assembly::ProblemAssemblyInput<
      std::decay_t<decltype(refBody)>,
      decltype(uRef), decltype(vRef)> refInput(refBody, uRef, vRef);
    Assembly::ProblemAssemblyInput<
      std::decay_t<decltype(idBody)>,
      decltype(uId), decltype(vId)> idInput(idBody, uId, vId);

    Assembler<LinearSystemType, ProblemType> assembler;
    assembler.execute(refLS, refInput);
    assembler.execute(idLS, idInput);

    auto& ARef = refLS.getOperator();
    auto& AId  = idLS.getOperator();
    auto& bRef = refLS.getVector();
    auto& bId  = idLS.getVector();

    for (PetscInt i = 0; i < n; i++)
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

  void expectSameStandaloneVector(Vec expected, Vec actual)
  {
    PetscErrorCode ierr;
    PetscInt nExpected = 0;
    PetscInt nActual = 0;
    ierr = VecGetSize(expected, &nExpected);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecGetSize(actual, &nActual);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(nExpected, nActual);

    for (PetscInt i = 0; i < nExpected; i++)
    {
      PetscScalar a = 0;
      PetscScalar b = 0;
      ierr = VecGetValues(expected, 1, &i, &a);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      ierr = VecGetValues(actual, 1, &i, &b);
      ASSERT_EQ(ierr, PETSC_SUCCESS);
      EXPECT_NEAR(PetscRealPart(a - b), 0.0, 1e-14) << "row " << i;
    }
  }

  void expectSameStandaloneMatrix(Mat expected, Mat actual)
  {
    PetscErrorCode ierr;
    PetscInt expectedRows = 0;
    PetscInt expectedCols = 0;
    PetscInt actualRows = 0;
    PetscInt actualCols = 0;
    ierr = MatGetSize(expected, &expectedRows, &expectedCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatGetSize(actual, &actualRows, &actualCols);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ASSERT_EQ(expectedRows, actualRows);
    ASSERT_EQ(expectedCols, actualCols);

    for (PetscInt i = 0; i < expectedRows; i++)
    {
      for (PetscInt j = 0; j < expectedCols; j++)
      {
        PetscScalar a = 0;
        PetscScalar b = 0;
        ierr = MatGetValues(expected, 1, &i, 1, &j, &a);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        ierr = MatGetValues(actual, 1, &i, 1, &j, &b);
        ASSERT_EQ(ierr, PETSC_SUCCESS);
        EXPECT_NEAR(PetscRealPart(a - b), 0.0, 1e-14)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

#ifdef RODIN_USE_OPENMP
  void checkPETScStandaloneOpenMPFormsMatchSequential()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction  v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(2.0), v);

    using LFType = decltype(lf);
    Vec seqVec = nullptr;
    Vec ompVec = nullptr;
    PetscErrorCode ierr = VecCreate(PETSC_COMM_SELF, &seqVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecCreate(PETSC_COMM_SELF, &ompVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PETSc::Assembly::Sequential<::Vec, LFType> seqLF;
    PETSc::Assembly::OpenMP<::Vec, LFType> ompLF;
    seqLF.execute(seqVec, { fes, lf.getIntegrators() });
    ompLF.execute(ompVec, { fes, lf.getIntegrators() });
    expectSameStandaloneVector(seqVec, ompVec);

    ierr = VecDestroy(&seqVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = VecDestroy(&ompVec);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v)) + Integral(u, v);

    using BFType = decltype(bf);
    Mat seqMat = nullptr;
    Mat ompMat = nullptr;
    ierr = MatCreate(PETSC_COMM_SELF, &seqMat);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatCreate(PETSC_COMM_SELF, &ompMat);
    ASSERT_EQ(ierr, PETSC_SUCCESS);

    PETSc::Assembly::Sequential<::Mat, BFType> seqBF;
    PETSc::Assembly::OpenMP<::Mat, BFType> ompBF;
    seqBF.execute(seqMat, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });
    ompBF.execute(ompMat, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });
    expectSameStandaloneMatrix(seqMat, ompMat);

    ierr = MatDestroy(&seqMat);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
    ierr = MatDestroy(&ompMat);
    ASSERT_EQ(ierr, PETSC_SUCCESS);
  }
#endif
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

  TEST(PETSc_Form, SequentialAffineIdentificationHasDefect)
  {
    checkPETScAffineIdentificationDefect<PETSc::Assembly::Sequential>();
  }

  TEST(PETSc_Form, SequentialSelfIdentificationMatchesZeroValueConstraint)
  {
    checkPETScSelfIdentificationMatchesZeroValueConstraint<PETSc::Assembly::Sequential>();
  }

  TEST(PETSc_Form, SequentialReassemblyKeepsExplicitZeroStructuralEntries)
  {
    checkPETScReassemblyKeepsExplicitZeroStructuralEntries<PETSc::Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  TEST(PETSc_Form, OpenMPStandaloneFormsMatchSequential)
  {
    checkPETScStandaloneOpenMPFormsMatchSequential();
  }

  TEST(PETSc_Form, OpenMPVectorMasterIdentificationProjectsMatrixAndVector)
  {
    checkPETScVectorMasterIdentification<PETSc::Assembly::OpenMP>();
  }

  TEST(PETSc_Form, OpenMPAffineIdentificationHasDefect)
  {
    checkPETScAffineIdentificationDefect<PETSc::Assembly::OpenMP>();
  }

  TEST(PETSc_Form, OpenMPSelfIdentificationMatchesZeroValueConstraint)
  {
    checkPETScSelfIdentificationMatchesZeroValueConstraint<PETSc::Assembly::OpenMP>();
  }

  TEST(PETSc_Form, OpenMPReassemblyKeepsExplicitZeroStructuralEntries)
  {
    checkPETScReassemblyKeepsExplicitZeroStructuralEntries<PETSc::Assembly::OpenMP>();
  }
#endif

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
