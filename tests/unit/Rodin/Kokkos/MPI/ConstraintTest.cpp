/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file ConstraintTest.cpp
 * @brief Periodic and identification conditions of distributed device
 * problems.
 *
 * The host distributed assembly does not consume periodic conditions, so
 * the periodic problems are compared node by node with the same problem
 * solved sequentially on the whole grid. The identification conditions are
 * compared with the host distributed assembly.
 */

#include "Common.h"

namespace Rodin::Tests::Unit
{
  /// @brief Grid node, by its integer coordinates.
  using Node = std::pair<size_t, size_t>;

  /**
   * @brief Distributed unit grid, its device space and a sequential copy,
   * with the global and sequential degree of freedom of each node.
   */
  struct Grid
  {
    using SequentialMesh = Geometry::Mesh<Context::Local>;

    /**
     * @brief Distributes the grid and numbers its nodes; collective.
     * @param ctx MPI context. @param type Cell type. @param size Nodes per
     * side.
     */
    Grid(const Context::MPI& ctx, Polytope::Type type, size_t size)
      : n(size),
        host(distributeFromRoot(ctx, type, size)),
        mesh(DeviceContext(ctx), host),
        fes(mesh),
        seq(sequential(type, size)),
        sfes(seq)
    {
      const auto& shard = host.getShard();
      std::vector<size_t> held;
      for (auto it = shard.getVertex(); !it.end(); ++it)
      {
        const auto node = position(it->getCoordinates());
        const Index g = fes.getGlobalIndex(fes.getShard().getDOFs(0, it->getIndex())[0]);
        nodes[g] = node;
        held.insert(held.end(), {node.first, node.second, static_cast<size_t>(g)});
      }
      std::vector<std::vector<size_t>> gathered;
      boost::mpi::all_gather(ctx.getCommunicator(), held, gathered);
      for (const auto& h : gathered)
      {
        for (size_t k = 0; k < h.size(); k += 3)
          globals[{h[k], h[k + 1]}] = static_cast<Index>(h[k + 2]);
      }
      for (auto it = seq.getVertex(); !it.end(); ++it)
        seqDOFs[position(it->getCoordinates())] = sfes.getDOFs(0, it->getIndex())[0];
    }

    /// @brief Integer coordinates of a point of the grid.
    template <class Point>
    Node position(const Point& x) const
    {
      return {static_cast<size_t>(std::llround(x.x() * Real(n - 1))),
        static_cast<size_t>(std::llround(x.y() * Real(n - 1)))};
    }

    /**
     * @brief Periodic pairs of the slaves the process holds, in global
     * indices, and of all the slaves, in sequential indices.
     * @param master Master node of a node, or the node itself when it is
     * not a slave.
     */
    template <class Master>
    std::pair<IndexMap<IndexSet>, IndexMap<IndexSet>> pairs(Master&& master) const
    {
      IndexMap<IndexSet> distributed, local;
      for (const auto& [g, node] : nodes)
      {
        const Node m = master(node);
        if (m != node)
          distributed[g].insert(globals.at(m));
      }
      for (const auto& [node, dof] : seqDOFs)
      {
        const Node m = master(node);
        if (m != node)
          local[dof].insert(seqDOFs.at(m));
      }
      return {std::move(distributed), std::move(local)};
    }

    /**
     * @brief Expects a distributed solution to match the sequential one at
     * each node; collective.
     * @param x Distributed solution of a field on the grid.
     * @param ref Sequential solution of the field.
     */
    void expectSame(::Vec x, const Math::Vector<Real>& ref) const
    {
      PetscInt begin = 0, end = 0;
      VecGetOwnershipRange(x, &begin, &end);
      const PetscScalar* values = nullptr;
      VecGetArrayRead(x, &values);
      int compared = 0;
      for (const auto& [g, node] : nodes)
      {
        const PetscInt global = static_cast<PetscInt>(g);
        if (global < begin || global >= end)
          continue;
        EXPECT_NEAR(values[global - begin], ref(seqDOFs.at(node)), Real(1e-9));
        ++compared;
      }
      VecRestoreArrayRead(x, &values);
      int total = 0;
      MPI_Allreduce(&compared, &total, 1, MPI_INT, MPI_SUM, *g_world);
      EXPECT_EQ(static_cast<size_t>(total), n * n);
    }

    static SequentialMesh sequential(Polytope::Type type, size_t size)
    {
      auto res = SequentialMesh::UniformGrid(type, {size, size});
      res.scale(Real(1) / static_cast<Real>(size - 1));
      res.getConnectivity().compute(1, 2);
      return res;
    }

    const size_t n;
    const HostMesh host;
    const DeviceMesh mesh;
    Variational::P1<Real, DeviceMesh> fes;
    SequentialMesh seq;
    Variational::P1<Real, SequentialMesh> sfes;

    /// @brief Node of each global degree of freedom the process holds.
    std::map<Index, Node> nodes;

    /// @brief Global degree of freedom of each node.
    std::map<Node, Index> globals;

    /// @brief Sequential degree of freedom of each node.
    std::map<Node, Index> seqDOFs;
  };

  /// @brief Host source term of the periodic problems.
  inline Real source(const Geometry::Point& p)
  {
    return 1 + std::sin(2 * M_PI * p.x()) * std::cos(2 * M_PI * p.y()) + p.x() * p.y();
  }

  /// @brief Device source term of the periodic problems.
  KOKKOS_INLINE_FUNCTION Real source(const DeviceGeometry::PointView& p)
  {
    return 1 + ::Kokkos::sin(2 * M_PI * p.x()) * ::Kokkos::cos(2 * M_PI * p.y()) +
      p.x() * p.y();
  }

  /**
   * @brief A doubly periodic problem matches the sequential host solution.
   *
   * The periodic pairs relate global degrees of freedom. Each process
   * gives those of the slaves it holds, whose masters another process may
   * hold.
   */
  void checkPeriodic(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    Grid grid(ctx, type, 13);
    const size_t n = grid.n;
    const auto [pairs, seqPairs] = grid.pairs([n](const Node& node) {
      return Node(node.first % (n - 1), node.second % (n - 1));
    });

    RealFunction hf([](const Geometry::Point& p) { return source(p); });
    TrialFunction su(grid.sfes);
    TestFunction sv(grid.sfes);
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) - Integral(hf, sv) +
      PeriodicBC(su, seqPairs);
    Solver::SparseLU(sp).solve();

    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return source(p);
    });
    PETSc::Variational::TrialFunction du(grid.fes);
    PETSc::Variational::TestFunction dv(grid.fes);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      PeriodicBC(du, pairs);

    solveDirect(dp);
    grid.expectSame(du.getSolution().getData(), su.getSolution().getData());
    // A second solve reuses the cached pattern and coordinates.
    solveDirect(dp);
    grid.expectSame(du.getSolution().getData(), su.getSolution().getData());
  }

  /**
   * @brief Periodic pairs combined with Dirichlet values match the
   * sequential host solution.
   *
   * The boundary is fixed. Interior slaves of the column before the last
   * one are identified with interior masters, with fixed boundary masters,
   * which fix them once flattened, or with other slaves, which chains the
   * identifications.
   */
  void checkPeriodicDirichlet(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    Grid grid(ctx, type, 13);
    const size_t n = grid.n;
    const auto [pairs, seqPairs] = grid.pairs([n](const Node& node) {
      const auto [i, j] = node;
      if (j == 0 || j == n - 1)
        return node;
      if (i == n - 2)
        return Node(j % 2 ? 1 : 0, j);
      if (i == n - 3 && j % 3 == 0)
        return Node(n - 2, j);
      return node;
    });

    RealFunction hf([](const Geometry::Point& p) { return source(p); });
    RealFunction hg([](const Geometry::Point& p) { return 1 + p.x() + 2 * p.y(); });
    TrialFunction su(grid.sfes);
    TestFunction sv(grid.sfes);
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) - Integral(hf, sv) +
      DirichletBC(su, hg) + PeriodicBC(su, seqPairs);
    Solver::SparseLU(sp).solve();

    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return source(p);
    });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() + 2 * p.y();
    });
    PETSc::Variational::TrialFunction du(grid.fes);
    PETSc::Variational::TestFunction dv(grid.fes);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg) + PeriodicBC(du, pairs);

    solveDirect(dp);
    grid.expectSame(du.getSolution().getData(), su.getSolution().getData());
    solveDirect(dp);
    grid.expectSame(du.getSolution().getData(), su.getSolution().getData());
  }

  /**
   * @brief A doubly periodic second field, coupled to a fixed first one,
   * matches the sequential host solution.
   *
   * The periodic pairs are global indices of the second field, which the
   * system offsets by the degrees of freedom of the first one.
   */
  void checkPeriodicSecondField(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    Grid grid(ctx, type, 13);
    const size_t n = grid.n;
    const auto [pairs, seqPairs] = grid.pairs([n](const Node& node) {
      return Node(node.first % (n - 1), node.second % (n - 1));
    });
    P1 sq(grid.seq);
    P1 dq(grid.mesh);
    EXPECT_EQ(dq.getSize(), grid.fes.getSize());

    RealFunction hf([](const Geometry::Point& p) { return source(p); });
    TrialFunction sa(grid.sfes);
    TestFunction sb(grid.sfes);
    TrialFunction sp(sq);
    TestFunction sr(sq);
    Problem sys(sa, sp, sb, sr);
    sys = Integral(Grad(sa), Grad(sb)) + Integral(sa, sb) - Integral(hf, sb) +
      Integral(Grad(sp), Grad(sr)) + Integral(sp, sr) + Integral(Real(0.5) * sa, sr) -
      Integral(hf, sr) + DirichletBC(sa, RealFunction(Real(0))) + PeriodicBC(sp, seqPairs);
    Solver::SparseLU(sys).solve();

    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return source(p);
    });
    PETSc::Variational::TrialFunction da(grid.fes);
    PETSc::Variational::TestFunction db(grid.fes);
    PETSc::Variational::TrialFunction dp(dq);
    PETSc::Variational::TestFunction dr(dq);
    Problem dev(da, dp, db, dr);
    dev = Integral(Grad(da), Grad(db)) + Integral(da, db) - Integral(df, db) +
      Integral(Grad(dp), Grad(dr)) + Integral(dp, dr) + Integral(Real(0.5) * da, dr) -
      Integral(df, dr) + DirichletBC(da, RealFunction(Real(0))) + PeriodicBC(dp, pairs);

    solveDirect(dev);
    grid.expectSame(da.getSolution().getData(), sa.getSolution().getData());
    grid.expectSame(dp.getSolution().getData(), sp.getSolution().getData());
    solveDirect(dev);
    grid.expectSame(da.getSolution().getData(), sa.getSolution().getData());
    grid.expectSame(dp.getSolution().getData(), sp.getSolution().getData());
  }

  /**
   * @brief An affine identification of the boundary of a field with the
   * other field matches the host distributed assembly.
   *
   * Each process assembles the identification of the boundary degrees of
   * freedom it holds, in global indices of both fields; the masters of the
   * second field may be held by other processes.
   */
  void checkIdentification(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 13);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1<Real, HostMesh> hah(host);
    P1<Real, HostMesh> hph(host);
    P1<Real, DeviceMesh> dah(mesh);
    P1<Real, DeviceMesh> dph(mesh);

    PETSc::Variational::TrialFunction ha(hah);
    PETSc::Variational::TestFunction hb(hah);
    PETSc::Variational::TrialFunction hp(hph);
    PETSc::Variational::TestFunction hq(hph);
    PETSc::Variational::TrialFunction da(dah);
    PETSc::Variational::TestFunction db(dah);
    PETSc::Variational::TrialFunction dp(dph);
    PETSc::Variational::TestFunction dq(dph);

    RealFunction hf([](const Geometry::Point& p) { return source(p); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return source(p);
    });

    Problem hs(ha, hp, hb, hq);
    hs = Integral(Grad(ha), Grad(hb)) + Integral(ha, hb) - Integral(hf, hb) +
      Integral(Grad(hp), Grad(hq)) + Integral(hp, hq) + Integral(Real(0.25) * ha, hq) -
      Integral(hf, hq) + DirichletBC(ha, Real(0.5) * hp, RealFunction(2));
    Problem ds(da, dp, db, dq);
    ds = Integral(Grad(da), Grad(db)) + Integral(da, db) - Integral(df, db) +
      Integral(Grad(dp), Grad(dq)) + Integral(dp, dq) + Integral(Real(0.25) * da, dq) -
      Integral(df, dq) + DirichletBC(da, Real(0.5) * dp, RealFunction(2));

    solveDirect(hs);
    solveDirect(ds);
    expectSameVector(ha.getSolution().getData(), da.getSolution().getData());
    expectSameVector(hp.getSolution().getData(), dp.getSolution().getData());
    solveDirect(ds);
    expectSameVector(ha.getSolution().getData(), da.getSolution().getData());
    expectSameVector(hp.getSolution().getData(), dp.getSolution().getData());
  }

  TEST(Kokkos_MPI_Periodic, Triangles)
  {
    checkPeriodic(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Periodic, Quadrilaterals)
  {
    checkPeriodic(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_PeriodicDirichlet, Triangles)
  {
    checkPeriodicDirichlet(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_PeriodicDirichlet, Quadrilaterals)
  {
    checkPeriodicDirichlet(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_PeriodicSecondField, Triangles)
  {
    checkPeriodicSecondField(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_PeriodicSecondField, Quadrilaterals)
  {
    checkPeriodicSecondField(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_Identification, Triangles)
  {
    checkIdentification(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Identification, Quadrilaterals)
  {
    checkIdentification(Polytope::Type::Quadrilateral);
  }
}
