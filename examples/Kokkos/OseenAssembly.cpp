/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Compares the device form language with a hand written Oseen kernel.
 *
 * Usage: `KokkosOseenAssembly [mesh.mesh | n] [repetitions] [renumber]`. A
 * number n builds an n x n x n tetrahedral grid, otherwise a MEDIT mesh is
 * loaded, whose cells must all be tetrahedra since the hand written kernel
 * knows no other cell. The word renumber orders the vertices and the cells of
 * the mesh along space before the assembly, see Mesh::renumber.
 *
 * Both sides assemble the same linearized Navier-Stokes (Oseen) step with
 * P1/P1 elements and a pressure stabilization,
 *
 *   rho/dt (u, v) + mu (Jac u, Jac v) + rho (Jac(u) w, v)
 *   + rho/2 (div(w) u, v) - (p, div v) + (div u, q) + eps (grad p, grad q)
 *   = rho/dt (w, v),
 *
 * into a PETSc Kokkos matrix and vector. The form language goes through
 * Problem and the device COO assembly, one kernel per integrator. The hand
 * written version computes the 16 x 16 local matrix of each cell in a single
 * kernel and hands the values to MatSetValuesCOO and VecSetValuesCOO.
 */
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Kokkos/PETSc.h>

namespace Rodin::Examples
{
  using namespace Rodin::Variational;

#if defined(KOKKOS_ENABLE_CUDA)
  using Backend = Context::CUDA;

#elif defined(KOKKOS_ENABLE_HIP)
  using Backend = Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
  using Backend = Context::SYCL;
#else
  using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

  using DeviceContext = Context::Kokkos<Context::Local, Backend>;
  using HostMesh = Geometry::Mesh<Context::Local>;
  using DeviceMesh = Geometry::Mesh<DeviceContext>;
  using MemorySpace = DeviceContext::MemorySpaceType;
  using ExecutionSpace = DeviceContext::ExecutionSpaceType;

  /// @brief Physical and stabilization parameters of the Oseen step.
  struct Parameters
  {
    Real rho = 1060;
    Real dt = 1e-3;
    Real mu = 0.0049;
    Real eps = 1e-3;
  };

  /// @brief Timing summary in milliseconds.
  struct Timing
  {
    double mean;
    double min;
  };

  /**
   * @brief Times a device operation, fencing before and after each call.
   * @param f Operation.
   * @param repetitions Number of timed calls, after two warm-up calls.
   */
  template <class F>
  Timing measure(F&& f, int repetitions)
  {
    using Clock = std::chrono::steady_clock;
    for (int i = 0; i < 2; ++i)
      f();
    std::vector<double> times;
    for (int i = 0; i < repetitions; ++i)
    {
      ::Kokkos::fence();
      const auto t0 = Clock::now();
      f();
      ::Kokkos::fence();
      times.push_back(
        std::chrono::duration<double, std::milli>(Clock::now() - t0).count());
    }
    double sum = 0;
    for (double t : times)
      sum += t;
    return {sum / times.size(), *std::min_element(times.begin(), times.end())};
  }

  /**
   * @brief Hand written assembly of the Oseen step.
   *
   * The coordinate list has 256 matrix entries and 16 vector entries per cell,
   * in the DOF numbering of Rodin: velocity component i of vertex a is
   * a + i * nv, and the pressure of vertex a is 3 * nv + a.
   */
  struct HandAssembly
  {
    using RealView = ::Kokkos::View<Real*, MemorySpace>;
    using IndexView = ::Kokkos::View<int*, MemorySpace>;

    int nv;
    int nt;
    RealView x;
    IndexView tets;
    RealView w;
    RealView values;
    RealView rhs;
    ::Mat A = nullptr;
    ::Vec b = nullptr;

    /// @param mesh Host mesh of tetrahedra.
    explicit HandAssembly(const HostMesh& mesh)
      : nv(static_cast<int>(mesh.getVertexCount())),
        nt(static_cast<int>(mesh.getCellCount()))
    {
      std::vector<Real> hx(3 * size_t(nv));
      for (int v = 0; v < nv; ++v)
      {
        const auto c = mesh.getVertexCoordinates(v);
        for (int i = 0; i < 3; ++i)
          hx[3 * v + i] = c(i);
      }
      std::vector<int> ht(4 * size_t(nt));
      for (int e = 0; e < nt; ++e)
      {
        const auto& key = mesh.getConnectivity().getPolytope(3, e);
        for (int a = 0; a < 4; ++a)
          ht[4 * e + a] = static_cast<int>(key(a));
      }
      x = upload<Real>(hx, "x");
      tets = upload<int>(ht, "tets");
      w = RealView("w", 3 * size_t(nv));
      values = RealView("values", 256 * size_t(nt));
      rhs = RealView("rhs", 16 * size_t(nt));

      const PetscInt n = 4 * nv;
      std::vector<PetscInt> rows(256 * size_t(nt)), cols(256 * size_t(nt));
      std::vector<PetscInt> entries(16 * size_t(nt));
      for (int e = 0; e < nt; ++e)
      {
        PetscInt dofs[16];
        for (int r = 0; r < 16; ++r)
        {
          dofs[r] =
            r < 12 ? ht[4 * e + r / 3] + (r % 3) * nv : 3 * nv + ht[4 * e + r - 12];
        }
        for (int r = 0; r < 16; ++r)
        {
          entries[16 * size_t(e) + r] = dofs[r];
          for (int c = 0; c < 16; ++c)
          {
            rows[256 * size_t(e) + 16 * r + c] = dofs[r];
            cols[256 * size_t(e) + 16 * r + c] = dofs[c];
          }
        }
      }
      MatCreate(PETSC_COMM_SELF, &A);
      MatSetSizes(A, n, n, n, n);
      MatSetType(A, MATAIJKOKKOS);
      MatSetPreallocationCOO(
        A, static_cast<PetscCount>(rows.size()), rows.data(), cols.data());
      VecCreate(PETSC_COMM_SELF, &b);
      VecSetSizes(b, n, n);
      VecSetType(b, VECKOKKOS);
      VecSetPreallocationCOO(
        b, static_cast<PetscCount>(entries.size()), entries.data());
    }

    ~HandAssembly()
    {
      MatDestroy(&A);
      VecDestroy(&b);
    }

    HandAssembly(const HandAssembly&) = delete;
    HandAssembly& operator=(const HandAssembly&) = delete;

    /**
     * @brief Copies a host array to the device.
     * @param data Host values.
     * @param label View label.
     */
    template <class T>
    static ::Kokkos::View<T*, MemorySpace> upload(
      const std::vector<T>& data, const char* label)
    {
      ::Kokkos::View<T*, MemorySpace> view(std::string(label), data.size());
      ::Kokkos::deep_copy(
        view, ::Kokkos::View<const T*, ::Kokkos::HostSpace>(data.data(), data.size()));
      return view;
    }

    /**
     * @brief Computes the local matrices and vectors of all cells.
     * @param P Parameters.
     * @param nv Vertex count.
     * @param nt Cell count.
     * @param x Vertex coordinates.
     * @param tets Cell vertices.
     * @param w Advection field, component-blocked.
     * @param values Matrix entries.
     * @param rhs Vector entries.
     */
    static void kernel(const Parameters P, int nv, int nt, RealView x, IndexView tets,
      RealView w, RealView values, RealView rhs)
    {
      ::Kokkos::parallel_for(
        "Hand Oseen cells", ::Kokkos::RangePolicy<ExecutionSpace>(0, nt),
        KOKKOS_LAMBDA(const int e) {
          int v[4];
          Real p[4][3], wa[4][3];
          for (int a = 0; a < 4; ++a)
          {
            v[a] = tets(4 * e + a);
            for (int i = 0; i < 3; ++i)
            {
              p[a][i] = x(3 * v[a] + i);
              wa[a][i] = w(v[a] + i * nv);
            }
          }
          Real J[3][3];
          for (int i = 0; i < 3; ++i)
          {
            for (int k = 0; k < 3; ++k)
              J[i][k] = p[k + 1][i] - p[0][i];
          }
          const Real det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
            J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
            J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
          Real G[4][3];
          // Rows of J^{-1} are the barycentric gradients of vertices 1..3.
          G[1][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
          G[1][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det;
          G[1][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
          G[2][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det;
          G[2][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
          G[2][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
          G[3][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
          G[3][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det;
          G[3][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
          for (int i = 0; i < 3; ++i)
            G[0][i] = -(G[1][i] + G[2][i] + G[3][i]);
          const Real V = ::Kokkos::fabs(det) / 6;

          Real divw = 0, ws[3] = {};
          for (int a = 0; a < 4; ++a)
          {
            for (int i = 0; i < 3; ++i)
            {
              divw += G[a][i] * wa[a][i];
              ws[i] += wa[a][i];
            }
          }
          // Exact P1 integrals: (l_a, l_b) = V (1 + d_ab) / 20 and
          // (w.G_b, l_a) = V G_b.(ws + w_a) / 20.
          Real M[4][4], C[4][4], K[4][4];
          for (int a = 0; a < 4; ++a)
          {
            for (int b = 0; b < 4; ++b)
            {
              M[a][b] = V * (a == b ? 2 : 1) / 20;
              C[a][b] = 0;
              K[a][b] = 0;
              for (int i = 0; i < 3; ++i)
              {
                C[a][b] += V * G[b][i] * (ws[i] + wa[a][i]) / 20;
                K[a][b] += V * G[a][i] * G[b][i];
              }
            }
          }
          const Real react = P.rho / P.dt + Real(0.5) * P.rho * divw;
          Real* out = &values(256 * size_t(e));
          for (int r = 0; r < 16; ++r)
          {
            for (int c = 0; c < 16; ++c)
            {
              Real val = 0;
              if (r < 12 && c < 12)
              {
                const int a = r / 3, i = r % 3, b = c / 3, j = c % 3;
                if (i == j)
                  val = react * M[a][b] + P.mu * K[a][b] + P.rho * C[a][b];
              }
              else if (r < 12)
                val = -V / 4 * G[r / 3][r % 3];
              else if (c < 12)
                val = V / 4 * G[c / 3][c % 3];
              else
                val = P.eps * K[r - 12][c - 12];
              out[16 * r + c] = val;
            }
          }
          Real* f = &rhs(16 * size_t(e));
          for (int a = 0; a < 4; ++a)
          {
            for (int i = 0; i < 3; ++i)
            {
              Real s = 0;
              for (int b = 0; b < 4; ++b)
                s += M[a][b] * wa[b][i];
              f[3 * a + i] = P.rho / P.dt * s;
            }
          }
          for (int a = 0; a < 4; ++a)
            f[12 + a] = 0;
        });
    }

    /**
     * @brief Assembles the system.
     * @param P Parameters.
     * @param data Advection field on the host, as stored by GridFunction.
     */
    void assemble(const Parameters& P, const Math::Vector<Real>& data)
    {
      ::Kokkos::deep_copy(
        w, ::Kokkos::View<const Real*, ::Kokkos::HostSpace>(data.data(), data.size()));
      kernel(P, nv, nt, x, tets, w, values, rhs);
      ::Kokkos::fence();
      MatSetValuesCOO(A, values.data(), INSERT_VALUES);
      VecSetValuesCOO(b, rhs.data(), INSERT_VALUES);
    }
  };

  /**
   * @brief Relative Frobenius distance between two matrices.
   * @param A Reference.
   * @param B Candidate.
   */
  Real distance(::Mat A, ::Mat B)
  {
    // The difference is taken on host copies: the device sum of Kokkos-Kernels
    // asks for more team scratch memory than RDNA GPUs have on large grids.
    ::Mat D = nullptr, H = nullptr;
    MatConvert(A, MATSEQAIJ, MAT_INITIAL_MATRIX, &D);
    MatConvert(B, MATSEQAIJ, MAT_INITIAL_MATRIX, &H);
    MatAXPY(D, -1, H, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&H);
    PetscReal d = 0, a = 0;
    MatNorm(D, NORM_FROBENIUS, &d);
    MatNorm(A, NORM_FROBENIUS, &a);
    MatDestroy(&D);
    return d / a;
  }

  /**
   * @brief Relative Euclidean distance between two vectors.
   * @param a Reference.
   * @param b Candidate.
   */
  Real distance(::Vec a, ::Vec b)
  {
    ::Vec d = nullptr;
    VecDuplicate(a, &d);
    VecWAXPY(d, -1, b, a);
    PetscReal nd = 0, na = 0;
    VecNorm(d, NORM_2, &nd);
    VecNorm(a, NORM_2, &na);
    VecDestroy(&d);
    return nd / na;
  }

  /**
   * @brief Loads the tetrahedra of a MEDIT mesh, or builds a grid.
   * @param source Path, or the number of cells per side.
   */
  HostMesh load(const std::string& source)
  {
    if (std::all_of(
          source.begin(), source.end(), [](char c) { return c >= '0' && c <= '9'; }))
    {
      const size_t n = std::stoul(source) + 1;
      return HostMesh::UniformGrid(Geometry::Polytope::Type::Tetrahedron, {n, n, n});
    }
    HostMesh mesh;
    mesh.load(source, IO::FileFormat::MEDIT);
    return mesh;
  }

  /**
   * @brief Runs the comparison.
   * @param argc Argument count.
   * @param argv Arguments.
   */
  int run(int argc, char** argv)
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    {
      const std::string source = argc > 1 ? argv[1] : "16";
      const int repetitions = argc > 2 ? std::atoi(argv[2]) : 20;
      const bool renumber = argc > 3 && std::string(argv[3]) == "renumber";
      const Parameters P;

      const HostMesh host = renumber ? load(source).renumber() : load(source);
      if (host.getPolytopeCount(Geometry::Polytope::Type::Tetrahedron) !=
        host.getCellCount())
      {
        std::fprintf(stderr, "%s has cells other than tetrahedra.\n", source.c_str());
        PetscFinalize();
        return EXIT_FAILURE;
      }
      const DeviceMesh mesh(DeviceContext{}, host);
      std::printf("mesh: %zu vertices, %zu tetrahedra, %zu unknowns\n",
        host.getVertexCount(), host.getCellCount(), 4 * host.getVertexCount());

      P1 hvh(host, 3);
      GridFunction hw(hvh);
      hw = VectorFunction{[](const Geometry::Point& p) { return std::sin(p.y()) + 0.5; },
        [](const Geometry::Point& p) { return std::cos(p.x()) * p.z(); },
        [](const Geometry::Point& p) { return p.x() - p.y(); }};

      P1 vh(mesh, 3);
      P1 qh(mesh);
      GridFunction w(vh);
      w.getData() = hw.getData();
      PETSc::Variational::TrialFunction u(vh);
      PETSc::Variational::TrialFunction p(qh);
      PETSc::Variational::TestFunction v(vh);
      PETSc::Variational::TestFunction q(qh);

      Problem oseen(u, p, v, q);
      oseen = Integral(P.rho / P.dt * u, v) + Integral(P.mu * Jacobian(u), Jacobian(v)) +
        Integral(P.rho * Jacobian(u) * w, v) + Integral(0.5 * P.rho * Div(w) * u, v) -
        Integral(p, Div(v)) + Integral(Div(u), q) + Integral(P.eps * Grad(p), Grad(q)) -
        Integral(P.rho / P.dt * w, v);

      // The hand written assembly is built after the first ones of the form
      // language, so that its arrays do not share the device memory with the
      // sort of the coordinate pattern.
      const Timing dsl = measure([&] { oseen.assemble(); }, repetitions);
      HandAssembly hand(host);
      const Timing manual = measure([&] { hand.assemble(P, hw.getData()); }, repetitions);

      auto& axb = oseen.getLinearSystem();
      std::printf("relative distance: matrix %.3e, vector %.3e\n",
        distance(hand.A, axb.getOperator()), distance(hand.b, axb.getVector()));
      std::printf("%-14s %10s %10s\n", "assembly (ms)", "mean", "min");
      std::printf("%-14s %10.3f %10.3f\n", "form language", dsl.mean, dsl.min);
      std::printf("%-14s %10.3f %10.3f\n", "hand written", manual.mean, manual.min);
      std::printf(
        "ratio %.2f (mean), %.2f (min)\n", dsl.mean / manual.mean, dsl.min / manual.min);
    }
    PetscFinalize();
    return 0;
  }
} // namespace Rodin::Examples

int main(int argc, char** argv)
{
  return Rodin::Examples::run(argc, argv);
}
