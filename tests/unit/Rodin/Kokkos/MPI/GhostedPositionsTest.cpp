/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks adapter cache ownership and numbering revisions for one space.
   * @tparam FES Distributed device space.
   * @param fes Space whose mesh remains alive throughout the test.
   */
  template <class FES>
  void checkGhostedPositionCache(const FES& fes)
  {
    using Cache = ::Rodin::PETSc::Internal::GhostedPositions<TestMemorySpace>;
    static_assert(requires(const Cache& adapter) { adapter.getView(); });
    static_assert(
      !requires(const Cache& adapter, const FES& space) { adapter.setSpace(space); });
    using Coefficients = ::Rodin::Kokkos::Variational::GridFunctionCoefficients<::Vec>;
    PETSc::Variational::GridFunction field(fes);
    field = Real(2);
    field.sync();
    const auto data = field.getData();
    const auto table = Cache::of(data).setSpace(fes).getView();
    const auto expected = fes.getGhostedPositions();
    EXPECT_EQ(table.size, expected.size);
    auto host =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, table.positions);
    const auto reference =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, expected.positions);
    ASSERT_EQ(host.extent(0), reference.extent(0));
    for (size_t i = 0; i < host.extent(0); ++i)
      EXPECT_EQ(host(i), reference(i));
    const auto& prepared = Cache::of(data);
    EXPECT_EQ(prepared.getView().positions.data(), table.positions.data());
    EXPECT_EQ(
      Cache::of(data).setSpace(fes).getView().positions.data(), table.positions.data());
    const auto original = Coefficients::getView(fes, data);
    EXPECT_EQ(Coefficients::getView(fes, data).data(), original.data());
    field = Real(3);
    field.sync();
    EXPECT_EQ(prepared.getView().positions.data(), table.positions.data());
    EXPECT_EQ(
      Cache::of(data).setSpace(fes).getView().positions.data(), table.positions.data());
    const auto updated = Coefficients::getView(fes, data);
    const auto updatedHost =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, updated);
    for (size_t i = 0; i < updatedHost.extent(0); ++i)
      EXPECT_EQ(updatedHost(i), Real(3));
    const auto oldHost =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, original);
    for (size_t i = 0; i < oldHost.extent(0); ++i)
      EXPECT_EQ(oldHost(i), Real(2));

    // Each PETSc vector owns its own map; destroying it must drop its cache.
    PETSc::Variational::GridFunction another(fes);
    another = Real(4);
    another.sync();
    const auto fresh = Coefficients::getView(fes, another.getData());
    EXPECT_EQ(Coefficients::getView(fes, another.getData()).data(), fresh.data());
    const auto separate = Cache::of(another.getData()).setSpace(fes).getView();
    EXPECT_NE(separate.positions.data(), table.positions.data());

    FES copied(fes);
    EXPECT_NE(copied.getNumberingRevision(), fes.getNumberingRevision());
    const auto copiedTable = Cache::of(data).setSpace(copied).getView();
    EXPECT_NE(copiedTable.positions.data(), table.positions.data());
    FES moved(std::move(copied));
    EXPECT_NE(moved.getNumberingRevision(), copied.getNumberingRevision());
    const auto movedTable = Cache::of(data).setSpace(moved).getView();
    EXPECT_NE(movedTable.positions.data(), copiedTable.positions.data());

    // Reconstruct at exactly the same address with identical vector states.
    // Coefficient copies must not be keyed by the address of this space.
    Optional<FES> reused(std::in_place, fes);
    const auto* address = &*reused;
    const auto revision = reused->getNumberingRevision();
    const auto firstTable = Cache::of(data).setSpace(*reused).getView();
    const auto firstValues = Coefficients::getView(*reused, data);
    reused.reset();
    reused.emplace(fes);
    EXPECT_EQ(&*reused, address);
    EXPECT_NE(reused->getNumberingRevision(), revision);
    const auto nextTable = Cache::of(data).setSpace(*reused).getView();
    EXPECT_NE(firstTable.positions.data(), nextTable.positions.data());
    const auto nextValues = Coefficients::getView(*reused, data);
    EXPECT_NE(firstValues.data(), nextValues.data());

    typename Cache::Table retained{};
    {
      PETSc::Variational::GridFunction temporary(fes);
      temporary = Real(1);
      temporary.sync();
      retained = Cache::of(temporary.getData()).setSpace(fes).getView();
      EXPECT_EQ(retained.positions.use_count(), 2);
    }
    EXPECT_EQ(retained.positions.use_count(), 1);
    // Earlier borrows remain valid after all cache replacements.
    const auto after =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, table.positions);
    for (size_t i = 0; i < after.extent(0); ++i)
      EXPECT_EQ(after(i), reference(i));
  }

  /// @brief P1, H1 and P0g retain adapter maps without mutating their spaces.
  TEST(Kokkos_MPI_PETSc, GhostedPositionCacheOwnership)
  {
    Context::MPI context(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(context, Polytope::Type::Triangle, 6);
    const DeviceMesh mesh(DeviceContext(context), host);
    const Variational::P1 scalar(mesh);
    const Variational::H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    const Variational::P0g constant(mesh);
    checkGhostedPositionCache(scalar);
    checkGhostedPositionCache(quadratic);
    checkGhostedPositionCache(constant);
  }
}
