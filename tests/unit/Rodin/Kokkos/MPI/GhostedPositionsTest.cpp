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
    using Cache = ::Rodin::PETSc::Internal::GridFunctionAdapter<TestExecutionSpace>;
    static_assert(requires(const Cache& adapter) { adapter.getGhostedPositions(); });
    static_assert(
      !requires(const Cache& adapter, const FES& space) { adapter.setSpace(space); });
    using DataAccess = ::Rodin::Kokkos::Variational::GridFunctionData<::Vec>;
    using Retention = ::Rodin::PETSc::Internal::CoefficientCache<TestExecutionSpace>;
    const auto& execution = fes.getMesh().getContext().getExecutionSpace();
    Retention retention;
    const auto prepare = [&](const auto& space, ::Vec data) {
      Retention::Scope scope(retention, execution);
      return DataAccess::prepareView(space, data);
    };
    PETSc::Variational::GridFunction field(fes);
    field = Real(2);
    field.sync();
    const auto data = field.getData();
    const auto table = Cache::of(data).setSpace(fes).getGhostedPositions();
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
    EXPECT_EQ(prepared.getGhostedPositions().positions.data(), table.positions.data());
    EXPECT_EQ(Cache::of(data).setSpace(fes).getGhostedPositions().positions.data(),
      table.positions.data());
    const auto original = prepare(fes, data);
    EXPECT_EQ(prepare(fes, data).data(), original.data());
    field = Real(3);
    field.sync();
    EXPECT_EQ(prepared.getGhostedPositions().positions.data(), table.positions.data());
    EXPECT_EQ(Cache::of(data).setSpace(fes).getGhostedPositions().positions.data(),
      table.positions.data());
    const auto updated = prepare(fes, data);
    const auto updatedHost =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, updated);
    for (size_t i = 0; i < updatedHost.extent(0); ++i)
      EXPECT_EQ(updatedHost(i), Real(3));
    const auto oldHost =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, original);
    for (size_t i = 0; i < oldHost.extent(0); ++i)
      EXPECT_EQ(oldHost(i), Real(2));

    // Changing just a ghost entry refreshes the gather with the same numbering.
    PetscInt owned = 0, count = 0;
    ::Vec form = nullptr;
    ASSERT_EQ(VecGetLocalSize(data, &owned), PETSC_SUCCESS);
    ASSERT_EQ(VecGhostGetLocalForm(data, &form), PETSC_SUCCESS);
    ASSERT_EQ(VecGetLocalSize(form, &count), PETSC_SUCCESS);
    if (owned < count)
    {
      PetscScalar* values = nullptr;
      ASSERT_EQ(VecGetArray(form, &values), PETSC_SUCCESS);
      values[owned] = 7;
      ASSERT_EQ(VecRestoreArray(form, &values), PETSC_SUCCESS);
    }
    ASSERT_EQ(VecGhostRestoreLocalForm(data, &form), PETSC_SUCCESS);
    const auto withGhostWrite = prepare(fes, data);
    if (owned < count)
    {
      EXPECT_NE(withGhostWrite.data(), updated.data());
    }
    EXPECT_EQ(prepared.getGhostedPositions().positions.data(), table.positions.data());
    const auto written =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, withGhostWrite);
    for (size_t i = 0; i < written.extent(0); ++i)
      EXPECT_EQ(written(i), owned < count && host(i) == owned ? Real(7) : Real(3));

    // VecDuplicate can inherit the source vector's composed numbering adapter.
    auto another = field;
    another = Real(4);
    another.sync();
    const auto fresh = prepare(fes, another.getData());
    EXPECT_EQ(prepare(fes, another.getData()).data(), fresh.data());
    const auto separate =
      Cache::of(another.getData()).setSpace(fes).getGhostedPositions();
    EXPECT_NE(separate.positions.data(), table.positions.data());

    FES copied(fes);
    EXPECT_NE(copied.getRevision(), fes.getRevision());
    const auto copiedTable = Cache::of(data).setSpace(copied).getGhostedPositions();
    EXPECT_NE(copiedTable.positions.data(), table.positions.data());
    FES moved(std::move(copied));
    EXPECT_NE(moved.getRevision(), copied.getRevision());
    const auto movedTable = Cache::of(data).setSpace(moved).getGhostedPositions();
    EXPECT_NE(movedTable.positions.data(), copiedTable.positions.data());

    // Reconstruct at exactly the same address with identical vector states.
    // Coefficient copies must not be keyed by the address of this space.
    Optional<FES> reused(std::in_place, fes);
    const auto* address = &*reused;
    const auto revision = reused->getRevision();
    const auto firstTable = Cache::of(data).setSpace(*reused).getGhostedPositions();
    const auto firstValues = prepare(*reused, data);
    reused.reset();
    reused.emplace(fes);
    EXPECT_EQ(&*reused, address);
    EXPECT_NE(reused->getRevision(), revision);
    const auto nextTable = Cache::of(data).setSpace(*reused).getGhostedPositions();
    EXPECT_NE(firstTable.positions.data(), nextTable.positions.data());
    const auto nextValues = prepare(*reused, data);
    EXPECT_NE(firstValues.data(), nextValues.data());

    typename Cache::Table retained{};
    {
      PETSc::Variational::GridFunction temporary(fes);
      temporary = Real(1);
      temporary.sync();
      retained = Cache::of(temporary.getData()).setSpace(fes).getGhostedPositions();
      EXPECT_EQ(retained.positions.use_count(), 2);
    }
    EXPECT_EQ(retained.positions.use_count(), 1);
    // Earlier borrows remain valid after all cache replacements.
    const auto after =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, table.positions);
    for (size_t i = 0; i < after.extent(0); ++i)
      EXPECT_EQ(after(i), reference(i));
  }

  TEST(Kokkos_MPI_PETSc, GatherRetentionIsBounded)
  {
    using Cache = ::Rodin::PETSc::Internal::CoefficientCache<TestExecutionSpace>;
    using Data = ::Rodin::Kokkos::Variational::GridFunctionData<::Vec>;
    Context::MPI context(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(context, Polytope::Type::Triangle, 6);
    const DeviceMesh mesh(DeviceContext(context), host);
    const Variational::P1 fes(mesh);
    const auto& execution = mesh.getContext().getExecutionSpace();
    const size_t bytes = fes.getLocalSize() * sizeof(Real);
    ASSERT_GT(bytes, 0);
    Cache cache(bytes);
    PETSc::Variational::GridFunction first(fes), second(fes);
    first = Real(1);
    first.sync();
    second = Real(2);
    second.sync();
    const Real* idle = nullptr;
    {
      Cache::Scope scope(cache, execution);
      const auto active = Data::prepareView(fes, first.getData());
      const auto temporary = Data::prepareView(fes, second.getData());
      EXPECT_NE(active.data(), temporary.data());
      EXPECT_EQ(Data::prepareView(fes, second.getData()).data(), temporary.data());
      EXPECT_EQ(cache.getRetainedBytes(), bytes);
      EXPECT_EQ(::Rodin::Kokkos::Internal::read(execution, active, 0), Real(1));
      EXPECT_EQ(::Rodin::Kokkos::Internal::read(execution, temporary, 0), Real(2));
      idle = active.data();
    }
    {
      Cache::Scope scope(cache, execution);
      const auto reused = Data::prepareView(fes, second.getData());
      EXPECT_EQ(reused.data(), idle);
      EXPECT_EQ(cache.getRetainedBytes(), bytes);
      EXPECT_EQ(::Rodin::Kokkos::Internal::read(execution, reused, 0), Real(2));
    }
    second = Real(3);
    second.sync();
    {
      Cache::Scope scope(cache, execution);
      const auto refreshed = Data::prepareView(fes, second.getData());
      EXPECT_EQ(refreshed.data(), idle);
      EXPECT_EQ(::Rodin::Kokkos::Internal::read(execution, refreshed, 0), Real(3));
    }
    Cache tooSmall(bytes - 1);
    {
      Cache::Scope scope(tooSmall, execution);
      const auto temporary = Data::prepareView(fes, first.getData());
      EXPECT_EQ(Data::prepareView(fes, first.getData()).data(), temporary.data());
      EXPECT_EQ(::Rodin::Kokkos::Internal::read(execution, temporary, 0), Real(1));
      EXPECT_EQ(tooSmall.getRetainedBytes(), 0);
    }
    EXPECT_EQ(tooSmall.getRetainedBytes(), 0);
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
