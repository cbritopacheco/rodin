/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_SEQUENTIAL_H
#define RODIN_ASSEMBLY_SEQUENTIAL_H

#include <unordered_map>
#include <type_traits>
#include <algorithm>

#include "Rodin/Context/Local.h"

#include "Rodin/Tuple.h"

#include "Rodin/Math/Traits.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Region.h"

#include "Rodin/Variational/ForwardDecls.h"

#include "Rodin/Assembly/AssemblyBase.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  /**
   * @brief Sequential mesh iteration for single-threaded assembly.
   *
   * SequentialIteration provides a single-threaded iteration strategy over
   * mesh elements for assembly operations. This class encapsulates the logic
   * for iterating through mesh polytopes in a specified region in a 
   * deterministic, sequential order.
   *
   * @tparam MeshType Type of mesh to iterate over (specialized for local context)
   */
  template <>
  class SequentialIteration<Geometry::Mesh<Context::Local>>
  {
    public:
      /// @brief Mesh type for local execution context
      using MeshType = Geometry::Mesh<Context::Local>;

      /**
       * @brief Constructs sequential iteration over a mesh region.
       *
       * @param mesh Mesh to iterate over
       */
      SequentialIteration(const MeshType& mesh, const Geometry::Region&);

      /**
       * @brief Gets an iterator for traversing mesh elements.
       * @return PolytopeIterator for sequential mesh traversal
       */
      Geometry::PolytopeIterator getIterator() const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;  ///< Reference to the mesh
      Geometry::Region m_region;                      ///< Region to iterate over
  };

  /// @brief Template argument deduction guide for SequentialIteration
  SequentialIteration(
      const Geometry::Mesh<Context::Local>& mesh, const Geometry::Region&)
    -> SequentialIteration<Geometry::Mesh<Context::Local>>;
}

namespace Rodin::Assembly
{
  /**
   * @brief Sequential assembly implementation for linear forms.
   *
   * This class provides a single-threaded assembly implementation for linear
   * forms @f$ l(v) @f$, computing the discrete vector representation
   * @f$ b_i = l(\psi_i) @f$ by sequentially iterating through mesh elements
   * and accumulating contributions.
   *
   * @tparam FES Finite element space type
   *
   * ## Assembly Algorithm
   * For each element @f$ K @f$ in the mesh:
   * 1. Retrieve element quadrature rule and basis functions
   * 2. Compute local element vector @f$ b_K @f$
   * 3. Accumulate contributions to global vector @f$ b @f$
   *
   * The assembly is deterministic and reproducible, making it suitable for
   * debugging and verification purposes.
   */
  template <class FES>
  class Sequential<
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
    Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>> final
    : public AssemblyBase<
        Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>,
        Variational::LinearForm<FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using VectorType = Math::Vector<ScalarType>;

      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      using Parent = AssemblyBase<VectorType, LinearFormType>;

      using InputType = typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      void execute(VectorType& res, const InputType& input) const override
      {
        res.resize(input.getFES().getSize());
        res.setZero();
        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              lfi.setPolytope(*it);
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.getFES().getDOFs(d, i);
              for (size_t l = 0; l < static_cast<size_t>(dofs.size()); l++)
                res(dofs(l)) += lfi.integrate(l);
            }
          }
        }
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Sequential<
    Math::Matrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
      Solution, TrialFES, TestFES,
      Math::Matrix<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>> final
    : public AssemblyBase<
        Math::Matrix<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Variational::BilinearForm<
          Solution, TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = Math::Matrix<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setZero();
        const auto& mesh = input.getTrialFES().getMesh();
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
              for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                  res(rows(l), cols(m)) += Math::conj(bfi.integrate(m, l));
            }
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                    for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                      res(rows(l), cols(m)) += Math::conj(bfi.integrate(m, l));
                }
              }
            }
          }
        }
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Sequential<
    Math::SparseMatrix<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Variational::BilinearForm<
      Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>> final
    : public AssemblyBase<
        Math::SparseMatrix<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Variational::BilinearForm<
          Solution, TrialFES, TestFES,
          Math::SparseMatrix<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        std::vector<Eigen::Triplet<ScalarType>> triplets;
        Sequential<
          std::vector<Eigen::Triplet<ScalarType>>,
          Variational::BilinearForm<
            Solution, TrialFES, TestFES,
            std::vector<Eigen::Triplet<ScalarType>>>> assembly;
        assembly.execute(triplets, {
          input.getTrialFES(), input.getTestFES(),
          input.getLocalBFIs(), input.getGlobalBFIs() });
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class Solution, class TrialFES, class TestFES>
  class Sequential<
    std::vector<Eigen::Triplet<
      typename FormLanguage::Dot<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
    Variational::BilinearForm<Solution, TrialFES, TestFES,
      std::vector<Eigen::Triplet<
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>> final
    : public AssemblyBase<
        std::vector<Eigen::Triplet<
          typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>,
        Variational::BilinearForm<Solution, TrialFES, TestFES,
          std::vector<Eigen::Triplet<
            typename FormLanguage::Dot<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>>>
  {
    public:
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      using InputType = typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      void execute(OperatorType& res, const InputType& input) const override
      {
        const auto& mesh = input.getTrialFES().getMesh();
        res.reserve(input.getTestFES().getSize() * std::log(input.getTrialFES().getSize()));
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
              for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
              {
                for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                {
                  const ScalarType s = Math::conj(bfi.integrate(m, l));
                  if (s != ScalarType(0))
                    res.emplace_back(rows(l), cols(m), s);
                }
              }
            }
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              SequentialIteration trialseq(mesh, bfi.getTrialRegion());
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                  {
                    for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                    {
                      const ScalarType s = Math::conj(bfi.integrate(m, l));
                      if (s != ScalarType(0))
                        res.emplace_back(rows(l), cols(m), s);
                    }
                  }
                }
              }
            }
          }
        }
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class ... Solution, class ... TrialFES, class ... TestFES>
  class Sequential<
    std::vector<Eigen::Triplet<Real>>,
    Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>> final
      : public AssemblyBase<
          std::vector<Eigen::Triplet<Real>>,
          Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using TupleType =
        Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>...>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using Parent = AssemblyBase<OperatorType, TupleType>;

      using InputType = typename Parent::InputType;

      using Offsets = typename InputType::Offsets;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      void execute(OperatorType& res, const InputType& input) const override
      {
        using AssemblyTuple =
          Tuple<Sequential<std::vector<Eigen::Triplet<Real>>,
          Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>...>;

        AssemblyTuple assembly;

        const auto& t = input.getTuple();

        // Compute each block of triplets
        std::array<std::vector<Eigen::Triplet<Real>>, AssemblyTuple::Size> ts;
        assembly.zip(t).iapply(
            [&](const Index i, auto& p)
            {
              const auto& as = p.first();
              const auto& in = p.second();
              as.execute(ts[i], in);
            });

        // Add the triplets with the offsets
        size_t capacity = 0;
        for (const auto& v : ts)
          capacity += v.size();
        res.reserve(capacity);

        const Offsets& offsets = input.getOffsets();
        for (size_t i = 0; i < ts.size(); i++)
        {
          for (const Eigen::Triplet<Real>& t : ts[i])
          {
            res.emplace_back(
                t.row() + offsets[i].second(), t.col() + offsets[i].first(), t.value());
          }
        }
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class ... Solution, class ... TrialFES, class ... TestFES>
  class Sequential<
    Math::SparseMatrix<Real>,
    Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>> final
      : public AssemblyBase<
          Math::SparseMatrix<Real>,
          Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>>
    {
      public:
        using Parent =
          AssemblyBase<
            Math::SparseMatrix<Real>,
            Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>>;

        using InputType = typename Parent::InputType;

        using OperatorType = Math::SparseMatrix<Real>;

        Sequential() = default;

        Sequential(const Sequential& other)
          : Parent(other)
        {}

        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        void execute(OperatorType& res, const InputType& input) const override
        {
          Sequential<
            std::vector<Eigen::Triplet<Real>>,
            Tuple<
              Variational::BilinearForm<Solution, TrialFES, TestFES,
              std::vector<Eigen::Triplet<Real>>>...>> assembly;
          res.resize(input.getRows(), input.getColumns());
          std::vector<Eigen::Triplet<Real>> triplets;
          assembly.execute(triplets, input);
          res.setFromTriplets(triplets.begin(), triplets.end());
        }

        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  template <class ... FES>
  class Sequential<
    Math::Vector<Real>,
    Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>> final
      : public AssemblyBase<
          Math::Vector<Real>,
          Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>>
    {
      public:
        using Parent =
          AssemblyBase<
            Math::Vector<Real>,
            Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>>;

        using InputType = typename Parent::InputType;

        using VectorType = Math::Vector<Real>;

        Sequential() = default;

        Sequential(const Sequential& other)
          : Parent(other)
        {}

        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        void execute(VectorType& res, const InputType& input) const override
        {
          using AssemblyTuple =
            Tuple<
              Sequential<Math::Vector<Real>,
              Variational::LinearForm<FES, Math::Vector<Real>>>...>;

          const auto& t = input.getTuple();
          const auto& offsets = input.getOffsets();

          res.resize(input.getSize());
          res.setZero();

          AssemblyTuple assembly;
          VectorType vec;

          assembly.zip(t).iapply(
              [&](const Index i, const auto& p)
              {
                const auto& as = p.first();
                const auto& in = p.second();
                as.execute(vec, in);
                res.segment(offsets[i], vec.size()) = vec;
              });
        }

        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };


  template <class LinearSystem, class U1, class U2, class U3, class ... Us>
  class Sequential<
      LinearSystem,
      Variational::Problem<LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      using LinearSystemType = LinearSystem;

      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, U1, U2, U3, Us...>>;

      using InputType = typename Parent::InputType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        auto&       us           = input.getTrialFunctions();
        auto&       vs           = input.getTestFunctions();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets  = input.getTestOffsets();
        auto&       trialUUIDMap = input.getTrialUUIDMap();
        auto&       testUUIDMap  = input.getTestUUIDMap();

        const size_t ncols = input.getTotalTrialSize();
        const size_t nrows = input.getTotalTestSize();

        b.resize(nrows);
        b.setZero();

        std::vector<Eigen::Triplet<ScalarType>> triplets;

        const auto findTrialBlock = [&](const auto& uuid) -> size_t
        {
          auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return it->second;
        };

        const auto findTestBlock = [&](const auto& uuid) -> size_t
        {
          auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return it->second;
        };

        const auto& getTrialFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTrialBlock(uuid);
          const void* addr = nullptr;
          us.iapply([&](size_t i, const auto& uref)
          {
            if (i == k)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        const auto& getTestFESByUUID = [&](const auto& uuid) -> const auto&
        {
          const size_t k = findTestBlock(uuid);
          const void* addr = nullptr;
          vs.iapply([&](size_t i, const auto& vref)
          {
            if (i == k)
              addr = static_cast<const void*>(&vref.get().getFiniteElementSpace());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(vs.template get<0>().get().getFiniteElementSpace())>*>(addr);
        };

        const auto& mesh = [&]() -> const auto&
        {
          const void* addr = nullptr;
          us.apply([&](const auto& uref)
          {
            if (!addr)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace().getMesh());
          });
          assert(addr);
          return *static_cast<const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>*>(addr);
        }();

        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& uFES = getTrialFESByUUID(uUUID);
          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            const size_t d = it->getDimension();
            const Index  p = it->getIndex();

            bfi.setPolytope(*it);

            const auto& rows = vFES.getDOFs(d, p);
            const auto& cols = uFES.getDOFs(d, p);

            for (size_t i = 0; i < static_cast<size_t>(rows.size()); ++i)
            {
              const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows[i]));
              for (size_t j = 0; j < static_cast<size_t>(cols.size()); ++j)
              {
                const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols[j]));
                const ScalarType val = Math::conj(bfi.integrate(j, i));
                if (val != ScalarType(0))
                  triplets.emplace_back(I, J, val);
              }
            }
          }
        }

        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& uFES = getTrialFESByUUID(uUUID);
          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty() && !testAttrs.count(teIt->getAttribute()))
              continue;

            const auto& rows = vFES.getDOFs(teIt.getDimension(), teIt->getIndex());

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty() && !trialAttrs.count(trIt->getAttribute()))
                continue;

              const auto& cols = uFES.getDOFs(trIt.getDimension(), trIt->getIndex());

              bfi.setPolytope(*trIt, *teIt);

              for (size_t i = 0; i < static_cast<size_t>(rows.size()); ++i)
              {
                const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows[i]));
                for (size_t j = 0; j < static_cast<size_t>(cols.size()); ++j)
                {
                  const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols[j]));
                  const ScalarType val = Math::conj(bfi.integrate(j, i));
                  if (val != ScalarType(0))
                    triplets.emplace_back(I, J, val);
                }
              }
            }
          }
        }

        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& vFES = getTestFESByUUID(vUUID);

          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());

          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            lfi.setPolytope(*it);

            const auto& dofs = vFES.getDOFs(it.getDimension(), it->getIndex());
            for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
            {
              const Index I = static_cast<Index>(vOff + static_cast<size_t>(dofs[l]));
              b.coeffRef(I) -= lfi.integrate(l);
            }
          }
        }

        std::unordered_map<Index, ScalarType> fixed;
        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
            fixed.emplace(static_cast<Index>(uOff + local), static_cast<ScalarType>(value));
        }

        if (!fixed.empty())
        {
          std::vector<Eigen::Triplet<ScalarType>> filtered;
          filtered.reserve(triplets.size());

          for (const auto& t : triplets)
          {
            auto colIt = fixed.find(t.col());
            if (colIt != fixed.end() && t.row() != t.col())
              b.coeffRef(t.row()) -= t.value() * colIt->second;

            const bool rowFixed = fixed.find(t.row()) != fixed.end();
            const bool colFixed = colIt != fixed.end();

            if (rowFixed || colFixed)
              continue;
            filtered.emplace_back(t);
          }

          for (const auto& [idx, value] : fixed)
          {
            filtered.emplace_back(idx, idx, ScalarType(1));
            b.coeffRef(idx) = value;
          }

          triplets.swap(filtered);
        }

        if constexpr (std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>)
        {
          A.resize(nrows, ncols);
          A.setFromTriplets(triplets.begin(), triplets.end());
        }
        else
        {
          A.resize(nrows, ncols);
          A.setZero();
          std::sort(
            triplets.begin(), triplets.end(),
            [](const auto& a, const auto& b)
            {
              return a.row() == b.row() ? a.col() < b.col() : a.row() < b.row();
            });
          for (const auto& t : triplets)
            A(t.row(), t.col()) += t.value();
        }
      }

      Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Sequential<
    LinearSystem,
    Variational::Problem<LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      using LinearSystemType = LinearSystem;

      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, TrialFunction, TestFunction>>;

      using InputType = typename Parent::InputType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using BilinearFormType =
        Variational::BilinearForm<
          typename FormLanguage::Traits<TrialFunction>::SolutionType,
          typename FormLanguage::Traits<TrialFunction>::FESType,
          typename FormLanguage::Traits<TestFunction>::FESType,
          OperatorType>;

      using LinearFormType =
        Variational::LinearForm<
          typename FormLanguage::Traits<TestFunction>::FESType,
          VectorType>;

      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        auto& A = axb.getOperator();
        auto& b = axb.getVector();

        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();

        const auto& trialFES = u.getFiniteElementSpace();
        const auto& testFES  = v.getFiniteElementSpace();
        const auto& mesh     = trialFES.getMesh();

        const size_t cols = static_cast<size_t>(trialFES.getSize());
        const size_t rows = static_cast<size_t>(testFES.getSize());

        b.resize(rows);
        b.setZero();

        std::vector<Eigen::Triplet<ScalarType>> triplets;

        // Local BFIs
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            const size_t d = it->getDimension();
            const Index  p = it->getIndex();

            bfi.setPolytope(*it);

            const auto& rowsDOF = testFES.getDOFs(d, p);
            const auto& colsDOF = trialFES.getDOFs(d, p);

            for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
            {
              for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
              {
                const ScalarType val = Math::conj(bfi.integrate(j, i));
                if (val != ScalarType(0))
                  triplets.emplace_back(rowsDOF[i], colsDOF[j], val);
              }
            }
          }
        }

        // Global BFIs
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty() && !testAttrs.count(teIt->getAttribute()))
              continue;

            const auto& rowsDOF = testFES.getDOFs(teIt->getDimension(), teIt->getIndex());

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty() && !trialAttrs.count(trIt->getAttribute()))
                continue;

              const auto& colsDOF = trialFES.getDOFs(trIt->getDimension(), trIt->getIndex());

              bfi.setPolytope(*trIt, *teIt);

              for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
              {
                for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                {
                  const ScalarType val = Math::conj(bfi.integrate(j, i));
                  if (val != ScalarType(0))
                    triplets.emplace_back(rowsDOF[i], colsDOF[j], val);
                }
              }
            }
          }
        }

        // Preassembled bilinear forms
        for (auto& bf : pb.getBFs())
        {
          const auto& op = bf.getOperator();
          for (int k = 0; k < op.outerSize(); ++k)
            for (typename OperatorType::InnerIterator it(op, k); it; ++it)
              triplets.emplace_back(it.row(), it.col(), it.value());
        }

        // Linear forms
        for (auto& lfi : pb.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty() && !attrs.count(it->getAttribute()))
              continue;

            lfi.setPolytope(*it);
            const auto& dofs = testFES.getDOFs(it->getDimension(), it->getIndex());
            for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
              b.coeffRef(dofs[l]) -= lfi.integrate(l);
          }
        }

        // Preassembled linear forms
        for (auto& lf : pb.getLFs())
          b -= lf.getVector();

        // Dirichlet BC elimination in triplets
        std::unordered_map<Index, ScalarType> fixed;
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
            fixed.emplace(static_cast<Index>(local), static_cast<ScalarType>(value));
        }

        if (!fixed.empty())
        {
          for (const auto& t : triplets)
          {
            auto colIt = fixed.find(t.col());
            if (colIt != fixed.end() && t.row() != t.col())
              b.coeffRef(t.row()) -= t.value() * colIt->second;
          }

          std::vector<Eigen::Triplet<ScalarType>> filtered;
          filtered.reserve(triplets.size());

          for (const auto& t : triplets)
          {
            const bool rowFixed = fixed.find(t.row()) != fixed.end();
            const bool colFixed = fixed.find(t.col()) != fixed.end();

            if (rowFixed || colFixed)
              continue;
            filtered.emplace_back(t);
          }

          for (const auto& [idx, value] : fixed)
          {
            filtered.emplace_back(idx, idx, ScalarType(1));
            b.coeffRef(idx) = value;
          }

          triplets.swap(filtered);
        }

        A.resize(rows, cols);
        A.setFromTriplets(triplets.begin(), triplets.end());
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class Scalar, class Solution, class FES, class ValueDerived>
  class Sequential<
    IndexMap<Scalar>,
    Variational::DirichletBC<
      Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>> final
    : public AssemblyBase<
        IndexMap<Scalar>,
        Variational::DirichletBC<
          Variational::TrialFunction<Solution, FES>, Variational::FunctionBase<ValueDerived>>>
  {
    public:
      using FESType = FES;

      using TrialFunctionType = Variational::TrialFunction<Solution, FES>;

      using ValueType = Variational::FunctionBase<ValueDerived>;

      using DirichletBCType = Variational::DirichletBC<TrialFunctionType, ValueType>;

      using Parent = AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      using FESRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using InputType = typename Parent::InputType;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      void execute(IndexMap<Scalar>& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        const auto& essBdr = input.getEssentialBoundary();
        const auto& value = input.getValue();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t faceCount = mesh.getFaceCount();
        const size_t faceDim = mesh.getDimension() - 1;
        res.clear();
        for (Index i = 0; i < faceCount; i++)
        {
          if (mesh.isBoundary(i))
          {
            if (essBdr.size() == 0 || essBdr.count(mesh.getAttribute(faceDim, i)))
            {
              const auto& fe = fes.getFiniteElement(faceDim, i);
              const auto& mapping = fes.getPullback({ faceDim, i }, value);
              for (Index local = 0; local < fe.getCount(); local++)
              {
                const Index global = fes.getGlobalIndex({ faceDim, i }, local);
                auto find = res.find(global);
                if (find == res.end())
                  res.insert(find, std::pair{ global, fe.getLinearForm(local)(mapping) });
              }
            }
          }
        }
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif
