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

#include "Rodin/Math/Common.h"
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
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            lfi.setPolytope(*it);
            const size_t d = it.getDimension();
            const size_t i = it->getIndex();
            const auto& dofs = input.getFES().getDOFs(d, i);
            for (size_t l = 0; l < static_cast<size_t>(dofs.size()); l++)
              res(dofs(l)) += lfi.integrate(l);
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
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
            bfi.setPolytope(*it);
            const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
            const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
            for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
              for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                res(rows(l), cols(m)) += Math::conj(bfi.integrate(m, l));
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
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }
            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }
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
        res.clear();
        res.reserve(input.getTestFES().getSize() * std::log(input.getTrialFES().getSize()));
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }
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

        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }
            SequentialIteration trialseq(mesh, bfi.getTrialRegion());
            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }
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
        res.clear();
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

        constexpr bool IsSparse =
          std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>;

        if constexpr (!IsSparse)
        {
          A.resize(nrows, ncols);
          A.setZero();
        }

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

        // Mesh (assumed common across spaces)
        const auto& mesh = [&]() -> const auto&
        {
          const void* addr = nullptr;
          us.apply([&](const auto& uref)
          {
            if (!addr)
              addr = static_cast<const void*>(&uref.get().getFiniteElementSpace().getMesh());
          });
          assert(addr);
          return *static_cast<
            const std::decay_t<decltype(us.template get<0>().get().getFiniteElementSpace().getMesh())>*
          >(addr);
        }();

        // ------------------------------------------------------------
        // Dirichlet BC elimination data (NaN sentinel, global indexing)
        // ------------------------------------------------------------
        const size_t ndofs = std::max(nrows, ncols);
        std::vector<ScalarType> fixed(ndofs, Math::nan<ScalarType>());

        auto isFixed = [&](Index i) -> bool
        {
          const size_t k = static_cast<size_t>(i);
          return k < fixed.size() && !Math::isNaN(fixed[k]);
        };

        auto fixedValue = [&](Index i) -> ScalarType
        {
          return fixed[static_cast<size_t>(i)];
        };

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
          {
            const Index g = static_cast<Index>(uOff + static_cast<size_t>(local));
            if (static_cast<size_t>(g) < fixed.size())
              fixed[static_cast<size_t>(g)] = static_cast<ScalarType>(value);
          }
        }

        // ------------------------------------------------------------
        // Sparse-only: eliminate during assembly; Dense: accumulate
        // ------------------------------------------------------------
        auto sparse_entry = [&](Index row, Index col, ScalarType val)
        {
          if (val == ScalarType(0))
            return;

          if constexpr (IsSparse)
          {
            const bool rowFixed = isFixed(row);
            const bool colFixed = isFixed(col);

            if (rowFixed)
              return;

            if (colFixed && row != col)
            {
              b.coeffRef(row) -= val * fixedValue(col);
              return;
            }

            triplets.emplace_back(row, col, val);
          }
          else
          {
            A(row, col) += val;
          }
        };

        // ------------------------------------------------------------
        // Local BFIs (type-safe per-block FES access)
        // ------------------------------------------------------------
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());

          us.iapply([&](size_t ui, const auto& uref)
          {
            if (ui != uBlock) return;
            const auto& uFES = uref.get().getFiniteElementSpace();

            vs.iapply([&](size_t vi, const auto& vref)
            {
              if (vi != vBlock) return;
              const auto& vFES = vref.get().getFiniteElementSpace();

              for (auto it = seq.getIterator(); it; ++it)
              {
                if (!attrs.empty())
                {
                  const auto a = it->getAttribute();
                  if (!a || !attrs.count(*a))
                    continue;
                }

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
                    sparse_entry(I, J, val);
                  }
                }
              }
            });
          });
        }

        // ------------------------------------------------------------
        // Global BFIs (type-safe per-block FES access)
        // ------------------------------------------------------------
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto uUUID = bfi.getTrialFunction().getUUID();
          const auto vUUID = bfi.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();

          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          us.iapply([&](size_t ui, const auto& uref)
          {
            if (ui != uBlock) return;
            const auto& uFES = uref.get().getFiniteElementSpace();

            vs.iapply([&](size_t vi, const auto& vref)
            {
              if (vi != vBlock) return;
              const auto& vFES = vref.get().getFiniteElementSpace();

              for (auto teIt = testseq.getIterator(); teIt; ++teIt)
              {
                if (!testAttrs.empty())
                {
                  const auto a = teIt->getAttribute();
                  if (!a || !testAttrs.count(*a))
                    continue;
                }

                const auto& rows = vFES.getDOFs(teIt->getDimension(), teIt->getIndex());

                for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
                {
                  if (!trialAttrs.empty())
                  {
                    const auto a = trIt->getAttribute();
                    if (!a || !trialAttrs.count(*a))
                      continue;
                  }

                  const auto& cols = uFES.getDOFs(trIt->getDimension(), trIt->getIndex());

                  bfi.setPolytope(*trIt, *teIt);

                  for (size_t i = 0; i < static_cast<size_t>(rows.size()); ++i)
                  {
                    const Index I = static_cast<Index>(vOff + static_cast<size_t>(rows[i]));
                    for (size_t j = 0; j < static_cast<size_t>(cols.size()); ++j)
                    {
                      const Index J = static_cast<Index>(uOff + static_cast<size_t>(cols[j]));
                      const ScalarType val = Math::conj(bfi.integrate(j, i));
                      sparse_entry(I, J, val);
                    }
                  }
                }
              }
            });
          });
        }

        // ------------------------------------------------------------
        // LFIs (type-safe test FES access)
        // ------------------------------------------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto vUUID = lfi.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());

          vs.iapply([&](size_t vi, const auto& vref)
          {
            if (vi != vBlock) return;
            const auto& vFES = vref.get().getFiniteElementSpace();

            for (auto it = seq.getIterator(); it; ++it)
            {
              if (!attrs.empty())
              {
                const auto a = it->getAttribute();
                if (!a || !attrs.count(*a))
                  continue;
              }

              lfi.setPolytope(*it);

              const auto& dofs = vFES.getDOFs(it->getDimension(), it->getIndex());
              for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
              {
                const Index I = static_cast<Index>(vOff + static_cast<size_t>(dofs[l]));
                if constexpr (IsSparse)
                {
                  if (isFixed(I))
                    continue;
                }
                b.coeffRef(I) -= lfi.integrate(l);
              }
            }
          });
        }

        // ------------------------------------------------------------
        // Finalize
        // ------------------------------------------------------------
        if constexpr (IsSparse)
        {
          for (Index i = 0; i < static_cast<Index>(nrows); ++i)
          {
            if (isFixed(i))
            {
              triplets.emplace_back(i, i, ScalarType(1));
              b.coeffRef(i) = fixedValue(i);
            }
          }

          A.resize(nrows, ncols);
          A.setFromTriplets(triplets.begin(), triplets.end());
        }
        else
        {
          for (Index idx = 0; idx < static_cast<Index>(nrows); ++idx)
          {
            if (!isFixed(idx))
              continue;

            const ScalarType value = fixedValue(idx);

            for (size_t r = 0; r < nrows; ++r)
            {
              if (r == static_cast<size_t>(idx))
                continue;
              b.coeffRef(r) -= A(r, idx) * value;
              A(r, idx) = ScalarType(0);
            }

            for (size_t c = 0; c < ncols; ++c)
            {
              if (c == static_cast<size_t>(idx))
                continue;
              A(idx, c) = ScalarType(0);
            }

            A(idx, idx) = ScalarType(1);
            b.coeffRef(static_cast<size_t>(idx)) = value;
          }
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

        constexpr bool IsSparse =
          std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>;

        // ------------------------------------------------------------
        // Dirichlet BC elimination data (NaN sentinel)
        // ------------------------------------------------------------
        std::vector<ScalarType> fixed(rows, Math::nan<ScalarType>());

        auto isFixed = [&](Index i) -> bool
        {
          return !Math::isNaN(fixed[static_cast<size_t>(i)]);
        };

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          for (const auto& [local, value] : dofs)
            fixed[static_cast<size_t>(local)] = static_cast<ScalarType>(value);
        }

        // ------------------------------------------------------------
        // Matrix init
        // ------------------------------------------------------------
        std::vector<Eigen::Triplet<ScalarType>> triplets;
        if constexpr (!IsSparse)
        {
          A.resize(rows, cols);
          A.setZero();
        }

        // ------------------------------------------------------------
        // Sparse-only: eliminate during assembly
        // Dense: do plain accumulation; eliminate afterwards
        // ------------------------------------------------------------
        auto sparse_entry = [&](Index row, Index col, ScalarType val)
        {
          if (val == ScalarType(0))
            return;

          if constexpr (IsSparse)
          {
            const bool rowFixed = isFixed(row);
            const bool colFixed = isFixed(col);

            if (rowFixed)
              return;

            if (colFixed && row != col)
            {
              b.coeffRef(row) -= val * fixed[static_cast<size_t>(col)];
              return;
            }

            triplets.emplace_back(row, col, val);
          }
          else
          {
            A(row, col) += val;
          }
        };

        // ------------------------------------------------------------
        // Local BFIs
        // ------------------------------------------------------------
        for (auto& bfi : pb.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

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
                sparse_entry(rowsDOF[i], colsDOF[j], val);
              }
            }
          }
        }

        // ------------------------------------------------------------
        // Global BFIs
        // ------------------------------------------------------------
        for (auto& bfi : pb.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs  = bfi.getTestAttributes();
          SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          SequentialIteration testseq(mesh, bfi.getTestRegion());

          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (!testAttrs.empty())
            {
              const auto a = teIt->getAttribute();
              if (!a || !testAttrs.count(*a))
                continue;
            }

            const auto& rowsDOF = testFES.getDOFs(teIt->getDimension(), teIt->getIndex());

            for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
            {
              if (!trialAttrs.empty())
              {
                const auto a = trIt->getAttribute();
                if (!a || !trialAttrs.count(*a))
                  continue;
              }

              const auto& colsDOF = trialFES.getDOFs(trIt->getDimension(), trIt->getIndex());

              bfi.setPolytope(*trIt, *teIt);

              for (size_t i = 0; i < static_cast<size_t>(rowsDOF.size()); ++i)
              {
                for (size_t j = 0; j < static_cast<size_t>(colsDOF.size()); ++j)
                {
                  const ScalarType val = Math::conj(bfi.integrate(j, i));
                  sparse_entry(rowsDOF[i], colsDOF[j], val);
                }
              }
            }
          }
        }

        // ------------------------------------------------------------
        // Preassembled bilinear forms
        // ------------------------------------------------------------
        for (auto& bf : pb.getBFs())
        {
          const auto& op = bf.getOperator();
          if constexpr (IsSparse)
          {
            for (int k = 0; k < op.outerSize(); ++k)
              for (typename OperatorType::InnerIterator it(op, k); it; ++it)
                sparse_entry(it.row(), it.col(), it.value());
          }
          else
          {
            A += op;
          }
        }

        // ------------------------------------------------------------
        // Linear forms (unchanged)
        // ------------------------------------------------------------
        for (auto& lfi : pb.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (!attrs.empty())
            {
              const auto a = it->getAttribute();
              if (!a || !attrs.count(*a))
                continue;
            }

            lfi.setPolytope(*it);
            const auto& dofs = testFES.getDOFs(it->getDimension(), it->getIndex());
            for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
              b.coeffRef(dofs[l]) -= lfi.integrate(l);
          }
        }

        for (auto& lf : pb.getLFs())
          b -= lf.getVector();

        // ------------------------------------------------------------
        // Finalize: Sparse build; Dense eliminate afterwards
        // ------------------------------------------------------------
        if constexpr (IsSparse)
        {
          // Inject identity rows for fixed dofs
          for (Index i = 0; i < static_cast<Index>(rows); ++i)
          {
            if (isFixed(i))
            {
              triplets.emplace_back(i, i, ScalarType(1));
              b.coeffRef(i) = fixed[static_cast<size_t>(i)];
            }
          }

          A.resize(rows, cols);
          A.setFromTriplets(triplets.begin(), triplets.end());
        }
        else
        {
          // Dense elimination after full assembly (your original approach, driven by fixedValue)
          for (Index idx = 0; idx < static_cast<Index>(rows); ++idx)
          {
            if (!isFixed(idx))
              continue;

            const ScalarType value = fixed[static_cast<size_t>(idx)];

            for (size_t r = 0; r < rows; ++r)
            {
              if (r == static_cast<size_t>(idx))
                continue;
              b.coeffRef(r) -= A(r, idx) * value;
              A(r, idx) = ScalarType(0);
            }

            for (size_t c = 0; c < cols; ++c)
            {
              if (c == static_cast<size_t>(idx))
                continue;
              A(idx, c) = ScalarType(0);
            }

            A(idx, idx) = ScalarType(1);
            b.coeffRef(static_cast<size_t>(idx)) = value;
          }
        }
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
            if (!essBdr.empty())
            {
              const auto a = mesh.getAttribute(faceDim, i);
              if (!a || !essBdr.count(*a))
                continue;
            }

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

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif
