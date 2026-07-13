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
#include "Rodin/Variational/IntegrationPoint.h"

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"

#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/ConstraintMap.h"

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
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Vector type of the linear system.
      using VectorType = Math::Vector<ScalarType>;

      /// @brief Linear form type assembled by this backend.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
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
      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Dot<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      /// @brief Assembled operator type.
      using OperatorType = Math::Matrix<ScalarType>;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Bilinear form type assembled by this backend.
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
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
      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      /// @brief Assembled operator type.
      using OperatorType = Math::SparseMatrix<ScalarType>;

      /// @brief Bilinear form type assembled by this backend.
      using BilinearFormType = Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential bilinear form assembly into Eigen triplets.
   */
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
      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Dot<
          typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type;

      /// @brief Assembled operator type.
      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      /// @brief Bilinear form type assembled by this backend.
      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType =
        Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType =
        Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential block bilinear form assembly into Eigen triplets.
   */
  template <class ... Solution, class ... TrialFES, class ... TestFES>
  class Sequential<
    std::vector<Eigen::Triplet<Real>>,
    Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>> final
      : public AssemblyBase<
          std::vector<Eigen::Triplet<Real>>,
          Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Assembled operator type.
      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      /// @brief Tuple of bilinear forms assembled by this backend.
      using TupleType =
        Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>...>;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Global bilinear form integrator base type.
      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<OperatorType, TupleType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Block offset array type.
      using Offsets = typename InputType::Offsets;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes block triplet assembly.
       * @param res Output triplet array.
       * @param input Block assembly input.
       */
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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential block bilinear form assembly into a sparse matrix.
   */
  template <class ... Solution, class ... TrialFES, class ... TestFES>
  class Sequential<
    Math::SparseMatrix<Real>,
    Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>> final
      : public AssemblyBase<
          Math::SparseMatrix<Real>,
          Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>>
    {
      public:
        /// @brief Parent class type.
        using Parent =
          AssemblyBase<
            Math::SparseMatrix<Real>,
            Tuple<Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>...>>;

        /// @brief Assembly input data type.
        using InputType = typename Parent::InputType;

        /// @brief Assembled operator type.
        using OperatorType = Math::SparseMatrix<Real>;

        /// @brief Default constructor.
        Sequential() = default;

        /// @brief Copy constructor.
        Sequential(const Sequential& other)
          : Parent(other)
        {}

        /// @brief Move constructor.
        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        /**
         * @brief Executes block sparse matrix assembly.
         * @param res Output sparse matrix.
         * @param input Block assembly input.
         */
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

        /**
         * @brief Creates a polymorphic copy.
         * @return Pointer to a new copy.
         */
        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  /**
   * @brief Sequential block linear form assembly into a vector.
   */
  template <class ... FES>
  class Sequential<
    Math::Vector<Real>,
    Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>> final
      : public AssemblyBase<
          Math::Vector<Real>,
          Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>>
    {
      public:
        /// @brief Parent class type.
        using Parent =
          AssemblyBase<
            Math::Vector<Real>,
            Tuple<Variational::LinearForm<FES, Math::Vector<Real>>...>>;

        /// @brief Assembly input data type.
        using InputType = typename Parent::InputType;

        /// @brief Vector type of the linear system.
        using VectorType = Math::Vector<Real>;

        /// @brief Default constructor.
        Sequential() = default;

        /// @brief Copy constructor.
        Sequential(const Sequential& other)
          : Parent(other)
        {}

        /// @brief Move constructor.
        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        /**
         * @brief Executes block vector assembly.
         * @param res Output vector.
         * @param input Block assembly input.
         */
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

        /**
         * @brief Creates a polymorphic copy.
         * @return Pointer to a new copy.
         */
        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };


  /**
   * @brief Sequential mixed problem assembly.
   */
  template <class LinearSystem, class U1, class U2, class U3, class ... Us>
  class Sequential<
      LinearSystem,
      Variational::Problem<LinearSystem, U1, U2, U3, Us...>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = LinearSystem;

      /// @brief Parent class type.
      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, U1, U2, U3, Us...>>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Assembled operator type.
      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      /// @brief Vector type of the linear system.
      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      /**
       * @brief Executes full mixed problem assembly.
       * @param axb Output linear system.
       * @param input Mixed problem input.
       */
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
        // Dirichlet BC elimination + identification data
        // ------------------------------------------------------------
        ConstraintMap<ScalarType> constraints(std::max(nrows, ncols));

        using DBCBaseType =
          Variational::DirichletBCBase<ScalarType>;
        using ValueDOFsType    = typename DBCBaseType::ValueDOFs;
        using IdentDOFsType    = typename DBCBaseType::IdentifiedDOFs;

        for (auto& dbc : pb.getDBCs())
        {
          const auto uUUID = dbc.getOperand().getUUID();
          const size_t uBlock = findTrialBlock(uUUID);
          const size_t uOff   = trialOffsets[uBlock];

          dbc.assemble();
          std::visit([&](auto&& dofs)
          {
            using T = std::decay_t<decltype(dofs)>;
            if constexpr (std::is_same_v<T, ValueDOFsType>)
            {
              for (const auto& [local, value] : dofs)
              {
                const Index g = static_cast<Index>(
                    uOff + static_cast<size_t>(local));
                constraints.setFixed(g, static_cast<ScalarType>(value));
              }
            }
            else if constexpr (std::is_same_v<T, IdentDOFsType>)
            {
              const auto vUUIDOpt = dbc.getValueUUID();
              assert(vUUIDOpt
                  && "Identification DBC missing value UUID");
              const size_t vBlock = findTrialBlock(*vUUIDOpt);
              const size_t vOff   = trialOffsets[vBlock];
              const auto& affineValues = dbc.getIdentificationValues();
              for (const auto& [slave, pair] : dofs)
              {
                const auto& masters = pair.first;
                const auto& coeffs  = pair.second;
                const Index gs = static_cast<Index>(
                    uOff + static_cast<size_t>(slave));
                const Index n =
                    static_cast<Index>(masters.size());
                std::vector<typename ConstraintMap<ScalarType>::Entry> entries;
                entries.reserve(static_cast<size_t>(n));
                for (Index k = 0; k < n; k++)
                {
                  const Index gm = static_cast<Index>(
                      vOff + static_cast<size_t>(masters[k]));
                  entries.push_back(
                      { gm, static_cast<ScalarType>(coeffs[k]) });
                }
                const auto valueIt = affineValues.find(slave);
                const ScalarType value =
                  valueIt == affineValues.end()
                    ? ScalarType(0)
                    : static_cast<ScalarType>(valueIt->second);
                constraints.setIdentification(gs, entries, value);
              }
            }
          }, dbc.getDOFs());
        }

        // ------------------------------------------------------------
        // Sparse-only: eliminate during assembly; Dense: accumulate
        // ------------------------------------------------------------
        auto matrixEntry = [&](Index row, Index col, ScalarType val) {
          if (val == ScalarType(0))
            return;

          const ScalarType colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : ScalarType(0);

          if constexpr (IsSparse)
          {
            for (const auto& r : constraints.expand(row))
            {
              if (colValue != ScalarType(0))
                b.coeffRef(r.index) -= r.coefficient * val * colValue;
              for (const auto& c : constraints.expand(col))
                triplets.emplace_back(
                    r.index, c.index, r.coefficient * val * c.coefficient);
            }
          }
          else
          {
            for (const auto& r : constraints.expand(row))
            {
              if (colValue != ScalarType(0))
                b.coeffRef(r.index) -= r.coefficient * val * colValue;
              for (const auto& c : constraints.expand(col))
                A(r.index, c.index) += r.coefficient * val * c.coefficient;
            }
          }
        };

        auto vectorEntry = [&](Index row, ScalarType val) {
          if (val == ScalarType(0))
            return;

          for (const auto& r : constraints.expand(row))
            b.coeffRef(r.index) += r.coefficient * val;
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
                    matrixEntry(I, J, val);
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
                      matrixEntry(I, J, val);
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
                vectorEntry(I, -static_cast<ScalarType>(lfi.integrate(l)));
              }
            }
          });
        }

        // ------------------------------------------------------------
        // Preassembled bilinear forms (with block offsets)
        // ------------------------------------------------------------
        for (auto& bf : pb.getBFs())
        {
          const auto uUUID = bf.getTrialFunction().getUUID();
          const auto vUUID = bf.getTestFunction().getUUID();

          const size_t uBlock = findTrialBlock(uUUID);
          const size_t vBlock = findTestBlock(vUUID);

          const size_t uOff = trialOffsets[uBlock];
          const size_t vOff = testOffsets[vBlock];

          const auto& op = bf.getOperator();
          if constexpr (IsSparse)
          {
            for (int k = 0; k < op.outerSize(); ++k)
              for (typename OperatorType::InnerIterator it(op, k); it; ++it)
                matrixEntry(static_cast<Index>(vOff) + it.row(),
                  static_cast<Index>(uOff) + it.col(), it.value());
          }
          else
          {
            const auto opRows = op.rows();
            const auto opCols = op.cols();
            for (Eigen::Index i = 0; i < opRows; ++i)
              for (Eigen::Index j = 0; j < opCols; ++j)
              {
                const auto val = op(i, j);
                if (val != ScalarType(0))
                  matrixEntry(
                    static_cast<Index>(vOff) + i, static_cast<Index>(uOff) + j, val);
              }
          }
        }

        // ------------------------------------------------------------
        // Preassembled linear forms (with block offsets)
        // ------------------------------------------------------------
        for (auto& lf : pb.getLFs())
        {
          const auto vUUID = lf.getTestFunction().getUUID();
          const size_t vBlock = findTestBlock(vUUID);
          const size_t vOff   = testOffsets[vBlock];

          const auto& vec = lf.getVector();
          for (Eigen::Index i = 0; i < vec.size(); ++i)
            vectorEntry(
              static_cast<Index>(vOff) + i, static_cast<ScalarType>(vec.coeff(i)));
        }

        // ------------------------------------------------------------
        // Finalize
        // ------------------------------------------------------------
        if constexpr (IsSparse)
        {
          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (static_cast<size_t>(gs) >= nrows)
              continue;
            triplets.emplace_back(gs, gs, ScalarType(1));
            for (const auto& e : constraints.expand(gs))
              triplets.emplace_back(gs, e.index, -e.coefficient);
            b.coeffRef(gs) = constraints.getIdentificationValue(gs);
          }

          std::vector<Eigen::Triplet<ScalarType>> filteredTriplets;
          filteredTriplets.reserve(triplets.size() + nrows);
          for (const auto& t : triplets)
          {
            const Index row = static_cast<Index>(t.row());
            const Index col = static_cast<Index>(t.col());
            const ScalarType val = t.value();
            if (constraints.isFixed(row))
              continue;
            if (constraints.isFixed(col) && row != col)
            {
              b.coeffRef(row) -= val * constraints.getFixedValue(col);
              continue;
            }
            filteredTriplets.push_back(t);
          }

          for (Index i = 0; i < static_cast<Index>(nrows); ++i)
          {
            if (constraints.isFixed(i))
            {
              filteredTriplets.emplace_back(i, i, ScalarType(1));
              b.coeffRef(i) = constraints.getFixedValue(i);
            }
          }

          A.resize(nrows, ncols);
          A.setFromTriplets(filteredTriplets.begin(), filteredTriplets.end());
        }
        else
        {
          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (static_cast<size_t>(gs) >= nrows)
              continue;
            for (Index c = 0; c < static_cast<Index>(ncols); ++c)
              A(gs, c) = ScalarType(0);
            A(gs, gs) = ScalarType(1);
            for (const auto& e : constraints.expand(gs))
              A(gs, e.index) -= e.coefficient;
            b.coeffRef(gs) = constraints.getIdentificationValue(gs);
          }

          for (Index idx = 0; idx < static_cast<Index>(nrows); ++idx)
          {
            if (!constraints.isFixed(idx))
              continue;

            const ScalarType value = constraints.getFixedValue(idx);

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

      // Targeted (LHS-only / RHS-only) assembly for the block Eigen backend:
      // assemble the full system into a scratch object and expose only the
      // requested side, leaving the other operand untouched (the targeted
      // contract). Keeps the block BC-elimination logic in one code path.
      /**
       * @brief Executes targeted mixed problem assembly.
       * @param axb Output linear system.
       * @param input Mixed problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
          LinearSystemType& axb,
          const InputType& input,
          Rodin::Variational::AssemblyTarget target) const
      {
        LinearSystemType scratch;
        execute(scratch, input);
        if (target == Rodin::Variational::AssemblyTarget::LHS)
          axb.getOperator() = std::move(scratch.getOperator());
        else
          axb.getVector() = std::move(scratch.getVector());
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  /**
   * @brief Sequential single-field problem assembly.
   */
  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Sequential<
    LinearSystem,
    Variational::Problem<LinearSystem, TrialFunction, TestFunction>> final
    : public AssemblyBase<
        LinearSystem,
        Variational::Problem<LinearSystem, TrialFunction, TestFunction>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = LinearSystem;

      /// @brief Parent class type.
      using Parent =
        AssemblyBase<
          LinearSystemType,
          Variational::Problem<LinearSystemType, TrialFunction, TestFunction>>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Assembled operator type.
      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      /// @brief Vector type of the linear system.
      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      /// @brief Scalar value type.
      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      /// @brief Bilinear form type assembled into the operator.
      using BilinearFormType =
        Variational::BilinearForm<
          typename FormLanguage::Traits<TrialFunction>::SolutionType,
          typename FormLanguage::Traits<TrialFunction>::FESType,
          typename FormLanguage::Traits<TestFunction>::FESType,
          OperatorType>;

      /// @brief Linear form type assembled into the vector.
      using LinearFormType =
        Variational::LinearForm<
          typename FormLanguage::Traits<TestFunction>::FESType,
          VectorType>;

      /**
       * @brief Executes full problem assembly.
       * @param axb Output linear system.
       * @param input Problem assembly input.
       */
      void execute(LinearSystemType& axb, const InputType& input) const override
      {
        execute(axb, input, AssemblyMode::Full);
      }

      /**
       * @brief Executes targeted problem assembly.
       * @param axb Output linear system.
       * @param input Problem assembly input.
       * @param target Side of the system to assemble.
       */
      void execute(
          LinearSystemType& axb,
          const InputType& input,
          Rodin::Variational::AssemblyTarget target) const
      {
        switch (target)
        {
          case Rodin::Variational::AssemblyTarget::LHS:
            execute(axb, input, AssemblyMode::LHS);
            break;
          case Rodin::Variational::AssemblyTarget::RHS:
            execute(axb, input, AssemblyMode::RHS);
            break;
        }
      }

    private:
      enum class AssemblyMode
      {
        Full,
        LHS,
        RHS
      };

      void execute(
          LinearSystemType& axb,
          const InputType& input,
          AssemblyMode mode) const
      {
        const bool doMatrix = mode != AssemblyMode::RHS;
        const bool doVector = mode != AssemblyMode::LHS;

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

        if (doVector)
        {
          b.resize(rows);
          b.setZero();
        }

        constexpr bool IsSparse =
          std::is_base_of_v<Eigen::SparseMatrixBase<OperatorType>, OperatorType>;

        // ------------------------------------------------------------
        // Dirichlet BC elimination + identification data
        // ------------------------------------------------------------
        ConstraintMap<ScalarType> constraints(std::max(rows, cols));

        using DBCBaseType =
          Variational::DirichletBCBase<ScalarType>;
        using ValueDOFsType = typename DBCBaseType::ValueDOFs;
        using IdentDOFsType = typename DBCBaseType::IdentifiedDOFs;

        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() != u.getUUID())
            continue;

          dbc.assemble();
          std::visit([&](auto&& dofs)
          {
            using T = std::decay_t<decltype(dofs)>;
            if constexpr (std::is_same_v<T, ValueDOFsType>)
            {
              for (const auto& [local, value] : dofs)
                constraints.setFixed(
                    static_cast<Index>(local),
                    static_cast<ScalarType>(value));
            }
            else if constexpr (std::is_same_v<T, IdentDOFsType>)
            {
              // Single-FES problem: master must live in the same block as
              // the slave (i.e. v.getLeaf() and u share the same FES); we
              // therefore use no offset.
              const auto& affineValues = dbc.getIdentificationValues();
              for (const auto& [slave, pair] : dofs)
              {
                const auto& masters = pair.first;
                const auto& coeffs  = pair.second;
                const Index n =
                    static_cast<Index>(masters.size());
                std::vector<typename ConstraintMap<ScalarType>::Entry> entries;
                entries.reserve(static_cast<size_t>(n));
                for (Index k = 0; k < n; k++)
                {
                  entries.push_back({
                      static_cast<Index>(masters[k]),
                      static_cast<ScalarType>(coeffs[k]) });
                }
                const auto valueIt = affineValues.find(slave);
                const ScalarType value =
                  valueIt == affineValues.end()
                    ? ScalarType(0)
                    : static_cast<ScalarType>(valueIt->second);
                constraints.setIdentification(static_cast<Index>(slave), entries, value);
              }
            }
          }, dbc.getDOFs());
        }

        if (mode != AssemblyMode::Full && !constraints.getIdentifiedRows().empty())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Targeted assembly is not implemented for identification DirichletBCs."
            << Alert::Raise;
        }

        // ------------------------------------------------------------
        // Matrix init
        // ------------------------------------------------------------
        std::vector<Eigen::Triplet<ScalarType>> triplets;
        if (doMatrix)
        {
          if constexpr (!IsSparse)
          {
            A.resize(rows, cols);
            A.setZero();
          }
        }

        // ------------------------------------------------------------
        // Sparse-only: eliminate during assembly
        // Dense: do plain accumulation; eliminate afterwards
        // ------------------------------------------------------------
        auto matrixEntry = [&](Index row, Index col, ScalarType val) {
          if (val == ScalarType(0))
            return;

          const ScalarType colValue =
            constraints.isIdentified(col)
              ? constraints.getIdentificationValue(col)
              : ScalarType(0);

          if constexpr (IsSparse)
          {
            for (const auto& r : constraints.expand(row))
            {
              if (colValue != ScalarType(0))
                b.coeffRef(r.index) -= r.coefficient * val * colValue;
              for (const auto& c : constraints.expand(col))
                triplets.emplace_back(
                    r.index, c.index, r.coefficient * val * c.coefficient);
            }
          }
          else
          {
            for (const auto& r : constraints.expand(row))
            {
              if (colValue != ScalarType(0))
                b.coeffRef(r.index) -= r.coefficient * val * colValue;
              for (const auto& c : constraints.expand(col))
                A(r.index, c.index) += r.coefficient * val * c.coefficient;
            }
          }
        };

        auto vectorEntry = [&](Index row, ScalarType val) {
          if (val == ScalarType(0))
            return;

          for (const auto& r : constraints.expand(row))
            b.coeffRef(r.index) += r.coefficient * val;
        };

        // ------------------------------------------------------------
        // Local BFIs
        // ------------------------------------------------------------
        if (doMatrix)
        {
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
                matrixEntry(rowsDOF[i], colsDOF[j], val);
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
                  matrixEntry(rowsDOF[i], colsDOF[j], val);
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
                matrixEntry(it.row(), it.col(), it.value());
          }
          else
          {
            const auto opRows = op.rows();
            const auto opCols = op.cols();
            for (Eigen::Index i = 0; i < opRows; ++i)
              for (Eigen::Index j = 0; j < opCols; ++j)
              {
                const auto val = op(i, j);
                if (val != ScalarType(0))
                  matrixEntry(static_cast<Index>(i), static_cast<Index>(j), val);
              }
          }
        }
        } // doMatrix

        // ------------------------------------------------------------
        // Linear forms (unchanged)
        // ------------------------------------------------------------
        if (doVector)
        {
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
              vectorEntry(dofs[l], -static_cast<ScalarType>(lfi.integrate(l)));
          }
        }

        for (auto& lf : pb.getLFs())
        {
          const auto& vec = lf.getVector();
          for (Eigen::Index i = 0; i < vec.size(); ++i)
            vectorEntry(static_cast<Index>(i), static_cast<ScalarType>(vec.coeff(i)));
        }
        } // doVector

        // ------------------------------------------------------------
        // Finalize (Full): Sparse build; Dense eliminate afterwards
        // ------------------------------------------------------------
        if (mode == AssemblyMode::Full)
        {
        if constexpr (IsSparse)
        {
          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (static_cast<size_t>(gs) >= rows)
              continue;
            triplets.emplace_back(gs, gs, ScalarType(1));
            for (const auto& e : constraints.expand(gs))
              triplets.emplace_back(gs, e.index, -e.coefficient);
            b.coeffRef(gs) = constraints.getIdentificationValue(gs);
          }

          std::vector<Eigen::Triplet<ScalarType>> filteredTriplets;
          filteredTriplets.reserve(triplets.size() + rows);
          for (const auto& t : triplets)
          {
            const Index row = static_cast<Index>(t.row());
            const Index col = static_cast<Index>(t.col());
            const ScalarType val = t.value();
            if (constraints.isFixed(row))
              continue;
            if (constraints.isFixed(col) && row != col)
            {
              b.coeffRef(row) -= val * constraints.getFixedValue(col);
              continue;
            }
            filteredTriplets.push_back(t);
          }

          for (Index i = 0; i < static_cast<Index>(rows); ++i)
          {
            if (constraints.isFixed(i))
            {
              filteredTriplets.emplace_back(i, i, ScalarType(1));
              b.coeffRef(i) = constraints.getFixedValue(i);
            }
          }

          A.resize(rows, cols);
          A.setFromTriplets(filteredTriplets.begin(), filteredTriplets.end());
        }
        else
        {
          for (const Index gs : constraints.getIdentifiedRows())
          {
            if (static_cast<size_t>(gs) >= rows)
              continue;
            for (size_t c = 0; c < cols; ++c)
              A(gs, c) = ScalarType(0);
            A(gs, gs) = ScalarType(1);
            for (const auto& e : constraints.expand(gs))
              A(gs, e.index) -= e.coefficient;
            b.coeffRef(gs) = constraints.getIdentificationValue(gs);
          }

          for (Index idx = 0; idx < static_cast<Index>(rows); ++idx)
          {
            if (!constraints.isFixed(idx))
              continue;

            const ScalarType value = constraints.getFixedValue(idx);

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
        } // mode == Full
        else if (doMatrix)
        {
          // LHS-only: assemble the operator; fixed rows -> identity (rows only,
          // columns kept), mirroring the PETSc MatZeroRows targeted path.
          if constexpr (IsSparse)
          {
            std::vector<Eigen::Triplet<ScalarType>> filteredTriplets;
            filteredTriplets.reserve(triplets.size() + rows);
            for (const auto& t : triplets)
            {
              if (constraints.isFixed(static_cast<Index>(t.row())))
                continue;
              filteredTriplets.push_back(t);
            }
            for (Index i = 0; i < static_cast<Index>(rows); ++i)
              if (constraints.isFixed(i))
                filteredTriplets.emplace_back(i, i, ScalarType(1));
            A.resize(rows, cols);
            A.setFromTriplets(filteredTriplets.begin(), filteredTriplets.end());
          }
          else
          {
            for (Index idx = 0; idx < static_cast<Index>(rows); ++idx)
            {
              if (!constraints.isFixed(idx))
                continue;
              for (size_t c = 0; c < cols; ++c)
                A(idx, c) = ScalarType(0);
              A(idx, idx) = ScalarType(1);
            }
          }
        }
        else
        {
          // RHS-only: assemble the vector; fixed entries -> prescribed value.
          for (Index idx = 0; idx < static_cast<Index>(rows); ++idx)
            if (constraints.isFixed(idx))
              b.coeffRef(static_cast<size_t>(idx)) = constraints.getFixedValue(idx);
        }
      }

    public:
      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential value Dirichlet boundary condition assembly.
   */
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
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Trial function type constrained by the boundary condition.
      using TrialFunctionType = Variational::TrialFunction<Solution, FES>;

      /// @brief Boundary value function type.
      using ValueType = Variational::FunctionBase<ValueDerived>;

      /// @brief Dirichlet condition type.
      using DirichletBCType = Variational::DirichletBC<TrialFunctionType, ValueType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<IndexMap<Scalar>, DirichletBCType>;

      /// @brief Range type of the finite element space.
      using FESRangeType = typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes value Dirichlet boundary condition assembly.
       * @param res Output map from constrained DOFs to values.
       * @param input Boundary condition input.
       */
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
          if (essBdr.empty() && !mesh.isBoundary(i))
            continue;

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

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembler for the identification Dirichlet BC
   *        `u = A(v)`.
   *
   * Iterates the default exterior boundary when no attributes are specified.
   * When attributes are specified, every tagged codimension-one face is
   * eligible, including interior interface faces. On each selected face, the
   * slave DOFs (from
   * @f$ u @f$'s FES) and the master DOFs (from @f$ v @f$'s FES) are obtained
   * exactly via @c FiniteElementSpaceBase::getDOFs(face) — no geometric
   * matching, no tolerance.
   *
   * For each slave DOF a pair @c (masters, coefficients) is emitted. Each
   * coefficient is computed as @f$ \ell_s^u(A(\phi_j^v)) @f$: the master
   * basis is evaluated through @f$ A(v) @f$ and then sampled by the slave
   * finite element's DOF functional.
   */
  template <class Scalar, class Sol1, class FES1,
            class Derived2, class FES2,
            Variational::ShapeFunctionSpaceType Sp>
  class Sequential<
    IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
    Variational::DirichletBC<
      Variational::TrialFunction<Sol1, FES1>,
      Variational::ShapeFunctionBase<Derived2, FES2, Sp>>> final
    : public AssemblyBase<
        IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>,
        Variational::DirichletBC<
          Variational::TrialFunction<Sol1, FES1>,
          Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
  {
    public:
      /// @brief Output map type for slave-to-master identifications.
      using OutputType = IndexMap<std::pair<IndexArray, Math::Vector<Scalar>>>;

      /// @brief Slave trial function type.
      using TrialFunctionType = Variational::TrialFunction<Sol1, FES1>;

      /// @brief Shape-function expression type on the right-hand side.
      using ValueType = Variational::ShapeFunctionBase<Derived2, FES2, Sp>;

      /// @brief Dirichlet condition type.
      using DirichletBCType =
        Variational::DirichletBC<TrialFunctionType, ValueType>;

      /// @brief Parent class type.
      using Parent = AssemblyBase<OutputType, DirichletBCType>;

      /// @brief Assembly input data type.
      using InputType = typename Parent::InputType;

      /// @brief Default constructor.
      Sequential() = default;

      /// @brief Copy constructor.
      Sequential(const Sequential& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes identification Dirichlet boundary condition assembly.
       * @param res Output map from slave DOFs to master DOF coefficients.
       * @param input Boundary condition input.
       */
      void execute(OutputType& res, const InputType& input) const override
      {
        const auto& u = input.getOperand();
        // setIntegrationPoint mutates internal evaluation state; the input
        // exposes Av as const, but we need to drive its IP cursor while
        // probing each basis. The mutation is purely evaluation state, not
        // semantic.
        auto& Av = const_cast<ValueType&>(input.getShapeFunction());
        const auto& essBdr = input.getEssentialBoundary();

        const auto& fesU = u.getFiniteElementSpace();
        const auto& fesV = Av.getLeaf().getFiniteElementSpace();
        const auto& mesh = fesU.getMesh();
        const size_t faceCount = mesh.getFaceCount();
        const size_t faceDim   = mesh.getDimension() - 1;

        res.clear();
        for (Index fi = 0; fi < faceCount; fi++)
        {
          if (essBdr.empty() && !mesh.isBoundary(fi))
            continue;

          if (!essBdr.empty())
          {
            const auto a = mesh.getAttribute(faceDim, fi);
            if (!a || !essBdr.count(*a))
              continue;
          }

          const auto& feU = fesU.getFiniteElement(faceDim, fi);
          const auto& feV = fesV.getFiniteElement(faceDim, fi);
          const auto& slaveDOFs = fesU.getDOFs(faceDim, fi);
          const auto& masterDOFs = fesV.getDOFs(faceDim, fi);

          const Index nMasters = static_cast<Index>(feV.getCount());

          for (Index s = 0; s < static_cast<Index>(feU.getCount()); s++)
          {
            const Index slave = slaveDOFs[s];
            if (res.find(slave) != res.end())
              continue;

            // Compute the constraint row C_{sj} = ℓ_s^u(A(φ_j^v)) for each
            // master j. We construct a callable that, given a Geometry::Point
            // p, evaluates A(φ_j^v)(p) through a pointwise IntegrationPoint
            // with no quadrature formula and pulls the j-th basis. The
            // slave-FES pullback drives this
            // through whatever evaluation pattern the slave DOF functional
            // uses (point evaluation for Lagrange, integral for moment-based
            // elements, etc.) — making this work for any FES whose
            // FiniteElement::LinearForm operates on a (Geometry::Point ->
            // value) callable through the FES pullback.
            std::vector<Index>  mIdx;
            std::vector<Scalar> mCoef;
            mIdx.reserve(static_cast<size_t>(nMasters));
            mCoef.reserve(static_cast<size_t>(nMasters));

            for (Index j = 0; j < nMasters; j++)
            {
              auto basisCallable = [&Av, j]
                                   (const Geometry::Point& p)
              {
                const Variational::IntegrationPoint ip(p);
                Av.setIntegrationPoint(ip);
                return Av.getBasis(static_cast<size_t>(j));
              };
              const auto mapping =
                fesU.getPullback({faceDim, fi}, std::move(basisCallable));
              const Scalar c = static_cast<Scalar>(feU.getLinearForm(s)(mapping));
              if (c != Scalar(0))
              {
                mIdx.push_back(masterDOFs[j]);
                mCoef.push_back(c);
              }
            }

            if (mIdx.empty())
              continue;

            const Index n = static_cast<Index>(mIdx.size());
            IndexArray masters(n);
            Math::Vector<Scalar> coeffs(n);
            for (Index k = 0; k < n; k++)
            {
              masters.coeffRef(k) = mIdx[static_cast<size_t>(k)];
              coeffs.coeffRef(k)  = mCoef[static_cast<size_t>(k)];
            }
            res.emplace(slave,
                std::pair{ std::move(masters), std::move(coeffs) });
          }
        }
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif
