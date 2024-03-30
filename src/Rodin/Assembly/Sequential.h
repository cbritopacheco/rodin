/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_SEQUENTIAL_H
#define RODIN_ASSEMBLY_SEQUENTIAL_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Kernels.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Variational/BilinearForm.h"

#include "Rodin/Utility/Repeat.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"

namespace Rodin::Assembly::Internal
{
  class SequentialIteration
  {
    public:
      SequentialIteration(const Geometry::MeshBase& mesh, Variational::Integrator::Region);

      Geometry::PolytopeIterator getIterator() const;

    private:
      std::reference_wrapper<const Geometry::MeshBase> m_mesh;
      Variational::Integrator::Region m_region;
  };
}

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES>
  class Sequential<
    std::vector<Eigen::Triplet<Scalar>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> final
    : public AssemblyBase<
        std::vector<Eigen::Triplet<Scalar>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
  {
    public:
      using Parent =
        AssemblyBase<
          std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>;
      using InputType = typename Parent::InputType;
      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

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
      OperatorType execute(const InputType& input) const override
      {
        std::vector<Eigen::Triplet<Scalar>> res;
        const auto& mesh = input.getTrialFES().getMesh();
        res.reserve(input.getTestFES().getSize() * std::log(input.getTrialFES().getSize()));
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.getTrialFES().getDOFs(d, i);
              const auto& testDOFs = input.getTestFES().getDOFs(d, i);
              bfi.assemble(*it);
              Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              Internal::SequentialIteration trialseq(mesh, bfi.getTrialRegion());
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  const auto& trialDOFs = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  const auto& testDOFs = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  bfi.assemble(*trIt, *teIt);
                  Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
            }
          }
        }
        return res;
      }

      inline
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::SparseMatrix,
    Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>> final
    : public AssemblyBase<
        Math::SparseMatrix,
        Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::SparseMatrix,
          Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>;
      using InputType = typename Parent::InputType;
      using OperatorType = Math::SparseMatrix;

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
      OperatorType execute(const InputType& input) const override
      {
        Sequential<
          std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> assembly;
        const auto triplets =
          assembly.execute({
            input.getTrialFES(), input.getTestFES(),
            input.getLocalBFIs(), input.getGlobalBFIs() });
        OperatorType res(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
      }

      inline
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  /**
   * @brief Sequential assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::Matrix,
    Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>> final
    : public AssemblyBase<
        Math::Matrix,
        Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::Matrix,
          Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>;
      using InputType = typename Parent::InputType;
      using OperatorType = Math::Matrix;

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
      OperatorType execute(const InputType& input) const override
      {
        Math::Matrix res(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setZero();
        const auto& mesh = input.getTrialFES().getMesh();
        for (auto& bfi : input.getLocalBFIs())
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.getTrialFES().getDOFs(d, i);
              const auto& testDOFs = input.getTestFES().getDOFs(d, i);
              bfi.assemble(*it);
              Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
        }
        for (auto& bfi : input.getGlobalBFIs())
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::SequentialIteration trialseq(mesh, bfi.getTrialRegion());
          Internal::SequentialIteration testseq(mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  const auto& trialDOFs = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  const auto& testDOFs = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  bfi.assemble(*trIt, *teIt);
                  Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
            }
          }
        }
        return res;
      }

      inline
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class ... TrialFES, class ... TestFES>
  class Sequential<
    std::vector<Eigen::Triplet<Scalar>>,
    Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>...>> final
      : public AssemblyBase<
          std::vector<Eigen::Triplet<Scalar>>,
          Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>...>>
  {
    public:
      using Parent =
        AssemblyBase<
          std::vector<Eigen::Triplet<Scalar>>,
          Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>...>>;
      using InputType = typename Parent::InputType;
      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;
      using Offsets = typename InputType::Offsets;

      Sequential() = default;

      Sequential(const Sequential& other)
        : Parent(other)
      {}

      Sequential(Sequential&& other)
        : Parent(std::move(other))
      {}

      OperatorType execute(const InputType& input) const override
      {
        using AssemblyTuple =
          Tuple<Sequential<std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>...>;

        AssemblyTuple assembly;

        const auto& t = input.getTuple();

        // Compute each block of triplets
        std::array<std::vector<Eigen::Triplet<Scalar>>, AssemblyTuple::Size> ts;
        assembly.zip(t)
                .map([](const auto& p)
                     { return p.first().execute(p.second()); })
                .iapply([&](const Index i, auto& v)
                        { ts[i] = std::move(v); });

        // Add the triplets with the offsets
        std::vector<Eigen::Triplet<Scalar>> res;
        size_t capacity = 0;
        for (const auto& v : ts)
          capacity += v.size();
        res.reserve(capacity);

        const Offsets& offsets = input.getOffsets();
        for (size_t i = 0; i < ts.size(); i++)
        {
          for (const Eigen::Triplet<Scalar>& t : ts[i])
          {
            res.emplace_back(
                t.row() + offsets[i].second(), t.col() + offsets[i].first(), t.value());
          }
        }
        return res;
      }

      inline
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };

  template <class ... TrialFES, class ... TestFES>
  class Sequential<
    Math::SparseMatrix,
    Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>...>> final
      : public AssemblyBase<
          Math::SparseMatrix,
          Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>...>>
    {
      public:
        using Parent =
          AssemblyBase<
            Math::SparseMatrix,
            Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>...>>;

        using InputType = typename Parent::InputType;

        using OperatorType = Math::SparseMatrix;

        Sequential() = default;

        Sequential(const Sequential& other)
          : Parent(other)
        {}

        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        OperatorType execute(const InputType& input) const override
        {
          Sequential<
            std::vector<Eigen::Triplet<Scalar>>,
            Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>...>> assembly;
          OperatorType res(input.getRows(), input.getColumns());
          const auto triplets = assembly.execute(input);
          res.setFromTriplets(triplets.begin(), triplets.end());
          return res;
        }

        inline
        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  template <class ... FES>
  class Sequential<
    Math::Vector,
    Tuple<Variational::LinearForm<FES, Math::Vector>...>> final
      : public AssemblyBase<
          Math::Vector,
          Tuple<Variational::LinearForm<FES, Math::Vector>...>>
    {
      public:
        using Parent =
          AssemblyBase<
            Math::Vector,
            Tuple<Variational::LinearForm<FES, Math::Vector>...>>;

        using InputType = typename Parent::InputType;

        using VectorType = Math::Vector;

        Sequential() = default;

        Sequential(const Sequential& other)
          : Parent(other)
        {}

        Sequential(Sequential&& other)
          : Parent(std::move(other))
        {}

        VectorType execute(const InputType& input) const override
        {
          using AssemblyTuple =
            Tuple<Sequential<Math::Vector, Variational::LinearForm<FES, Math::Vector>>...>;

          AssemblyTuple assembly;

          const auto& t = input.getTuple();

          // Compute each block of triplets
          auto vs = assembly.zip(t)
                            .map([](const auto& p) { return p.first().execute(p.second()); });

          Math::Vector res = Math::Vector::Zero(input.getSize());
          const auto& offsets = input.getOffsets();
          vs.iapply(
              [&](size_t i, const auto& v)
              {
                res.segment(offsets[i], v.size()) = v;
              });
          return res;
        }

        inline
        Sequential* copy() const noexcept override
        {
          return new Sequential(*this);
        }
    };

  /**
   * @brief %Sequential assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <class FES>
  class Sequential<Math::Vector, Variational::LinearForm<FES, Math::Vector>> final
    : public AssemblyBase<Math::Vector, Variational::LinearForm<FES, Math::Vector>>
  {
    public:
      using Parent = AssemblyBase<Math::Vector, Variational::LinearForm<FES, Math::Vector>>;
      using InputType = typename Parent::InputType;
      using VectorType = Math::Vector;

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
      VectorType execute(const InputType& input) const override
      {
        VectorType res(input.getFES().getSize());
        res.setZero();
        const auto& mesh = input.getFES().getMesh();
        for (auto& lfi : input.getLFIs())
        {
          const auto& attrs = lfi.getAttributes();
          Internal::SequentialIteration seq(mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.getFES().getDOFs(d, i);
              lfi.assemble(*it);
              Math::Kernels::add(res, lfi.getVector(), dofs);
            }
          }
        }
        return res;
      }

      inline
      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif

