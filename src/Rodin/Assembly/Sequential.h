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
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<
        std::vector<Eigen::Triplet<Scalar>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
  {
    public:
      using Parent =
        AssemblyBase<
          std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>;
      using Input = typename Parent::Input;
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
      OperatorType execute(const Input& input) const override
      {
        std::vector<Eigen::Triplet<Scalar>> res;
        res.reserve(input.testFES.getSize() * std::log(input.trialFES.getSize()));
        for (auto& bfi : input.lbfis)
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(input.mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
        }
        for (auto& bfi : input.gbfis)
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::SequentialIteration testseq(input.mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              Internal::SequentialIteration trialseq(input.mesh, bfi.getTrialRegion());
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  const auto& trialDOFs = input.trialFES.getDOFs(trIt.getDimension(), trIt->getIndex());
                  const auto& testDOFs = input.testFES.getDOFs(teIt.getDimension(), teIt->getIndex());
                  bfi.assemble(*trIt, *teIt);
                  Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
            }
          }
        }
        return res;
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
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::SparseMatrix,
    Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
    : public AssemblyBase<
        Math::SparseMatrix,
        Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::SparseMatrix,
          Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>;
      using Input = typename Parent::Input;
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
      OperatorType execute(const Input& input) const override
      {
        Sequential<
          std::vector<Eigen::Triplet<Scalar>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> assembly;
        const auto triplets =
          assembly.execute({
            input.mesh,
            input.trialFES, input.testFES,
            input.lbfis, input.gbfis });
        OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
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
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::Matrix,
    Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
    : public AssemblyBase<
        Math::Matrix,
        Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
  {
    public:
      using Parent =
        AssemblyBase<
          Math::Matrix,
          Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>;
      using Input = typename Parent::Input;
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
      OperatorType execute(const Input& input) const override
      {
        Math::Matrix res(input.testFES.getSize(), input.trialFES.getSize());
        res.setZero();
        for (auto& bfi : input.lbfis)
        {
          const auto& attrs = bfi.getAttributes();
          Internal::SequentialIteration seq(input.mesh, bfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
        }
        for (auto& bfi : input.gbfis)
        {
          const auto& trialAttrs = bfi.getTrialAttributes();
          const auto& testAttrs = bfi.getTestAttributes();
          Internal::SequentialIteration trialseq(input.mesh, bfi.getTrialRegion());
          Internal::SequentialIteration testseq(input.mesh, bfi.getTestRegion());
          for (auto teIt = testseq.getIterator(); teIt; ++teIt)
          {
            if (testAttrs.size() == 0 || testAttrs.count(teIt->getAttribute()))
            {
              for (auto trIt = trialseq.getIterator(); trIt; ++trIt)
              {
                if (trialAttrs.size() == 0 || trialAttrs.count(trIt->getAttribute()))
                {
                  const auto& trialDOFs = input.trialFES.getDOFs(trIt.getDimension(), trIt->getIndex());
                  const auto& testDOFs = input.testFES.getDOFs(teIt.getDimension(), teIt->getIndex());
                  bfi.assemble(*trIt, *teIt);
                  Math::Kernels::add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
            }
          }
        }
        return res;
      }

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
  class Sequential<Math::Vector, Variational::LinearForm<FES, Math::Vector>>
    : public AssemblyBase<Math::Vector, Variational::LinearForm<FES, Math::Vector>>
  {
    public:
      using Parent = AssemblyBase<Math::Vector, Variational::LinearForm<FES, Math::Vector>>;
      using Input = typename Parent::Input;
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
      VectorType execute(const Input& input) const override
      {
        VectorType res(input.fes.getSize());
        res.setZero();
        for (auto& lfi : input.lfis)
        {
          const auto& attrs = lfi.getAttributes();
          Internal::SequentialIteration seq(input.mesh, lfi.getRegion());
          for (auto it = seq.getIterator(); it; ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              lfi.assemble(*it);
              Math::Kernels::add(res, lfi.getVector(), dofs);
            }
          }
        }
        return res;
      }

      Sequential* copy() const noexcept override
      {
        return new Sequential(*this);
      }
  };
}

#endif

