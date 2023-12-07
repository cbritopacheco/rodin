/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_SERIAL_H
#define RODIN_ASSEMBLY_SERIAL_H

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Variational/BilinearForm.h"

#include "ForwardDecls.h"
#include "AssemblyBase.h"

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES>
  class Serial<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>
  {
    /**
     * @internal
     */
    static void add(
        std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols)
    {
      assert(rows.size() >= 0);
      assert(cols.size() >= 0);
      assert(in.rows() == rows.size());
      assert(in.cols() == cols.size());
      for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
      {
        for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
        {
          const Scalar s = in(i, j);
          if (s != Scalar(0))
            out.emplace_back(rows(i), cols(j), s);
        }
      }
    }

    public:
      using Parent =
        AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>>;
      using Input = typename Parent::Input;
      using OperatorType = std::vector<Eigen::Triplet<Scalar>>;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
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
        for (auto& bfi : input.bfis)
        {
          const auto& attrs = bfi.getAttributes();
          switch (bfi.getRegion())
          {
            case Variational::Integrator::Region::Domain:
            {
              for (auto it = input.mesh.getCell(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Faces:
            {
              for (auto it = input.mesh.getFace(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Boundary:
            {
              for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Interface:
            {
              for (auto it = input.mesh.getInterface(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
          }
        }
        return res;
      }

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief Serial assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Serial<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>
  {
    public:
      using Parent = AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix>>;
      using Input = typename Parent::Input;
      using OperatorType = Math::SparseMatrix;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const Input& input) const override
      {
        Serial<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Scalar>>>> assembly;
        const auto triplets = assembly.execute({ input.mesh, input.trialFES, input.testFES, input.bfis });
        OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
        return res;
      }

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief Serial assembly of the Math::SparseMatrix associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Serial<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
    : public AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>
  {
    static void add(
        Math::Matrix& out, const Math::Matrix& in,
        const IndexArray& rows, const IndexArray& cols)
    {
      assert(rows.size() >= 0);
      assert(cols.size() >= 0);
      assert(in.rows() == rows.size());
      assert(in.cols() == cols.size());
      for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
      {
        for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
          out(rows(i), cols(j)) += in(i, j);
      }
    }

    public:
      using Parent = AssemblyBase<Variational::BilinearForm<TrialFES, TestFES, Math::Matrix>>;
      using Input = typename Parent::Input;
      using OperatorType = Math::Matrix;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the linear operator
       * associated to the bilinear form.
       */
      OperatorType execute(const Input& input) const override
      {
        Math::Matrix res(input.testFES.getSize(), input.trialFES.getSize());;
        res.setZero();
        for (auto& bfi : input.bfis)
        {
          const auto& attrs = bfi.getAttributes();
          switch (bfi.getRegion())
          {
            case Variational::Integrator::Region::Domain:
            {
              for (auto it = input.mesh.getCell(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Faces:
            {
              for (auto it = input.mesh.getFace(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Boundary:
            {
              for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
            case Variational::Integrator::Region::Interface:
            {
              for (auto it = input.mesh.getInterface(); !it.end(); ++it)
              {
                if (attrs.size() == 0 || attrs.count(it->getAttribute()))
                {
                  const size_t d = it->getDimension();
                  const size_t i = it->getIndex();
                  const auto& trialDOFs = input.trialFES.getDOFs(d, i);
                  const auto& testDOFs = input.testFES.getDOFs(d, i);
                  bfi.assemble(*it);
                  add(res, bfi.getMatrix(), testDOFs, trialDOFs);
                }
              }
              break;
            }
          }
        }
        return res;
      }

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };

  /**
   * @brief %Serial assembly of the Math::Vector associated to a LinearFormBase
   * object.
   */
  template <>
  class Serial<Variational::LinearFormBase<Math::Vector>>
    : public AssemblyBase<Variational::LinearFormBase<Math::Vector>>
  {

    static void add(Math::Vector& out, const Math::Vector& in, const IndexArray& s);

    public:
      using Parent = AssemblyBase<Variational::LinearFormBase<Math::Vector>>;
      using VectorType = Math::Vector;

      Serial() = default;

      Serial(const Serial& other)
        : Parent(other)
      {}

      Serial(Serial&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Executes the assembly and returns the vector associated to the
       * linear form.
       */
      VectorType execute(const Input& input) const override;

      Serial* copy() const noexcept override
      {
        return new Serial(*this);
      }
  };
}

#endif

