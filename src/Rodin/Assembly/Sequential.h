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
    std::vector<Eigen::Triplet<Real>>,
    Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>> final
    : public AssemblyBase<
        std::vector<Eigen::Triplet<Real>>,
        Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

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
      OperatorType execute(const InputType& input) const override
      {
        OperatorType res;
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
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
              for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
              {
                for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                {
                  const ScalarType s = bfi.integrate(m, l);
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
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                  {
                    for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                    {
                      const ScalarType s = bfi.integrate(m, l);
                      if (s != ScalarType(0))
                        res.emplace_back(rows(l), cols(m), s);
                    }
                  }
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
   * @brief Sequential assembly of the Math::SparseMatrix<Real> associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::SparseMatrix<Real>,
    Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>> final
    : public AssemblyBase<
        Math::SparseMatrix<Real>,
        Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

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
      OperatorType execute(const InputType& input) const override
      {
        Sequential<
          std::vector<Eigen::Triplet<Real>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>> assembly;
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
   * @brief Sequential assembly of the Math::SparseMatrix<Real> associated to a
   * BilinearFormBase object.
   */
  template <class TrialFES, class TestFES>
  class Sequential<
    Math::Matrix<Real>,
    Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>> final
    : public AssemblyBase<
        Math::Matrix<Real>,
        Variational::BilinearForm<TrialFES, TestFES, Math::Matrix<Real>>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = Math::Matrix<ScalarType>;

      using LocalBilinearFormIntegratorBaseType = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType = Variational::GlobalBilinearFormIntegratorBase<ScalarType>;

      using BilinearFormType = Variational::BilinearForm<TrialFES, TestFES, OperatorType>;

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
      OperatorType execute(const InputType& input) const override
      {
        OperatorType res(input.getTestFES().getSize(), input.getTrialFES().getSize());
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
              bfi.setPolytope(*it);
              const auto& rows = input.getTestFES().getDOFs(it.getDimension(), it->getIndex());
              const auto& cols = input.getTrialFES().getDOFs(it.getDimension(), it->getIndex());
              for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                  res(rows(l), cols(m)) += bfi.integrate(m, l);
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
                  bfi.setPolytope(*trIt, *teIt);
                  const auto& rows = input.getTestFES().getDOFs(teIt.getDimension(), teIt->getIndex());
                  const auto& cols = input.getTrialFES().getDOFs(trIt.getDimension(), trIt->getIndex());
                  for (size_t l = 0; l < static_cast<size_t>(rows.size()); l++)
                    for (size_t m = 0; m < static_cast<size_t>(cols.size()); m++)
                      res(rows(l), cols(m)) += bfi.integrate(m, l);
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
    std::vector<Eigen::Triplet<Real>>,
    Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>> final
      : public AssemblyBase<
          std::vector<Eigen::Triplet<Real>>,
          Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>>
  {
    public:
      using ScalarType = Real;

      using OperatorType = std::vector<Eigen::Triplet<ScalarType>>;

      using TupleType =
        Tuple<Variational::BilinearForm<TrialFES, TestFES, OperatorType>...>;

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

      OperatorType execute(const InputType& input) const override
      {
        using AssemblyTuple =
          Tuple<Sequential<std::vector<Eigen::Triplet<Real>>,
          Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>>...>;

        AssemblyTuple assembly;

        const auto& t = input.getTuple();

        // Compute each block of triplets
        std::array<std::vector<Eigen::Triplet<Real>>, AssemblyTuple::Size> ts;
        assembly.zip(t)
                .map([](const auto& p)
                     { return p.first().execute(p.second()); })
                .iapply([&](const Index i, auto& v)
                        { ts[i] = std::move(v); });

        // Add the triplets with the offsets
        std::vector<Eigen::Triplet<Real>> res;
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
    Math::SparseMatrix<Real>,
    Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>...>> final
      : public AssemblyBase<
          Math::SparseMatrix<Real>,
          Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>...>>
    {
      public:
        using Parent =
          AssemblyBase<
            Math::SparseMatrix<Real>,
            Tuple<Variational::BilinearForm<TrialFES, TestFES, Math::SparseMatrix<Real>>...>>;

        using InputType = typename Parent::InputType;

        using OperatorType = Math::SparseMatrix<Real>;

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
            std::vector<Eigen::Triplet<Real>>,
            Tuple<Variational::BilinearForm<TrialFES, TestFES, std::vector<Eigen::Triplet<Real>>>...>> assembly;
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

        VectorType execute(const InputType& input) const override
        {
          using AssemblyTuple =
            Tuple<Sequential<Math::Vector<Real>, Variational::LinearForm<FES, Math::Vector<Real>>>...>;

          AssemblyTuple assembly;

          const auto& t = input.getTuple();

          // Compute each block of triplets
          auto vs = assembly.zip(t)
                            .map([](const auto& p) { return p.first().execute(p.second()); });

          Math::Vector<Real> res = Math::Vector<Real>::Zero(input.getSize());
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
   * @brief %Sequential assembly of the Math::Vector<Real> associated to a LinearFormBase
   * object.
   */
  template <class FES>
  class Sequential<Math::Vector<Real>, Variational::LinearForm<FES, Math::Vector<Real>>> final
    : public AssemblyBase<Math::Vector<Real>, Variational::LinearForm<FES, Math::Vector<Real>>>
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
              lfi.setPolytope(*it);
              const size_t d = it.getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.getFES().getDOFs(d, i);
              for (size_t l = 0; l < static_cast<size_t>(dofs.size()); l++)
                res(dofs(l)) += lfi.integrate(l);
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

