/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_NAMEDFORMKERNEL_H
#define RODIN_ASSEMBLY_NAMEDFORMKERNEL_H

#include <algorithm>
#include <functional>
#include <vector>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Assembly
{
  template <class Scalar, class TrialFES, class TestFES>
  class MassFormKernel
  {
    public:
      using ScalarType = Scalar;
      using MatrixType = Math::Matrix<ScalarType>;

      MassFormKernel(const TrialFES& trialFES, const TestFES& testFES)
        : m_trialFES(trialFES), m_testFES(testFES)
      {}

      void compute(MatrixType& out, const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
        const auto& testFE = m_testFES.get().getFiniteElement(d, i);
        const size_t order = trialFE.getOrder() + testFE.getOrder();
        const auto& qf = QF::PolytopeQuadratureFormula::get(
          order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);

        const bool rebuild =
          !m_cached || m_geometry != polytope.getGeometry()
          || m_order != order
          || m_trialCount != trialFE.getCount()
          || m_testCount != testFE.getCount();
        if (rebuild)
        {
          m_cached = true;
          m_geometry = polytope.getGeometry();
          m_order = order;
          m_trialCount = trialFE.getCount();
          m_testCount = testFE.getCount();
          m_reference.resize(qf.getSize());
          for (size_t qp = 0; qp < qf.getSize(); ++qp)
          {
            auto& matrix = m_reference[qp];
            matrix.resize(testFE.getCount(), trialFE.getCount());
            const auto& reference = qf.getPoint(qp);
            for (size_t te = 0; te < testFE.getCount(); ++te)
            {
              const auto testValue = testFE.getBasis(te)(reference);
              for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
                matrix(te, tr) = Math::dot(
                  trialFE.getBasis(tr)(reference), testValue);
            }
          }
        }

        out.resize(testFE.getCount(), trialFE.getCount());
        out.setZero();
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const ScalarType weight = static_cast<ScalarType>(
            qf.getWeight(qp) * point.getDistortion());
          out += weight * m_reference[qp];
        }
      }

    private:
      std::reference_wrapper<const TrialFES> m_trialFES;
      std::reference_wrapper<const TestFES> m_testFES;
      mutable bool m_cached = false;
      mutable Geometry::Polytope::Type m_geometry =
        Geometry::Polytope::Type::Point;
      mutable size_t m_order = 0;
      mutable size_t m_trialCount = 0;
      mutable size_t m_testCount = 0;
      mutable std::vector<MatrixType> m_reference;
  };

  template <class Scalar, class TrialFES, class TestFES>
  class DiffusionFormKernel
  {
    public:
      using ScalarType = Scalar;
      using MatrixType = Math::Matrix<ScalarType>;

      DiffusionFormKernel(const TrialFES& trialFES, const TestFES& testFES)
        : m_trialFES(trialFES), m_testFES(testFES)
      {}

      void compute(MatrixType& out, const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const Index i = polytope.getIndex();
        const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
        const auto& testFE = m_testFES.get().getFiniteElement(d, i);
        const size_t order =
          trialFE.getOrder() == 0 || testFE.getOrder() == 0
            ? 0
            : trialFE.getOrder() + testFE.getOrder() - 2;
        const auto& qf = QF::PolytopeQuadratureFormula::get(
          order, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);

        const bool rebuild =
          !m_cached || m_geometry != polytope.getGeometry()
          || m_order != order
          || m_trialCount != trialFE.getCount()
          || m_testCount != testFE.getCount();
        if (rebuild)
        {
          m_cached = true;
          m_geometry = polytope.getGeometry();
          m_order = order;
          m_trialCount = trialFE.getCount();
          m_testCount = testFE.getCount();
          m_trialGradients.assign(qf.getSize(), {});
          m_testGradients.assign(qf.getSize(), {});
          for (size_t qp = 0; qp < qf.getSize(); ++qp)
          {
            const auto& reference = qf.getPoint(qp);
            auto& trialGradients = m_trialGradients[qp];
            auto& testGradients = m_testGradients[qp];
            trialGradients.reserve(trialFE.getCount());
            testGradients.reserve(testFE.getCount());
            for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
              trialGradients.emplace_back(
                trialFE.getBasis(tr).getGradient()(reference));
            for (size_t te = 0; te < testFE.getCount(); ++te)
              testGradients.emplace_back(
                testFE.getBasis(te).getGradient()(reference));
          }
        }

        out.resize(testFE.getCount(), trialFE.getCount());
        out.setZero();
        // The pullback of the physical gradients is folded into the metric
        // G = J^{-1} J^{-T}, so that only the trial gradients are mapped.
        m_scratch.resize(trialFE.getCount());
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const ScalarType weight = static_cast<ScalarType>(
            qf.getWeight(qp) * point.getDistortion());
          const auto& jacobianInverse = point.getJacobianInverse();
          const auto metric = jacobianInverse * jacobianInverse.transpose();
          const auto& trialGradients = m_trialGradients[qp];
          const auto& testGradients = m_testGradients[qp];
          for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
            m_scratch[tr] = metric * trialGradients[tr];
          for (size_t te = 0; te < testFE.getCount(); ++te)
          {
            for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
              out(te, tr) +=
                weight * Math::dot(testGradients[te], m_scratch[tr]);
          }
        }
      }

    private:
      std::reference_wrapper<const TrialFES> m_trialFES;
      std::reference_wrapper<const TestFES> m_testFES;
      mutable bool m_cached = false;
      mutable Geometry::Polytope::Type m_geometry =
        Geometry::Polytope::Type::Point;
      mutable size_t m_order = 0;
      mutable size_t m_trialCount = 0;
      mutable size_t m_testCount = 0;
      mutable std::vector<std::vector<Math::SpatialVector<ScalarType>>>
        m_trialGradients;
      mutable std::vector<std::vector<Math::SpatialVector<ScalarType>>>
        m_testGradients;
      mutable std::vector<Math::SpatialVector<ScalarType>> m_scratch;
  };
}

#endif
