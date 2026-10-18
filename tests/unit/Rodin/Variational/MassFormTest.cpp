/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Counts scalar H1 basis evaluations through the named-kernel interface.
   *
   * ## Architecture
   *
   * The space retains the production H1 mesh and basis definitions. Borrowed
   * element and basis views only count evaluations; they do not change values,
   * orders, geometry, or quadrature. Each space owns its independent count.
   *
   * @tparam K Field degree.
   * @tparam Scalar Real or complex field scalar.
   */
  template <size_t K, class Scalar>
  class CountingMassSpace final : public Variational::H1<K, Scalar>
  {
    public:
      /// @brief Production scalar space.
      using Parent = Variational::H1<K, Scalar>;
      /// @brief Production scalar finite element.
      using BaseElement = Variational::H1Element<K, Scalar>;

      /// @brief Borrowed finite-element view with counted basis evaluations.
      class Element
      {
        public:
          /// @brief A production basis with a borrowed evaluation counter.
          class Basis
          {
            public:
              /**
               * @brief Borrows a production basis by its local index.
               * @param element Production finite element.
               * @param index Local basis index.
               * @param evaluations Owning space's evaluation count.
               */
              Basis(const BaseElement& element, Index index, size_t& evaluations)
                : m_element(element),
                  m_index(index),
                  m_evaluations(evaluations)
              {}

              /**
               * @brief Counts and evaluates the unchanged production basis.
               * @param reference Reference coordinates.
               * @returns Production scalar basis value.
               */
              Scalar operator()(const Math::SpatialVector<Real>& reference) const
              {
                ++m_evaluations.get();
                return m_element.get().getBasis(m_index)(reference);
              }

            private:
              std::reference_wrapper<const BaseElement> m_element;
              Index m_index;
              std::reference_wrapper<size_t> m_evaluations;
          };

          /**
           * @brief Borrows an element and its space's count.
           * @param element Production finite element.
           * @param evaluations Owning space's evaluation count.
           */
          Element(const BaseElement& element, size_t& evaluations)
            : m_element(element),
              m_evaluations(evaluations)
          {}

          /// @brief Returns the unchanged local basis count.
          size_t getCount() const
          {
            return m_element.get().getCount();
          }
          /// @brief Returns the unchanged total polynomial order.
          size_t getOrder() const
          {
            return m_element.get().getOrder();
          }

          /**
           * @brief Borrows the indexed basis without evaluating it.
           * @param index Local basis index.
           * @returns Counted basis view.
           */
          Basis getBasis(Index index) const
          {
            return Basis(m_element.get(), index, m_evaluations.get());
          }

        private:
          std::reference_wrapper<const BaseElement> m_element;
          std::reference_wrapper<size_t> m_evaluations;
      };

      /**
       * @brief Constructs the production space with an independent zero count.
       * @param mesh Mesh retained by the production space.
       */
      explicit CountingMassSpace(const Geometry::LocalMesh& mesh)
        : Parent(std::integral_constant<size_t, K>{}, mesh),
          m_evaluations(0)
      {}

      /**
       * @brief Borrows a production element through its counted view.
       * @param dimension Entity dimension.
       * @param index Local entity index.
       * @returns Counted finite-element view.
       */
      Element getFiniteElement(size_t dimension, Index index) const
      {
        return Element(Parent::getFiniteElement(dimension, index), m_evaluations);
      }

      /// @brief Returns the number of actual scalar basis evaluations.
      size_t getEvaluationCount() const
      {
        return m_evaluations;
      }

    private:
      mutable size_t m_evaluations;
  };
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Retains production H1 metadata for the counted test space.
   * @tparam K Field degree.
   * @tparam Scalar Real or complex field scalar.
   */
  template <size_t K, class Scalar>
  struct Traits<Tests::Unit::CountingMassSpace<K, Scalar>>
    : Traits<Variational::H1<K, Scalar>>
  {};
}

namespace Rodin::Variational
{
  /**
   * @brief Retains production leaf types for the counted kernel test space.
   *
   * The named form exposes its trial/test leaf types even when only its kernel
   * is instantiated. Their shape-function interface is the unchanged H1 one;
   * only direct finite-element calls through CountingMassSpace are observed.
   *
   * @tparam Derived Trial/test leaf type.
   * @tparam K Field degree.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Space Trial or test grading.
   */
  template <class Derived, size_t K, class Scalar, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, Tests::Unit::CountingMassSpace<K, Scalar>, Space>
    : public ShapeFunction<Derived, H1<K, Scalar>, Space>
  {
    public:
      /// @brief Production H1 shape-function interface.
      using Parent = ShapeFunction<Derived, H1<K, Scalar>, Space>;
      /// @brief Counted space retained by the named leaf interface.
      using FESType = Tests::Unit::CountingMassSpace<K, Scalar>;

      /**
       * @brief Borrows the counted space and its production shape interface.
       * @param space Counted finite-element space.
       */
      explicit ShapeFunction(const FESType& space)
        : Parent(space),
          m_space(space)
      {}

      /// @brief Returns the counted space without erasing its observer type.
      const FESType& getFiniteElementSpace() const
      {
        return m_space.get();
      }

    private:
      std::reference_wrapper<const FESType> m_space;
  };
}

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks linear basis-evaluation counts and kernel reference lifetimes.
   *
   * Cold preparation requires @f$N_qN_r@f$ trial and @f$N_qN_t@f$ test
   * evaluations. Repeated computation uses the retained reference matrices;
   * copying a kernel constructs an independent reference cache. Matrix equality
   * is exact because all comparisons retain the production summation order.
   *
   * @tparam TrialDegree Trial H1 degree.
   * @tparam TestDegree Test H1 degree.
   * @tparam Scalar Real or complex scalar field.
   */
  template <size_t TrialDegree, size_t TestDegree, class Scalar>
  void checkMassBasisReuse()
  {
    using namespace Geometry;
    using namespace Variational;
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = dimension == 0
        ? LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize()
        : LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      for (size_t d = 1; d <= dimension; ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      CountingMassSpace<TrialDegree, Scalar> trial(mesh);
      CountingMassSpace<TestDegree, Scalar> test(mesh);
      using Solution = GridFunction<decltype(trial), Math::Vector<Scalar>>;
      using Form = MassForm<Solution, decltype(trial), decltype(test),
        Math::SparseMatrix<Scalar>, void>;
      typename Form::KernelType kernel(trial, test);
      const auto cell = *mesh.getCell();
      const auto trialElement = trial.getFiniteElement(dimension, cell.getIndex());
      const auto testElement = test.getFiniteElement(dimension, cell.getIndex());
      const auto& formula = QF::PolytopeQuadratureFormula::get(
        trialElement.getOrder() + testElement.getOrder(), geometry);
      const size_t trialCount = formula.getSize() * trialElement.getCount();
      const size_t testCount = formula.getSize() * testElement.getCount();
      Math::Matrix<Scalar> first, repeated, transferred;
      kernel.compute(first, cell);
      const auto& originalTrial =
        static_cast<const typename decltype(trial)::Parent&>(trial).getFiniteElement(
          dimension, cell.getIndex());
      const auto& originalTest =
        static_cast<const typename decltype(test)::Parent&>(test).getFiniteElement(
          dimension, cell.getIndex());
      Math::Matrix<Scalar> expected =
        Math::Matrix<Scalar>::Zero(originalTest.getCount(), originalTrial.getCount());
      Math::Matrix<Scalar> reference(expected.rows(), expected.cols());
      const auto& quadrature = cell.getQuadrature(formula);
      for (size_t qp = 0; qp < formula.getSize(); ++qp)
      {
        const auto& coordinates = formula.getPoint(qp);
        for (size_t te = 0; te < originalTest.getCount(); ++te)
        {
          const auto value = originalTest.getBasis(te)(coordinates);
          for (size_t tr = 0; tr < originalTrial.getCount(); ++tr)
            reference(te, tr) = Math::dot(originalTrial.getBasis(tr)(coordinates), value);
        }
        const Scalar weight = static_cast<Scalar>(
          formula.getWeight(qp) * quadrature.getPoint(qp).getDistortion());
        expected += weight * reference;
      }
      EXPECT_EQ(first, expected);
      EXPECT_EQ(trial.getEvaluationCount(), trialCount);
      EXPECT_EQ(test.getEvaluationCount(), testCount);
      kernel.compute(repeated, cell);
      EXPECT_EQ(trial.getEvaluationCount(), trialCount);
      EXPECT_EQ(test.getEvaluationCount(), testCount);
      EXPECT_EQ(first, repeated);
      typename Form::KernelType copy(kernel);
      copy.compute(transferred, cell);
      EXPECT_EQ(trial.getEvaluationCount(), 2 * trialCount);
      EXPECT_EQ(test.getEvaluationCount(), 2 * testCount);
      EXPECT_EQ(first, transferred);
    }
  }

  /// @brief Real equal-degree preparation evaluates each basis once per sample.
  TEST(Rodin_Variational_MassForm, RealEqualDegreeBasisReuse)
  {
    checkMassBasisReuse<1, 1, Real>();
  }

  /// @brief Complex equal-degree preparation preserves the same count contract.
  TEST(Rodin_Variational_MassForm, ComplexEqualDegreeBasisReuse)
  {
    checkMassBasisReuse<2, 2, Complex>();
  }

  /// @brief Rectangular trial/test reference matrices retain independent counts.
  TEST(Rodin_Variational_MassForm, MixedDegreeBasisReuse)
  {
    checkMassBasisReuse<1, 2, Real>();
    checkMassBasisReuse<2, 1, Complex>();
  }
}
