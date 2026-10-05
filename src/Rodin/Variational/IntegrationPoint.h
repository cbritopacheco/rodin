/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRATIONPOINT_H
#define RODIN_VARIATIONAL_INTEGRATIONPOINT_H

#include <cassert>
#include <optional>
#include <type_traits>

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Variational
{
  /**
   * @brief Evaluation point used by variational operators.
   *
   * An IntegrationPoint always stores a geometric point. It may optionally be
   * associated to a quadrature formula and quadrature index:
   * - `getQuadratureFormula() != nullptr`: the point is part of a quadrature
   *   loop, and `getIndex()` identifies the quadrature sample in that formula.
   * - `getQuadratureFormula() == nullptr`: the point represents a generic
   *   pointwise evaluation (outside a quadrature tabulation context). In this
   *   case, consumers are expected to evaluate using
   *   `getPoint().getReferenceCoordinates()` rather than quadrature tabulation.
   *
   * @note The quadrature formula pointer is non-owning. The referenced
   * quadrature object must outlive all IntegrationPoint uses that dereference
   * this pointer.
   */
  class IntegrationPoint
  {
    public:
      /// @brief Constructs an integration point from a geometric point.
      IntegrationPoint(const Geometry::Point& p)
        : m_p(p),
          m_qf(nullptr),
          m_qp(0)
      {}

      /**
       * @brief Constructs a quadrature-associated integration point.
       * @param[in] p Geometric point (usually built from quadrature sample coordinates)
       * @param[in] qf Non-owning pointer to the quadrature formula that owns
       *                the sample indexing, or @c nullptr for direct
       *                pointwise evaluation
       * @param[in] qp Quadrature sample index in @p qf
       */
      IntegrationPoint(
        const Geometry::Point& p, const QF::QuadratureFormulaBase* qf, size_t qp)
        : m_p(p),
          m_qf(qf),
          m_qp(qp)
      {}

      /**
       * @returns Underlying geometric point context.
       */
      const Geometry::Point& getPoint() const
      {
        return m_p;
      }

      /**
       * @returns Non-owning pointer to the quadrature formula, or @c nullptr
       * when this integration point is not in a quadrature context.
       */
      const QF::QuadratureFormulaBase* getQuadratureFormula() const
      {
        return m_qf;
      }

      /**
       * @returns Quadrature index associated with this point.
       * @note Meaningful only when getQuadratureFormula() is not @c nullptr.
       */
      size_t getIndex() const
      {
        return m_qp;
      }

    private:
      std::reference_wrapper<const Geometry::Point> m_p;
      const QF::QuadratureFormulaBase* m_qf;
      size_t m_qp;
  };

  namespace Internal
  {
    /**
     * @brief Whether a coefficient value of type @p Value owns its data and can
     * be held across basis evaluations: the declared range type of the
     * function, or a plain Eigen matrix. Lazy Eigen expressions, which refer
     * to temporaries, are excluded.
     */
    template <class Value, class Range>
    inline constexpr bool IsOwningCoefficient =
      std::is_same_v<Value, Range> ||
      std::is_base_of_v<Eigen::PlainObjectBase<Value>, Value>;

    /**
     * @brief Value of a coefficient function at the current quadrature point.
     *
     * A product of a function and a shape function evaluates the function at
     * every basis index. This holds the function value at the quadrature
     * point last passed to refresh(), so that it is evaluated once per point.
     * The value is used only while the point it was taken at is still the
     * current one (same mesh entity, quadrature formula and node); otherwise
     * the caller evaluates the function directly. With @p Plain false
     * (see IsOwningCoefficient) nothing is held.
     */
    template <class Value, bool Plain>
    class CoefficientCache
    {
      public:
        static constexpr bool Enabled = Plain &&
          std::is_copy_constructible_v<Value> && std::is_copy_assignable_v<Value>;

        CoefficientCache() = default;

        /// @brief Copies start empty: the value belongs to one evaluation pass.
        CoefficientCache(const CoefficientCache&) : CoefficientCache() {}

        CoefficientCache(CoefficientCache&&) = default;

        CoefficientCache& operator=(const CoefficientCache&)
        {
          m_value.reset();
          return *this;
        }

        CoefficientCache& operator=(CoefficientCache&&) = default;

        /// @brief Evaluates @p f at @p ip, if @p ip is a quadrature node.
        template <class F>
        void refresh(const F& f, const IntegrationPoint& ip)
        {
          if constexpr (Enabled)
          {
            if (!ip.getQuadratureFormula())
            {
              m_value.reset();
              return;
            }
            m_key = Key(ip);
            if (m_value)
              *m_value = f.getValue(ip);
            else
              m_value.emplace(f.getValue(ip));
          }
        }

        /// @brief The held value, or nullptr if it was not taken at @p ip.
        const Value* find(const IntegrationPoint& ip) const
        {
          if constexpr (Enabled)
          {
            if (m_value && ip.getQuadratureFormula() && m_key == Key(ip))
              return &*m_value;
          }
          return nullptr;
        }

      private:
        struct Key
        {
            Key() = default;

            explicit Key(const IntegrationPoint& ip)
            {
              const auto& polytope = ip.getPoint().getPolytope();
              mesh = static_cast<const void*>(&polytope.getMesh());
              dimension = polytope.getDimension();
              index = polytope.getIndex();
              qf = ip.getQuadratureFormula();
              qp = ip.getIndex();
            }

            bool operator==(const Key&) const = default;

            const void* mesh = nullptr;
            size_t dimension = 0;
            Index index = 0;
            const QF::QuadratureFormulaBase* qf = nullptr;
            size_t qp = 0;
        };

        Key m_key;
        std::optional<Value> m_value;
    };
  }
}

#endif
