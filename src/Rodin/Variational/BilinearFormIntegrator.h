#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  class BilinearFormIntegratorBase : public Integrator
  {
    public:
      using Parent = Integrator;

      BilinearFormIntegratorBase(
          const ShapeFunctionBase<TrialSpace>& u,
          const ShapeFunctionBase<TestSpace>& v)
        :  m_u(u.copy()), m_v(v.copy())
      {}

      BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
        :  Parent(other),
          m_u(other.m_u->copy()), m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
        :  Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      const std::set<Geometry::Attribute>& getAttributes() const
      {
        return m_attrs;
      }

      /**
       * @brief Gets reference to trial function.
       */
      const ShapeFunctionBase<TrialSpace>& getTrialFunction() const
      {
        assert(m_u);
        return *m_u;
      }

      /**
       * @brief Gets reference to test function.
       */
      const ShapeFunctionBase<TestSpace>& getTestFunction() const
      {
        assert(m_v);
        return *m_v;
      }

      /**
       * @brief Specifies the material reference over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material reference over which the integration should
       * take place.
       */
      BilinearFormIntegratorBase& over(Geometry::Attribute attr)
      {
        return over(std::set<Geometry::Attribute>{attr});
      }

      /**
       * @brief Specifies the material references over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material references over which the integration should
       * take place.
       */
      BilinearFormIntegratorBase& over(const std::set<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      Integrator::Type getType() const override
      {
        return Integrator::Type::Bilinear;
      }

      virtual ~BilinearFormIntegratorBase() = default;

      /**
       * @brief Performs the assembly of the element matrix for the given
       * element.
       */
      virtual Math::Matrix getMatrix(const Geometry::Simplex& element) const = 0;

      virtual BilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_u;
      std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_v;
      std::set<Geometry::Attribute> m_attrs;
  };
}

#endif
