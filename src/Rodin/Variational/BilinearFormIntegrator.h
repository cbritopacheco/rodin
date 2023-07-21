#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  /**
   * @brief Abstract base class for bilinear form integrators.
   *
   * This class provides the base functionality for bilinear form integrator
   * objects.
   */
  class BilinearFormIntegratorBase : public Integrator
  {
    public:
      /// Parent class
      using Parent = Integrator;

      /**
       * @brief Constructs the object given a TrialFunction and a TestFunction.
       */
      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v)
      {}

      /**
       * @brief Deleted constructor.
       */
      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(TrialFunction<TrialFES>&& u, const TestFunction<TestFES>& v) = delete;

      /**
       * @brief Deleted constructor.
       */
      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(const TrialFunction<TrialFES>& u, TestFunction<TestFES>&& v) = delete;

      /**
       * @brief Deleted constructor.
       */
      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(TrialFunction<TrialFES>&& u, TestFunction<TestFES>&& v) = delete;

      /**
       * @brief Copy constructor.
       */
      BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_attrs(other.m_attrs),
          m_matrix(other.m_matrix)
      {}

      /**
       * @brief Move constructor.
       */
      BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs)),
          m_matrix(std::move(other.m_matrix))
      {}

      virtual
      ~BilinearFormIntegratorBase() = default;

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      inline
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_attrs;
      }


      /**
       * @brief Specifies the material reference over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material reference over which the integration should
       * take place.
       */
      inline
      BilinearFormIntegratorBase& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      /**
       * @brief Specifies the material references over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material references over which the integration should
       * take place.
       */
      inline
      BilinearFormIntegratorBase& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      inline
      Integrator::Type getType() const
      final override
      {
        return Integrator::Type::Bilinear;
      }

      /**
       * @brief Gets a constant reference to trial function object.
       */
      inline
      const FormLanguage::Base& getTrialFunction() const
      {
        return m_u.get();
      }

      /**
       * @brief Gets a constant reference to test function object.
       */
      inline
      const FormLanguage::Base& getTestFunction() const
      {
        return m_v.get();
      }

      inline
      const Math::Matrix& getMatrix() const
      {
        return m_matrix;
      }

      /**
       * @returns The element matrix of size of @f$ m \times n @f$ where @f$ n
       * @f$ (resp. @f$ m @f$) denotes the number of degrees of freedom on the
       * polytope for the test (trial) space.
       */
      inline
      Math::Matrix& getMatrix()
      {
        return m_matrix;
      }

      /**
       * @brief Performs the assembly of the element matrix for the given
       * element.
       *
       * Assembles the stiffness matrix of the given element.
       *
       */
      virtual
      void assemble(const Geometry::Polytope& polytope) = 0;

      virtual
      BilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const FormLanguage::Base> m_u;
      std::reference_wrapper<const FormLanguage::Base> m_v;
      FlatSet<Geometry::Attribute> m_attrs;
      Math::Matrix m_matrix;
  };
}

#endif
