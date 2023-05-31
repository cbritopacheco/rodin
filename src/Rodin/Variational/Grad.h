/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  // /**
  //  * @ingroup GradSpecializations
  //  */
  // template <class ... Ts>
  // class Grad<GridFunction<H1<Scalar, Ts...>>> final
  //   : public VectorFunctionBase<Grad<GridFunction<H1<Scalar, Ts...>>>>
  // {
  //   public:
  //     using Operand = GridFunction<H1<Scalar, Ts...>>;
  //     using Parent = VectorFunctionBase<Grad<GridFunction<H1<Scalar, Ts...>>>>;

  //     /**
  //      * @brief Constructs the gradient of an @f$ H^1 @f$ function
  //      * @f$ u @f$.
  //      * @param[in] u Grid function to be differentiated
  //      */
  //     Grad(const Operand& u)
  //       : m_u(u)
  //     {
  //       assert(u.getFiniteElementSpace().getVectorDimension() == 1);
  //     }

  //     Grad(const Grad& other)
  //       : Parent(other),
  //         m_u(other.m_u)
  //     {}

  //     Grad(Grad&& other)
  //       : Parent(std::move(other)),
  //         m_u(std::move(other.m_u))
  //     {}

  //     inline
  //     constexpr
  //     size_t getDimension() const
  //     {
  //       return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
  //     }

  //     inline
  //     constexpr
  //     Grad& traceOf(Geometry::Attribute attr)
  //     {
  //       Parent::traceOf(attr);
  //       return *this;
  //     }

  //     inline
  //     constexpr
  //     Grad& traceOf(const std::set<Geometry::Attribute>& attrs)
  //     {
  //       Parent::traceOf(attrs);
  //       return *this;
  //     }

  //     Math::Vector getValue(const Geometry::Point& p) const
  //     {
  //       assert(false);
  //       // const auto& simplex = p.getSimplex();
  //       // const auto& simplexMesh = simplex.getMesh();
  //       // const auto& fes = getOperand().getFiniteElementSpace();
  //       // const auto& fesMesh = fes.getMesh();
  //       // Math::Vector grad(getDimension());
  //       // mfem::Vector tmp(grad.data(), grad.size());
  //       // if (simplex.getDimension() == fesMesh.getDimension())
  //       // {
  //       //   const auto& trans = p.getTransformation();
  //       //   m_u.get().getHandle().GetGradient(trans.getHandle(), tmp);
  //       //   return grad;
  //       // }
  //       // else if (simplex.getDimension() == fesMesh.getDimension() - 1)
  //       // {
  //       //   assert(dynamic_cast<const Geometry::Face*>(&p.getSimplex()));
  //       //   const auto& face = static_cast<const Geometry::Face&>(p.getSimplex());
  //       //   mfem::FaceElementTransformations* ft =
  //       //     simplexMesh.getHandle().GetFaceElementTransformations(face.getIndex());
  //       //   if (simplexMesh.isSubMesh())
  //       //   {
  //       //     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(simplexMesh);
  //       //     assert(submesh.getParent() == fesMesh);
  //       //     if (ft->Elem1 && this->getTraceDomain().count(ft->Elem1->Attribute))
  //       //     {
  //       //       Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
  //       //       ft->Elem1->ElementNo = parentIdx;
  //       //       ft->Elem1No = parentIdx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (ft->Elem2 && this->getTraceDomain().count(ft->Elem2->Attribute))
  //       //     {
  //       //       Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem2No);
  //       //       ft->Elem2->ElementNo = parentIdx;
  //       //       ft->Elem2No = parentIdx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem2, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (face.isBoundary())
  //       //     {
  //       //       assert(ft->Elem1);
  //       //       Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
  //       //       ft->Elem1->ElementNo = parentIdx;
  //       //       ft->Elem1No = parentIdx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else
  //       //     {
  //       //       assert(false);
  //       //       grad.setConstant(NAN);
  //       //       return grad;
  //       //     }
  //       //   }
  //       //   else if (fesMesh.isSubMesh())
  //       //   {
  //       //     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(fesMesh);
  //       //     assert(submesh.getParent() == simplexMesh);
  //       //     const auto& s2pe = submesh.getElementMap();
  //       //     if (ft->Elem1 && s2pe.right.count(ft->Elem1No) &&
  //       //         this->getTraceDomain().count(ft->Elem1->Attribute))
  //       //     {
  //       //       Geometry::Index idx = s2pe.right.at(ft->Elem1No);
  //       //       ft->Elem1->ElementNo = idx;
  //       //       ft->Elem1No = idx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (ft->Elem2 && s2pe.right.count(ft->Elem2No) &&
  //       //         this->getTraceDomain().count(ft->Elem2->Attribute))
  //       //     {
  //       //       Geometry::Index idx = s2pe.right.at(ft->Elem2No);
  //       //       ft->Elem2->ElementNo = idx;
  //       //       ft->Elem2No = idx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem2, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (face.isBoundary())
  //       //     {
  //       //       assert(ft->Elem1);
  //       //       Geometry::Index idx = s2pe.right.at(ft->Elem1No);
  //       //       ft->Elem1->ElementNo = idx;
  //       //       ft->Elem1No = idx;
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else
  //       //     {
  //       //       assert(false);
  //       //       grad.setConstant(NAN);
  //       //       return grad;
  //       //     }
  //       //   }
  //       //   else
  //       //   {
  //       //     if (ft->Elem1 && this->getTraceDomain().count(ft->Elem1->Attribute))
  //       //     {
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (ft->Elem2 && this->getTraceDomain().count(ft->Elem2->Attribute))
  //       //     {
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem2, tmp);
  //       //       return grad;
  //       //     }
  //       //     else if (face.isBoundary())
  //       //     {
  //       //       ft->SetAllIntPoints(&p.getIntegrationPoint());
  //       //       assert(ft->Elem1);
  //       //       m_u.get().getHandle().GetGradient(*ft->Elem1, tmp);
  //       //       return grad;
  //       //     }
  //       //     else
  //       //     {
  //       //       assert(false);
  //       //       grad.setConstant(NAN);
  //       //       return grad;
  //       //     }
  //       //   }
  //       // }
  //       // else
  //       // {
  //       //   assert(false);
  //       //   grad.setConstant(NAN);
  //       //   return grad;
  //       // }
  //     }

  //     inline
  //     constexpr
  //     const Operand& getOperand() const
  //     {
  //       return m_u.get();
  //     }

  //     inline Grad* copy() const noexcept override
  //     {
  //       return new Grad(*this);
  //     }

  //   private:
  //     std::reference_wrapper<const Operand> m_u;
  // };

  // template <class ... Ts>
  // Grad(const GridFunction<H1<Ts...>>&) -> Grad<GridFunction<H1<Ts...>>>;

  /**
   * @ingroup GradSpecializations
   */
  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  class Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>, P1<Scalar, Ps...>, Space>
  {
    public:
      using FES = P1<Scalar, Ps...>;
      using Operand = ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>;
      using Parent = ShapeFunctionBase<Grad<Operand>, FES, Space>;

      Grad(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { getOperand().getFiniteElementSpace().getMesh().getSpaceDimension(), 1 };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      TensorBasis<Math::Vector> getTensorBasis(const Geometry::Point& p) const
      {
        using OperandRange =
          typename FormLanguage::Traits<ShapeFunctionBase<Operand, P1<Scalar, Ps...>, Space>>::RangeType;
        static_assert(std::is_same_v<OperandRange, Scalar>);
        const size_t d = p.getPolytope().getDimension();
        const Index i = p.getPolytope().getIndex();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(d, i);
        const Math::Vector& coords = p.getCoordinates(Geometry::Point::Coordinates::Reference);
        return (fe.getGradient(coords) * p.getJacobianInverse()).transpose();
      }

      inline Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class NestedDerived, class ... Ps, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>&)
    -> Grad<ShapeFunction<NestedDerived, P1<Scalar, Ps...>, Space>>;
}

#endif
