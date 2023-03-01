/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "H1.h"
#include "GridFunction.h"
#include "TensorBasis.h"
#include "ShapeFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JacobianSpecializations Jacobian Template Specializations
   * @brief Template specializations of the Jacobian class.
   * @see Jacobian
   */

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 GridFunction object.
   */
  template <class ... Ps>
  class Jacobian<GridFunction<H1<Math::Vector, Ps...>>> final
    : public MatrixFunctionBase<Jacobian<GridFunction<H1<Math::Vector, Ps...>>>>
  {
    public:
      using Operand = GridFunction<H1<Math::Vector, Ps...>>;
      using Parent = MatrixFunctionBase<Jacobian<GridFunction<H1<Math::Vector, Ps...>>>>;

      /**
       * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
       * @f$ u @f$.
       * @param[in] u Grid function to be differentiated
       */
      constexpr
      Jacobian(const Operand& u)
        :  m_u(u)
      {}

      constexpr
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      constexpr
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      inline
      constexpr
      size_t getRows() const
      {
        return m_u.get().getFiniteElementSpace().getVectorDimension();
      }

      inline
      constexpr
      size_t getColumns() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getDimension();
      }

      auto getValue(const Geometry::Point& p) const
      {
        const auto& simplex = p.getSimplex();
        const auto& simplexMesh = simplex.getMesh();
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        const size_t sdim = simplex.getMesh().getSpaceDimension();
        const size_t vdim = fes.getVectorDimension();
        Math::Matrix jacobian(vdim, sdim);
        mfem::DenseMatrix tmp(jacobian.data(), jacobian.rows(), jacobian.cols());
        if (simplex.getDimension() == fesMesh.getDimension())
        {
          const auto& trans = p.getTransformation();
          m_u.get().getHandle().GetVectorGradient(trans.getHandle(), tmp);
          return jacobian;
        }
        else if (simplex.getDimension() == fesMesh.getDimension() - 1)
        {
          assert(dynamic_cast<const Geometry::Face*>(&p.getSimplex()));
          const auto& face = static_cast<const Geometry::Face&>(p.getSimplex());
          mfem::FaceElementTransformations* ft =
            simplexMesh.getHandle().GetFaceElementTransformations(face.getIndex());
          if (simplexMesh.isSubMesh())
          {
            const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(simplexMesh);
            assert(submesh.getParent() == fesMesh);
            if (ft->Elem1 && this->getTraceDomain().count(ft->Elem1->Attribute))
            {
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
              ft->Elem1->ElementNo = parentIdx;
              ft->Elem1No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else if (ft->Elem2 && this->getTraceDomain().count(ft->Elem2->Attribute))
            {
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem2No);
              ft->Elem2->ElementNo = parentIdx;
              ft->Elem2No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem2, tmp);
              return jacobian;
            }
            else if (face.isBoundary())
            {
              assert(ft->Elem1);
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
              ft->Elem1->ElementNo = parentIdx;
              ft->Elem1No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else
            {
              assert(false);
              jacobian.setConstant(NAN);
              return jacobian;
            }
          }
          else if (fesMesh.isSubMesh())
          {
            const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(fesMesh);
            assert(submesh.getParent() == simplexMesh);
            const auto& s2pe = submesh.getElementMap();
            if (ft->Elem1 && s2pe.right.count(ft->Elem1No) &&
                this->getTraceDomain().count(ft->Elem1->Attribute))
            {
              Geometry::Index idx = s2pe.right.at(ft->Elem1No);
              ft->Elem1->ElementNo = idx;
              ft->Elem1No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else if (ft->Elem2 && s2pe.right.count(ft->Elem2No) &&
                this->getTraceDomain().count(ft->Elem2->Attribute))
            {
              Geometry::Index idx = s2pe.right.at(ft->Elem2No);
              ft->Elem2->ElementNo = idx;
              ft->Elem2No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem2, tmp);
              return jacobian;
            }
            else if (face.isBoundary())
            {
              assert(ft->Elem1);
              Geometry::Index idx = s2pe.right.at(ft->Elem1No);
              ft->Elem1->ElementNo = idx;
              ft->Elem1No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else
            {
              assert(false);
              jacobian.setConstant(NAN);
              return jacobian;
            }
          }
          else
          {
            if (ft->Elem1 && this->getTraceDomain().count(ft->Elem1->Attribute))
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else if (ft->Elem2 && this->getTraceDomain().count(ft->Elem2->Attribute))
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.get().getHandle().GetVectorGradient(*ft->Elem2, tmp);
              return jacobian;
            }
            else if (face.isBoundary())
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              assert(ft->Elem1);
              m_u.get().getHandle().GetVectorGradient(*ft->Elem1, tmp);
              return jacobian;
            }
            else
            {
              assert(false);
              jacobian.setConstant(NAN);
              return jacobian;
            }
          }
        }
        else
        {
          assert(false);
          jacobian.setConstant(NAN);
          return jacobian;
        }
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class ... Ps>
  Jacobian(const GridFunction<H1<Math::Vector, Ps...>>&)
    -> Jacobian<GridFunction<H1<Math::Vector, Ps...>>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 ShapeFunction object.
   */
  template <class ShapeFunctionDerived, ShapeFunctionSpaceType Space, class ... Ts>
  class Jacobian<ShapeFunction<ShapeFunctionDerived, H1<Math::Vector, Ts...>, Space>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<ShapeFunctionDerived, H1<Math::Vector, Ts...>, Space>>, H1<Math::Vector, Ts...>, Space>
  {
    public:
      using FES = H1<Math::Vector, Ts...>;
      using Operand = ShapeFunction<ShapeFunctionDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Jacobian<Operand>, FES, Space>;

      constexpr
      Jacobian(const Operand& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u)
      {}

      constexpr
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      constexpr
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(other.m_u)
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
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
        return { getOperand().getFiniteElementSpace().getMesh().getSpaceDimension(),
                 getOperand().getFiniteElementSpace().getVectorDimension() };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      TensorBasis<Math::Matrix> getTensorBasis(const Geometry::Point& p) const
      {
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(p.getSimplex());
        const auto& inv = p.getJacobianInverse();
        const Eigen::TensorMap<const Eigen::Tensor<Scalar, 2>> lift(inv.data(), inv.rows(), inv.cols());
        static constexpr const Eigen::array<Eigen::IndexPair<int>, 1> dims = { Eigen::IndexPair<int>(2, 0) };
        return fe.getJacobian(p.getReference()).contract(lift, dims)
                                               .shuffle(Eigen::array<int, 3>{2, 1, 0});
      }

      inline Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };

  template <class ShapeFunctionDerived, class FES, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, FES, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, FES, Space>>;
}

#endif
