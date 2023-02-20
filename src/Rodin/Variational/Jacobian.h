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
#include "BasisOperator.h"
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
  template <class ... Ts>
  class Jacobian<GridFunction<H1<Ts...>>> final
    : public MatrixFunctionBase<Jacobian<GridFunction<H1<Ts...>>>>
  {
    public:
      using Operand = GridFunction<H1<Ts...>>;
      using Parent = MatrixFunctionBase<Jacobian<Operand>>;

      /**
       * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
       * @f$ u @f$.
       * @param[in] u Grid function to be differentiated
       */
      constexpr
      Jacobian(Operand& u)
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
        return m_u.getFiniteElementSpace().getMesh().getDimension();
      }

      inline
      constexpr
      size_t getColumns() const
      {
        return m_u.getFiniteElementSpace().getVectorDimension();
      }

      auto getValue(const Geometry::Point& p) const
      {
        mfem::DenseMatrix grad;
        const auto& simplex = p.getSimplex();
        const auto& simplexMesh = simplex.getMesh();
        const auto& fesMesh = m_u.getFiniteElementSpace().getMesh();
        if (simplex.getDimension() == fesMesh.getDimension())
        {
          assert(dynamic_cast<const Geometry::Element*>(&p.getSimplex()));
          const auto& element = p.getSimplex();
          auto& trans = element.getTransformation();
          m_u.getHandle().GetVectorGradient(trans, grad);
          Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
          return res;
        }
        else if (simplex.getDimension() == fesMesh.getDimension() - 1)
        {
          assert(dynamic_cast<const Geometry::Face*>(&p.getSimplex()));
          const auto& face = static_cast<const Geometry::Face&>(p.getSimplex());
          mfem::FaceElementTransformations* ft =
            const_cast<Geometry::MeshBase&>(simplexMesh).getHandle()
            .GetFaceElementTransformations(face.getIndex());
          if (simplexMesh.isSubMesh())
          {
            const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(simplexMesh);
            assert(submesh.getParent() == fesMesh);
            if (ft->Elem1 && this->getTraceDomain() == ft->Elem1->Attribute)
            {
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
              ft->Elem1->ElementNo = parentIdx;
              ft->Elem1No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (ft->Elem2 && this->getTraceDomain() == ft->Elem2->Attribute)
            {
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem2No);
              ft->Elem2->ElementNo = parentIdx;
              ft->Elem2No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (face.isBoundary())
            {
              assert(ft->Elem1);
              Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
              ft->Elem1->ElementNo = parentIdx;
              ft->Elem1No = parentIdx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else
            {
              assert(false);
              return Math::Matrix(0, 0);
            }
          }
          else if (fesMesh.isSubMesh())
          {
            const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(fesMesh);
            assert(submesh.getParent() == simplexMesh);
            const auto& s2pe = submesh.getElementMap();
            if (ft->Elem1 && s2pe.right.count(ft->Elem1No) &&
                this->getTraceDomain() == ft->Elem1->Attribute)
            {
              Geometry::Index idx = s2pe.right.at(ft->Elem1No);
              ft->Elem1->ElementNo = idx;
              ft->Elem1No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (ft->Elem2 && s2pe.right.count(ft->Elem2No) &&
                this->getTraceDomain() == ft->Elem2->Attribute)
            {
              Geometry::Index idx = s2pe.right.at(ft->Elem2No);
              ft->Elem2->ElementNo = idx;
              ft->Elem2No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (face.isBoundary())
            {
              assert(ft->Elem1);
              Geometry::Index idx = s2pe.right.at(ft->Elem1No);
              ft->Elem1->ElementNo = idx;
              ft->Elem1No = idx;
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else
            {
              assert(false);
              return Math::Matrix(0, 0);
            }
          }
          else
          {
            if (ft->Elem1 && this->getTraceDomain() == ft->Elem1->Attribute)
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (ft->Elem2 && this->getTraceDomain() == ft->Elem2->Attribute)
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else if (face.isBoundary())
            {
              ft->SetAllIntPoints(&p.getIntegrationPoint());
              assert(ft->Elem1);
              m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
              Math::Matrix res = Eigen::Map<Math::Matrix>(grad.GetData(), grad.NumRows(), grad.NumCols());
              return res;
            }
            else
            {
              assert(false);
              return Math::Matrix(0, 0);
            }
          }
        }
        else
        {
          assert(false);
          return Math::Matrix(0, 0);
        }
      }

    private:
      std::reference_wrapper<Operand> m_u;
  };

  template <class FES>
  Jacobian(GridFunction<FES>&) -> Jacobian<GridFunction<FES>>;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of an H1 ShapeFunction object.
   */
  template <class NestedDerived, ShapeFunctionSpaceType Space, class ... Ts>
  class Jacobian<ShapeFunction<NestedDerived, H1<Ts...>, Space>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<NestedDerived, H1<Ts...>, Space>>, Space>
  {
    public:
      using FES = H1<Ts...>;
      using Operand = ShapeFunction<NestedDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Jacobian<Operand>, Space>;

      constexpr
      Jacobian(Operand& u)
        : m_u(u)
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
      FES& getFiniteElementSpace()
      {
        return m_u.getFiniteElementSpace();
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_u.getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return m_u.get().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return { m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension(),
                 m_u.getFiniteElementSpace().getVectorDimension() };
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return m_u.get().getDOFs(element);
      }

      auto getOperator(ShapeComputator& compute, const Geometry::Point& p) const
      {
        const auto& element = p.getSimplex();
        auto& trans = element.getTransformation();
        const auto& fe = getFiniteElementSpace().getFiniteElement(element);
        const auto& dshape = compute.getPhysicalDShape(fe, trans, trans.GetIntPoint());
        const size_t n = dshape.NumRows();
        const size_t sdim = trans.GetSpaceDim();
        const size_t vdim = m_u.getFiniteElementSpace().getVectorDimension();
        assert(false);
        return TensorBasis(vdim * n,
            [sdim, vdim](size_t) -> Math::Matrix { return Math::Matrix::Zero(sdim, vdim); });
        // TensorBasis res(static_cast<int>(vdim * n), static_cast<int>(sdim), static_cast<int>(vdim));
        // res.setZero();
        // for (size_t i = 0; i < vdim; i++)
        // {
        //   for (size_t j = 0; j < n; j++)
        //   {
        //     for (size_t k = 0; k < sdim; k++)
        //     {
        //       res(static_cast<int>(j + i * n), static_cast<int>(k), static_cast<int>(i)) = dshape(j, k);
        //     }
        //   }
        // }
      }

    private:
      std::reference_wrapper<Operand> m_u;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Jacobian(ShapeFunction<NestedDerived, FES, Space>&)
    -> Jacobian<ShapeFunction<NestedDerived, FES, Space>>;
}

#endif
