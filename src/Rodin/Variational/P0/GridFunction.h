/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_GRIDFUNCTION_H
#define RODIN_VARIATIONAL_P0_GRIDFUNCTION_H

#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Geometry/SubMesh.h"

#include "P0.h"

// namespace Rodin::FormLanguage
// {
//   /**
//    * @ingroup TraitsSpecializations
//    */
//   template <class Range, class Mesh>
//   struct Traits<Variational::GridFunction<Variational::P0<Range, Mesh>>>
//   {
//     using FESType = Variational::P0<Range, Mesh>;
//     using MeshType = typename FormLanguage::Traits<FESType>::MeshType;
//     using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
//     using ContextType = typename FormLanguage::Traits<FESType>::ContextType;
//     using ElementType = typename FormLanguage::Traits<FESType>::ElementType;
//   };
// }
// 
// namespace Rodin::Variational
// {
//   /**
//    * @ingroup GridFunctionSpecializations
//    * @brief P0 GridFunction
//    */
//   template <class Range, class Mesh>
//   class GridFunction<P0<Range, Mesh>> final
//     : public GridFunctionBase<P0<Range, Mesh>, GridFunction<P0<Range, Mesh>>>
//   {
//     public:
//       /// Type of finite element space to which the GridFunction belongs to
//       using FESType = P0<Range, Mesh>;
// 
//       using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
// 
//       using MeshType = typename FormLanguage::Traits<FESType>::MeshType;
// 
//       using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
// 
//       using ContextType = typename FormLanguage::Traits<FESType>::ContextType;
// 
//       using ElementType = typename FormLanguage::Traits<FESType>::ElementType;
// 
//       /// Parent class
//       using Parent = GridFunctionBase<FESType, GridFunction<FESType>>;
// 
//       using Parent::getValue;
//       using Parent::operator=;
//       using Parent::operator+=;
//       using Parent::operator-=;
//       using Parent::operator*=;
//       using Parent::operator/=;
// 
//       /**
//        * @brief Constructs a grid function on a finite element space.
//        * @param[in] fes Finite element space to which the function belongs
//        * to.
//        */
//       GridFunction(const FESType& fes)
//         : Parent(fes)
//       {}
// 
//       /**
//        * @brief Copies the grid function.
//        * @param[in] other Other grid function to copy.
//        */
//       GridFunction(const GridFunction& other)
//         : Parent(other)
//       {}
// 
//       /**
//        * @brief Move constructs the grid function.
//        * @param[in] other Other grid function to move.
//        */
//       GridFunction(GridFunction&& other)
//         : Parent(std::move(other))
//       {}
// 
//       /**
//        * @brief Move assignment operator.
//        */
//       constexpr
//       GridFunction& operator=(GridFunction&& other)
//       {
//         Parent::operator=(std::move(other));
//         return *this;
//       }
// 
//       GridFunction& operator=(const GridFunction&)  = delete;
// 
//       void interpolate(RangeType& res, const Geometry::Point& p) const
//       {
//         static_assert(std::is_same_v<RangeType, Real>);
//         const auto& fes = this->getFiniteElementSpace();
//         const auto& mesh = fes.getMesh();
//         const auto& polytope = p.getPolytope();
//         assert(mesh == polytope.getMesh());
//         const size_t d = mesh.getDimension();
//         assert(d == polytope.getDimension());
//         const Index i = polytope.getIndex();
//         assert(fes.getFiniteElement(d, i).getCount() == 1);
//         if constexpr (std::is_same_v<RangeType, ScalarType>)
//         {
//           res = getValue({ d, i }, 0);
//         }
//         else if constexpr (std::is_same_v<RangeType, Math::Vector<ScalarType>>)
//         {
//           assert(false);
//         }
//         else
//         {
//           assert(false);
//         }
//       }
// 
//       GridFunction& setWeights()
//       {
//         auto& data = this->getData();
//         auto& w = this->getWeights().emplace(this->getFiniteElementSpace().getSize());
//         if (!this->getWeights().has_value())
//         {
//           auto& weights = this->getWeights().emplace(this->getFiniteElementSpace().getSize());
//           if constexpr (std::is_same_v<RangeType, Real>)
//           {
//             assert(data.rows() == 1);
//             w = data.transpose();
//           }
//           else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
//           {
//             assert(false);
//             weights.setConstant(NAN);
//           }
//           else
//           {
//             assert(false);
//             weights.setConstant(NAN);
//           }
//         }
//         return *this;
//       }
// 
//       template <class Vector>
//       GridFunction& setWeights(Vector&& weights)
//       {
//         auto& data = this->getData();
//         auto& w = this->getWeights().emplace();
//         Math::duplicate(weights, w);
//         Math::copy(weights, w);
//         if constexpr (std::is_same_v<RangeType, Real>)
//         {
//           assert(data.rows() == 1);
//           data = w.transpose();
//         }
//         else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
//         {
//           assert(false);
//           data.setConstant(NAN);
//         }
//         else
//         {
//           assert(false);
//           data.setConstant(NAN);
//         }
//         return *this;
//       }
// 
//   };
// 
//   template <class Range, class Mesh>
//   GridFunction(const P0<Range, Mesh>&) -> GridFunction<P0<Range, Mesh>>;
// }

#endif

