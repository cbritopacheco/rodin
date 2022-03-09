#ifndef RODIN_VARIATIONAL_INNER_H
#define RODIN_VARIATIONAL_INNER_H

#include "Integral.h"

namespace Rodin::Variational
{
   template <>
   class Inner<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>
      : public Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
   {
      public:
         Inner(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>(
                  Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>(lhs, rhs))
         {}

         Inner(const Inner& other)
            : Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>(other)
         {}

         Inner(Inner&& other)
            : Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>(std::move(other))
         {}
   };
   Inner(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
      -> Inner<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

   template <>
   class Inner<ScalarCoefficientBase, ShapeFunctionBase<Test>>
      : public Integral<ShapeFunctionBase<Test>>
   {
      public:
         Inner(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral<ShapeFunctionBase<Test>>(
                  Dot<ScalarCoefficientBase, ShapeFunctionBase<Test>>(lhs, rhs))
         {
            assert(rhs.getFiniteElementSpace().getVectorDimension() == 1);
         }

         Inner(const Inner& other)
            : Integral<ShapeFunctionBase<Test>>(other)
         {}

         Inner(Inner&& other)
            : Integral<ShapeFunctionBase<Test>>(std::move(other))
         {}
   };
   Inner(const ScalarCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> Inner<ScalarCoefficientBase, ShapeFunctionBase<Test>>;

   template <>
   class Inner<VectorCoefficientBase, ShapeFunctionBase<Test>>
      : public Integral<ShapeFunctionBase<Test>>
   {
      public:
         Inner(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral<ShapeFunctionBase<Test>>(
                  Dot<VectorCoefficientBase, ShapeFunctionBase<Test>>(lhs, rhs))
         {}

         Inner(const Inner& other)
            : Integral<ShapeFunctionBase<Test>>(other)
         {}

         Inner(Inner&& other)
            : Integral<ShapeFunctionBase<Test>>(std::move(other))
         {}
   };
   Inner(const VectorCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> Inner<VectorCoefficientBase, ShapeFunctionBase<Test>>;

   template <>
   class BoundaryInner<ScalarCoefficientBase, ShapeFunctionBase<Test>>
      : public BoundaryIntegral<ShapeFunctionBase<Test>>
   {
      public:
         BoundaryInner(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(
                  Dot<ScalarCoefficientBase, ShapeFunctionBase<Test>>(lhs, rhs))
         {}

         BoundaryInner(const BoundaryInner& other)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(other)
         {}

         BoundaryInner(BoundaryInner&& other)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(std::move(other))
         {}
   };
   BoundaryInner(const ScalarCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> BoundaryInner<ScalarCoefficientBase, ShapeFunctionBase<Test>>;

   template <>
   class BoundaryInner<VectorCoefficientBase, ShapeFunctionBase<Test>>
      : public BoundaryIntegral<ShapeFunctionBase<Test>>
   {
      public:
         BoundaryInner(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(
                  Dot<VectorCoefficientBase, ShapeFunctionBase<Test>>(lhs, rhs))
         {}

         BoundaryInner(const BoundaryInner& other)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(other)
         {}

         BoundaryInner(BoundaryInner&& other)
            : BoundaryIntegral<ShapeFunctionBase<Test>>(std::move(other))
         {}
   };
   BoundaryInner(const VectorCoefficientBase&, const ShapeFunctionBase<Test>&)
      -> BoundaryInner<VectorCoefficientBase, ShapeFunctionBase<Test>>;
}

#endif
