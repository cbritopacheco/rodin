#ifndef RODIN_VARIATIONAL_GAUSSIANQUADRATURE_H
#define RODIN_VARIATIONAL_GAUSSIANQUADRATURE_H

#include "Dot.h"
#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{

  /**
   * @defgroup GaussianQuadratureSpecializations GaussianQuadrature Template Specializations
   * @brief Template specializations of the GaussianQuadrature class.
   *
   * @see GaussianQuadrature
   */

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class GaussianQuadrature<
    Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Integrand = Dot<LHS, RHS>;
      using Parent = BilinearFormIntegratorBase;

      constexpr
      GaussianQuadrature(const LHS& lhs, const RHS& rhs)
        : GaussianQuadrature(Dot(lhs, rhs))
      {}

      constexpr
      GaussianQuadrature(const Integrand& prod)
        : BilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
          m_prod(prod.copy())
      {}

      constexpr
      GaussianQuadrature(const GaussianQuadrature& other)
        : BilinearFormIntegratorBase(other),
          m_prod(other.m_prod->copy())
      {}

      constexpr
      GaussianQuadrature(GaussianQuadrature&& other)
        : BilinearFormIntegratorBase(std::move(other)),
          m_prod(std::move(other.m_prod))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_prod);
        return *m_prod;
      }

      Math::Matrix getMatrix(const Geometry::Simplex& simplex) const final override
      {
        const auto& integrand = getIntegrand();
        const auto& trial = integrand.getLHS();
        const auto& test = integrand.getRHS();
        const auto& trans = simplex.getTransformation();
        const size_t order =
          trial.getFiniteElementSpace().getOrder(simplex) +
          test.getFiniteElementSpace().getOrder(simplex) +
          simplex.getTransformation().getHandle().OrderW();
        const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(static_cast<mfem::Geometry::Type>(simplex.getGeometry()), order);
        Math::Matrix res = Math::Matrix::Zero(test.getDOFs(simplex), trial.getDOFs(simplex));
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
          const auto& ip = ir->IntPoint(i);
          Geometry::Point p(simplex, trans, Internal::ip2vec(ir->IntPoint(i), simplex.getDimension()));
          res += ip.weight * p.getDistortion() * integrand.getMatrix(p);
        }
        return res;
      }

      virtual Region getRegion() const override = 0;

      virtual GaussianQuadrature* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_prod;
  };

  template <class NestedDerived, class FES>
  class GaussianQuadrature<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = LinearFormIntegratorBase;

      template <class LHSDerived, class RHSDerived>
      constexpr
      GaussianQuadrature(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : GaussianQuadrature(Dot(lhs, rhs))
      {}

      constexpr
      GaussianQuadrature(const Integrand& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      GaussianQuadrature(const GaussianQuadrature& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      GaussianQuadrature(GaussianQuadrature&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      Math::Vector getVector(const Geometry::Simplex& simplex) const final override
      {
        const auto& integrand = getIntegrand();
        assert(integrand.getRangeType() == RangeType::Scalar);
        const auto& trans = simplex.getTransformation();
        const size_t order =
          integrand.getFiniteElementSpace().getOrder(simplex) + trans.getHandle().OrderW();
        const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(static_cast<mfem::Geometry::Type>(simplex.getGeometry()), order);
        Math::Vector res = Math::Vector::Zero(integrand.getDOFs(simplex));
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
          const auto& ip = ir->IntPoint(i);
          Geometry::Point p(simplex, trans, Internal::ip2vec(ip, simplex.getDimension()));
          res += ip.weight * p.getDistortion() * integrand.getTensorBasis(p).getVector();
        }
        return res;
      }

      virtual Region getRegion() const override = 0;

      virtual GaussianQuadrature* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  /* ||-- OPTIMIZATIONS -----------------------------------------------------
   * GaussianQuadrature<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   * ---------------------------------------------------------------------->>
   */

  /**
   * @ingroup GaussianQuadratureSpecializations
   *
   * @f[
   * \int \nabla u \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class GaussianQuadrature<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;
      using LHS = ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>;
      using Integrand = Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>;

      constexpr
      GaussianQuadrature(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      GaussianQuadrature(const GaussianQuadrature& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      GaussianQuadrature(GaussianQuadrature&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      Math::Matrix getMatrix(const Geometry::Simplex& simplex) const override
      {
        const auto& fe = getIntegrand().getLHS()
                                       .getFiniteElementSpace()
                                       .getFiniteElement(simplex);
        Math::Matrix res(fe.getDOFs(), fe.getDOFs());
        mfem::DenseMatrix tmp(res.data(), res.rows(), res.cols());
        mfem::ConstantCoefficient one(1.0);
        mfem::DiffusionIntegrator bfi(one);
        bfi.AssembleElementMatrix(fe.getHandle(), simplex.getTransformation().getHandle(), tmp);
        return res;
      }

      virtual Region getRegion() const override = 0;

      virtual GaussianQuadrature* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  GaussianQuadrature(const Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>&)
    -> GaussianQuadrature<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, TrialFES, TrialSpace>>, TrialFES, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, TestFES, TestSpace>>, TestFES, TestSpace>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f u) \cdot v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>& fu,
  //         const ShapeFunction<FES, TestSpace>& v)
  //       : Integral(Dot(fu, v))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Math::Matrix getMatrix(
  //         const Geometry::Simplex& element) const override
  //     {
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //           getIntegrand().getLHS().getFiniteElementSpace(),
  //           getIntegrand().getRHS().getFiniteElementSpace(),
  //           element);

  //       if (&trial == &test)
  //       {
  //         mfem::DenseMatrix mat;
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             switch (getIntegrand().getLHS().getRHS().getRangeType())
  //             {
  //               case RangeType::Scalar:
  //               {
  //                 mfem::MassIntegrator bfi(q.template get<RangeType::Scalar>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Vector:
  //               {
  //                 mfem::VectorMassIntegrator bfi(q.template get<RangeType::Scalar>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Matrix:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //             }
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             switch (getIntegrand().getLHS().getRHS().getRangeType())
  //             {
  //               case RangeType::Scalar:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //               case RangeType::Vector:
  //               {
  //                 mfem::VectorMassIntegrator bfi(q.template get<RangeType::Matrix>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Matrix:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //             }
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
  //     const ShapeFunction<FES, TestSpace>&)
  //   -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f \nabla u) \cdot \nabla v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<
  //   Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>,
  //   Grad<ShapeFunction<FES, TestSpace>>>>
  //     : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>& fgu,
  //         const Grad<ShapeFunction<FES, TestSpace>>& gv)
  //       : Integral(Dot(fgu, gv))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual Math::Matrix getMatrix(
  //         const Geometry::Simplex& element) const override
  //     {
  //       mfem::DenseMatrix mat;
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //           getIntegrand().getLHS().getFiniteElementSpace(),
  //           getIntegrand().getRHS().getFiniteElementSpace(),
  //           element);
  //       if (&trial == &test)
  //       {
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             mfem::DiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             mfem::DiffusionIntegrator bfi(q.template get<RangeType::Matrix>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //       Math::Matrix res = Eigen::Map<Math::Matrix>(mat.GetData(), mat.NumRows(), mat.NumCols());
  //       return res;
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
  //     const Grad<ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f \nabla u) \cdot \nabla v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<
  //   Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>,
  //   Jacobian<ShapeFunction<FES, TestSpace>>>>
  //     : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>& fgu,
  //         const Jacobian<ShapeFunction<FES, TestSpace>>& gv)
  //       : Integral(Dot(fgu, gv))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual Math::Matrix getMatrix(const Geometry::Simplex& element) const override
  //     {
  //       mfem::DenseMatrix mat;
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //             getIntegrand().getLHS().getFiniteElementSpace(),
  //             getIntegrand().getRHS().getFiniteElementSpace(),
  //             element);
  //       if (&trial == &test)
  //       {
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             mfem::VectorDiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             assert(false); // Unimplemented
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //       Math::Matrix res = Eigen::Map<Math::Matrix>(mat.GetData(), mat.NumRows(), mat.NumCols());
  //       return res;
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>&,
  //     const Jacobian<ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  //  * ----------------------------------------------------------------------||
  //  */

  // /* ||-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<ShapeFunctionBase<TestSpace>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega f \cdot v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar, vector or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<ShapeFunctionBase<TestSpace>>
  // {
  //   public:
  //     using IntegrationOrder =
  //       std::function<int(const FiniteElementSpaceBase&, const Geometry::Simplex&)>;
  //     using Parent    = Integral<ShapeFunctionBase<TestSpace>>;
  //     using Integrand  = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;

  //     constexpr
  //     Integral(const FunctionBase& f, const ShapeFunction<FES, TestSpace>& v)
  //       : Integral(Dot(f, v))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& fes, const Geometry::Simplex& element)
  //           {
  //             return 2 * fes.getFiniteElement(element).GetOrder();
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     Integral& setIntegrationOrder(IntegrationOrder order)
  //     {
  //       return static_cast<Integral&>(Parent::setIntegrationOrder(order));
  //     }

  //     int getIntegrationOrder(
  //         const FiniteElementSpaceBase& fes, const Geometry::Simplex& element) const
  //     {
  //       return Parent::getIntegrationOrder(fes, element);
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Math::Vector getVector(const Geometry::Simplex& element) const override
  //     {
  //       const FunctionBase& f = getIntegrand().getLHS();

  //       const auto& fe = getIntegrand().getFiniteElementSpace()
  //                            .getFiniteElement(element);
  //       const mfem::IntegrationRule *ir =
  //         &mfem::IntRules.Get(
  //             fe.GetGeomType(),
  //             getIntegrationOrder(getIntegrand().getFiniteElementSpace(), element));
  //       auto q = f.build(element.getMesh());

  //       mfem::Vector vec;
  //       switch (f.getRangeType())
  //       {
  //         case RangeType::Scalar:
  //         {
  //           mfem::DomainLFIntegrator lfi(q.get<RangeType::Scalar>());
  //           lfi.SetIntRule(ir);
  //           lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
  //           Math::Vector res = Eigen::Map<Math::Vector>(vec.GetData(), vec.Size());
  //           return res;
  //         }
  //         case RangeType::Vector:
  //         {
  //           mfem::VectorDomainLFIntegrator lfi(q.get<RangeType::Vector>());
  //           lfi.SetIntRule(ir);
  //           lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
  //           Math::Vector res = Eigen::Map<Math::Vector>(vec.GetData(), vec.Size());
  //           return res;
  //         }
  //         case RangeType::Matrix:
  //         {
  //           assert(false); // Unsupported
  //           Math::Vector res = Eigen::Map<Math::Vector>(nullptr, 0);
  //           return res;
  //         }
  //       }
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
  //   -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // Integral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<ShapeFunctionBase<TestSpace>>
  //  * ----------------------------------------------------------------------||
  //  */
}

#endif
