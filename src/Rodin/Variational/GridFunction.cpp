/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GridFunction.h"

#include "Rodin/IO/GridFunctionLoader.h"

#include "Exceptions.h"

#include "Sum.h"
#include "Mult.h"
#include "Minus.h"
#include "Division.h"
#include "Component.h"
#include "BooleanFunction.h"
#include "FiniteElementSpace.h"


namespace Rodin::Variational
{
  GridFunctionBase& GridFunctionBase::update()
  {
    getHandle().Update();
    return *this;
  }

  double GridFunctionBase::max() const
  {
    return getHandle().Max();
  }

  double GridFunctionBase::min() const
  {
    return getHandle().Min();
  }

  std::pair<const double*, int> GridFunctionBase::getData() const
  {
    return {static_cast<const double*>(getHandle().GetData()), getHandle().Size()};
  }

  GridFunctionBase& GridFunctionBase::setData(std::unique_ptr<double[]> data, int size)
  {
    getHandle().SetDataAndSize(data.release(), size);
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator+=(double t)
  {
    getHandle() += t;
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator+=(const GridFunctionBase& rhs)
  {
    if (this == &rhs)
    {
      operator*=(2.0);
    }
    else
    {
      assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
      getHandle() += rhs.getHandle();
    }
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator-=(double t)
  {
    getHandle() -= t;
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator-=(const GridFunctionBase& rhs)
  {
    if (this == &rhs)
    {
      operator=(0.0);
    }
    else
    {
      assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
      getHandle() -= rhs.getHandle();
    }
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator*=(double t)
  {
    getHandle() *= t;
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator*=(const GridFunctionBase& rhs)
  {
    if (this == &rhs)
    {
      for (auto& v : getHandle())
        v *= v;
    }
    else
    {
      assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
      getHandle() *= rhs.getHandle();
    }
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator/=(double t)
  {
    getHandle() /= t;
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator/=(const GridFunctionBase& rhs)
  {
    if (this == &rhs)
    {
      operator=(1.0);
    }
    else
    {
      assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
      getHandle() /= rhs.getHandle();
    }
    return *this;
  }

  GridFunctionBase& GridFunctionBase::operator=(double v)
  {
    getHandle() = v;
    return *this;
  }

  GridFunctionBase& GridFunctionBase::project(
      const FunctionBase& s, const std::set<int>& attrs)
  {
    auto va = s.build(getFiniteElementSpace().getMesh());
    switch (s.getRangeType())
    {
      case RangeType::Scalar:
      {
        if (attrs.size() == 0)
          getHandle().ProjectCoefficient(va.get<RangeType::Scalar>());
        else
        {
          mfem::Array<int> vdofs;
          const auto& fes = getFiniteElementSpace().getHandle();
          for (int i = 0; i < fes.GetNE(); i++)
          {
            if (attrs.count(fes.GetAttribute(i)) > 0)
            {
              fes.GetElementVDofs(i, vdofs);
              getHandle().ProjectCoefficient(va.get<RangeType::Scalar>(), vdofs);
            }
          }
        }
        break;
      }
      case RangeType::Vector:
      {
        if (attrs.size() == 0)
          getHandle().ProjectCoefficient(va.get<RangeType::Vector>());
        else
        {
          mfem::Array<int> vdofs;
          const auto& fes = getFiniteElementSpace().getHandle();
          for (int i = 0; i < fes.GetNE(); i++)
          {
            if (attrs.count(fes.GetAttribute(i)) > 0)
            {
              fes.GetElementVDofs(i, vdofs);
              getHandle().ProjectCoefficient(va.get<RangeType::Vector>(), vdofs);
            }
          }
        }
        break;
      }
      case RangeType::Matrix:
      {
        UnexpectedRangeTypeException(
            {RangeType::Scalar, RangeType::Vector}, RangeType::Matrix).raise();
        break;
      }
    }
    return *this;
  }

  GridFunctionBase& GridFunctionBase::project(const Restriction<FunctionBase>& s)
  {
    assert(false); // Not implemented
    // assert(getFiniteElementSpace().getVectorDimension() == 1);
    // auto iv = s.getScalarFunction().build();
    // getHandle() = std::numeric_limits<double>::quiet_NaN();
    // mfem::Array<int> vdofs;
    // mfem::Vector vals;
    // const auto& fes = getFiniteElementSpace().getHandle();
    // const auto& attrs = s.getAttributes();
    // for (int i = 0; i < fes.GetNE(); i++)
    // {
    //   if (attrs.count(fes.GetAttribute(i)) > 0)
    //   {
    //     fes.GetElementVDofs(i, vdofs);
    //     vals.SetSize(vdofs.Size());
    //     fes.GetFE(i)->Project(
    //         *iv, *fes.GetElementTransformation(i), vals);
    //     getHandle().SetSubVector(vdofs, vals);
    //   }
    // }
    return *this;
  }

  FunctionValue GridFunctionBase::getValue(const Geometry::Point& p) const
  {
    auto& trans = p.getSimplex().getTransformation();
    switch (getRangeType())
    {
      case RangeType::Scalar:
      {
        return getHandle().GetValue(trans, trans.GetIntPoint());
      }
      case RangeType::Vector:
      {
        FunctionValue::Vector value;
        getHandle().GetVectorValue(trans, trans.GetIntPoint(), value);
        return value;
      }
      case RangeType::Matrix:
      {
        assert(false);
        return 0.0;
      }
      default:
      {
        assert(false);
        return NAN;
      }
    }
    return NAN;
  }

  RangeShape GridFunctionBase::getRangeShape() const
  {
    return {getFiniteElementSpace().getVectorDimension(), 1};
  }

  Component<FunctionBase> GridFunctionBase::x() const
  {
    assert(getFiniteElementSpace().getVectorDimension() >= 1);
    return Component(*this, 0);
  }

  Component<FunctionBase> GridFunctionBase::y() const
  {
    assert(getFiniteElementSpace().getVectorDimension() >= 2);
    return Component(*this, 1);
  }

  Component<FunctionBase> GridFunctionBase::z() const
  {
    assert(getFiniteElementSpace().getVectorDimension() >= 3);
    return Component(*this, 2);
  }

  template <>
  GridFunction<L2<Context::Serial>>&
  GridFunction<L2<Context::Serial>>::load(
      const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    mfem::named_ifgzstream input(filename.c_str());
    if (!input)
    {
      Alert::Exception()
        << "Failed to open " << filename << " for reading."
        << Alert::Raise;
    }

    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::GridFunctionLoader<IO::FileFormat::MFEM, L2<Context::Serial>> loader(*this);
        loader.load(input);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::GridFunctionLoader<IO::FileFormat::MEDIT, L2<Context::Serial>> loader(*this);
        loader.load(input);
        break;
      }
      default:
      {
        Alert::Exception()
          << "Loading from \"" << fmt << "\" format unssuported."
          << Alert::Raise;
      }
    }
    return *this;
  }

  template <>
  void GridFunction<L2<Context::Serial>>
  ::save(const boost::filesystem::path& filename, IO::FileFormat fmt, int precision)
  const
  {
    std::ofstream output(filename.c_str());
    if (!output)
    {
      Alert::Exception()
        << "Failed to open " << filename << " for writing."
        << Alert::Raise;
    }

    output.precision(precision);
    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::GridFunctionPrinter<IO::FileFormat::MFEM, L2<Context::Serial>> printer(*this);
        printer.print(output);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::GridFunctionPrinter<IO::FileFormat::MEDIT, L2<Context::Serial>> printer(*this);
        printer.print(output);
        break;
      }
      default:
      {
        Alert::Exception()
          << "Saving to \"" << fmt << "\" format unssuported."
          << Alert::Raise;
      }
    }
  }

  template <>
  GridFunction<H1<Context::Serial>>&
  GridFunction<H1<Context::Serial>>::load(
      const boost::filesystem::path& filename, IO::FileFormat fmt)
  {
    mfem::named_ifgzstream input(filename.c_str());
    if (!input)
    {
      Alert::Exception()
        << "Failed to open " << filename << " for reading."
        << Alert::Raise;
    }

    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::GridFunctionLoader<IO::FileFormat::MFEM, H1<Context::Serial>> loader(*this);
        loader.load(input);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::GridFunctionLoader<IO::FileFormat::MEDIT, H1<Context::Serial>> loader(*this);
        loader.load(input);
        break;
      }
      default:
      {
        Alert::Exception()
          << "Loading from \"" << fmt << "\" format unssuported."
          << Alert::Raise;
      }
    }
    return *this;
  }

  template <>
  void GridFunction<H1<Context::Serial>>
  ::save(const boost::filesystem::path& filename, IO::FileFormat fmt, int precision)
  const
  {
    std::ofstream output(filename.c_str());
    if (!output)
    {
      Alert::Exception()
        << "Failed to open " << filename << " for writing."
        << Alert::Raise;
    }

    output.precision(precision);
    switch (fmt)
    {
      case IO::FileFormat::MFEM:
      {
        IO::GridFunctionPrinter<IO::FileFormat::MFEM, H1<Context::Serial>> printer(*this);
        printer.print(output);
        break;
      }
      case IO::FileFormat::MEDIT:
      {
        IO::GridFunctionPrinter<IO::FileFormat::MEDIT, H1<Context::Serial>> printer(*this);
        printer.print(output);
        break;
      }
      default:
      {
        Alert::Exception()
          << "Saving to \"" << fmt << "\" format unssuported."
          << Alert::Raise;
      }
    }
  }

  std::set<Geometry::Point> GridFunctionBase::where(
      const BooleanFunctionBase& pred,
      const std::set<int>& attrs,
      std::function<int(mfem::ElementTransformation&)> order) const
  {
    // std::set<Geometry::Point> result;
    // const auto& fes = getFiniteElementSpace();
    // const auto& mesh = fes.getMesh();
    assert(false);
    // for (int i = 0; i < mesh.count<Geometry::Element>(); i++)
    // {
    //   if (attrs.size() == 0 || attrs.count(mesh.get<Geometry::Element>(i).getAttribute()))
    //   {
    //     const auto& element = mesh.get<Geometry::Element>(i);
    //     auto& trans = element.getTransformation();
    //     for (const auto& p : element.getIntegrationRule(order(trans)))
    //     {
    //       const auto& ip = trans.GetIntPoint();
    //       trans.SetIntPoint(&ip);
    //       if (pred.getValue(p)) result.insert(p);
    //     }
    //   }
    // }
    // return result;
    return {};
  }

  int GridFunctionBase::getDimension() const
  {
    return getFiniteElementSpace().getVectorDimension();
  }


  RangeType GridFunctionBase::getRangeType() const
  {
    auto shape = getRangeShape();
    if (shape.height() == 1 && shape.width() == 1)
    {
      return RangeType::Scalar;
    }
    else if (shape.height() > 1 && shape.width() == 1)
    {
      return RangeType::Vector;
    }
    else
    {
      assert(false);
      return RangeType::Matrix;
    }
  }

  VectorFunctionBase* GridFunctionBase::copy() const noexcept
  {
    return new Internal::GridFunctionEvaluator(*this);
  }
}

