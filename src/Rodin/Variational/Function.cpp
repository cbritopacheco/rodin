#include "Rodin/Alert.h"
#include "Rodin/Geometry/SubMesh.h"

#include "Function.h"
#include "Transpose.h"

namespace Rodin::Variational
{
  Transpose<FunctionBase> FunctionBase::T() const
  {
    return Transpose<FunctionBase>(*this);
  }

  RangeType FunctionBase::getRangeType() const
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
      return RangeType::Matrix;
    }
  }

  Internal::MFEMFunction FunctionBase::build(const Geometry::MeshBase& mesh) const
  {
    switch (getRangeType())
    {
      case RangeType::Scalar:
        return Internal::MFEMFunction(new Internal::ScalarProxyFunction(mesh, *this));
      case RangeType::Vector:
        return Internal::MFEMFunction(new Internal::VectorProxyFunction(mesh, *this));
      case RangeType::Matrix:
        return Internal::MFEMFunction(new Internal::MatrixProxyFunction(mesh, *this));
    }
    // The following return is needed to prevent compiler warnings/errors
    return Internal::MFEMFunction(nullptr);
  }
}

namespace Rodin::Variational::Internal
{}

