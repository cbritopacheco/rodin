/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Concrete geometry-family reconstruction for native HDF5 mesh files.
#include "HDF5.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/Variational/P1/P1Element.h"

namespace Rodin::IO::HDF5
{
  namespace
  {
    static_assert(RODIN_IO_HDF5_MAXIMUM_GEOMETRY_ORDER > 0);

    /**
     * @brief Extracts one of the configured real H1 chart families.
     * @param transformation Chart being persisted.
     * @returns Family tag and retained control points.
     */
    template <size_t Degree = 1>
    TransformationData getH1TransformationData(
      const Geometry::PolytopeTransformation& transformation)
    {
      using Chart =
        Geometry::ParametricTransformation<Variational::RealH1Element<Degree>>;
      if (const auto* chart = dynamic_cast<const Chart*>(&transformation))
        return {static_cast<I32>(Degree + 1), &chart->getPointMatrix()};
      if constexpr (Degree < RODIN_IO_HDF5_MAXIMUM_GEOMETRY_ORDER)
        return getH1TransformationData<Degree + 1>(transformation);
      else
      {
        Alert::Exception() << "Unsupported HDF5 geometric transformation family or order."
                           << Alert::Raise;
        return {};
      }
    }

    /**
     * @brief Constructs a chart after validating its control-point count.
     * @param geometry Reference geometry.
     * @param points Physical control coordinates in the element's node order.
     * @returns A newly owned transformation.
     */
    template <class Element>
    std::unique_ptr<Geometry::PolytopeTransformation> makeChart(
      Geometry::Polytope::Type geometry, Geometry::PointCloud points)
    {
      Element element(geometry);
      if (points.cols() != element.getCount())
      {
        Alert::Exception()
          << "HDF5 geometry control-point count does not match its element."
          << Alert::Raise;
      }
      return std::make_unique<Geometry::ParametricTransformation<Element>>(
        std::move(points), std::move(element));
    }

    /**
     * @brief Restores the configured real H1 family selected by its persisted tag.
     * @param kind Chart-family tag.
     * @param geometry Reference geometry.
     * @param points Physical control coordinates.
     * @returns A newly owned transformation.
     */
    template <size_t Degree = 1>
    std::unique_ptr<Geometry::PolytopeTransformation> makeH1Transformation(
      I32 kind, Geometry::Polytope::Type geometry, Geometry::PointCloud points)
    {
      if (kind == static_cast<I32>(Degree + 1))
        return makeChart<Variational::RealH1Element<Degree>>(geometry, std::move(points));
      if constexpr (Degree < RODIN_IO_HDF5_MAXIMUM_GEOMETRY_ORDER)
        return makeH1Transformation<Degree + 1>(kind, geometry, std::move(points));
      else
      {
        Alert::Exception() << "Unsupported HDF5 geometric transformation tag: " << kind
                           << Alert::Raise;
        return {};
      }
    }
  }

  TransformationData getTransformationData(
    const Geometry::PolytopeTransformation& transformation)
  {
    using Chart = Geometry::ParametricTransformation<Variational::RealP1Element>;
    if (const auto* chart = dynamic_cast<const Chart*>(&transformation))
      return {1, &chart->getPointMatrix()};
    return getH1TransformationData(transformation);
  }

  std::unique_ptr<Geometry::PolytopeTransformation> makeTransformation(
    I32 kind, Geometry::Polytope::Type geometry, Geometry::PointCloud points)
  {
    if (kind == 1)
      return makeChart<Variational::RealP1Element>(geometry, std::move(points));
    return makeH1Transformation(kind, geometry, std::move(points));
  }
}
