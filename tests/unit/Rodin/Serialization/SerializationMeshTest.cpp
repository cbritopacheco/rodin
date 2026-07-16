#include <gtest/gtest.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "Rodin/Serialization/Export.h"
#include "Rodin/Serialization/Optional.h"

#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies to stringstream for serialization mesh by checking exact expected values, serialization round trips.
  TEST(Rodin_Serialization_Mesh, to_stringstream)
  {
    Mesh original;
    original = original.UniformGrid(Polytope::Type::Triangle, {16, 16});

    // Serialize the mesh into an in-memory string stream.
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      oa << original;
    }

    Rodin::Geometry::Mesh<Rodin::Context::Local> deserialized;
    {
      boost::archive::text_iarchive ia(ss);
      ia >> deserialized;
    }

    EXPECT_EQ(original.getDimension(), deserialized.getDimension());
    EXPECT_EQ(original.getSpaceDimension(), deserialized.getSpaceDimension());
    EXPECT_EQ(original.getVertexCount(), deserialized.getVertexCount());
    EXPECT_EQ(original.getCellCount(), deserialized.getCellCount());
  }

  /// @brief Verifies polymorphic serialization of a parametric transformation.
  TEST(Rodin_Serialization_Mesh, ParametricTransformationRoundTrip)
  {
    PointCloud points(2, 3);
    points(0, 0) = 0;
    points(0, 1) = 1;
    points(0, 2) = 0;
    points(1, 0) = 0;
    points(1, 1) = 0;
    points(1, 2) = 1;

    std::unique_ptr<PolytopeTransformation> original(
      new ParametricTransformation(points, RealP1Element(Polytope::Type::Triangle)));

    std::stringstream stream;
    {
      boost::archive::text_oarchive archive(stream);
      archive << original.get();
    }

    PolytopeTransformation* loadedRaw = nullptr;
    {
      boost::archive::text_iarchive archive(stream);
      archive >> loadedRaw;
    }
    std::unique_ptr<PolytopeTransformation> loaded(loadedRaw);

    Math::SpatialPoint reference(2);
    reference[0] = 0.25;
    reference[1] = 0.5;
    Math::SpatialPoint expected;
    Math::SpatialPoint actual;
    original->transform(expected, reference);
    loaded->transform(actual, reference);

    EXPECT_EQ(loaded->getReferenceDimension(), 2);
    EXPECT_EQ(loaded->getPhysicalDimension(), 2);
    EXPECT_EQ(loaded->getOrder(), 1);
    EXPECT_NEAR((actual - expected).norm(), 0, RODIN_FUZZY_CONSTANT);
  }
}
