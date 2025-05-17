#include <gtest/gtest.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Serialization_Mesh, to_stringstream)
  {
    Mesh original;
    original = original.UniformGrid(Polytope::Type::Triangle, { 16, 16 });

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
}

