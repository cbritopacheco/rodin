#include "Rodin/Geometry/Types.h"
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main(int argc, char** argv)
{

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <mesh-file>" << std::endl;
    return 1;
  }

  Mesh mesh;
  mesh.load(argv[1], IO::FileFormat::MEDIT);

  std::set<Geometry::Attribute> attrs;
  for (auto it = mesh.getCell(); it; ++it)
  {
    auto attr = it->getAttribute();
    if (attr)
      attrs.insert(*attr);
  }

  std::cout << "Cell attributes in mesh:" << std::endl;
  for (const auto& attr : attrs)
    std::cout << attr << ", ";

  attrs.clear();

  for (auto it = mesh.getFace(); it; ++it)
  {
    auto attr = it->getAttribute();
    if (attr)
      attrs.insert(*attr);
  }

  std::cout << "\nFace attributes in mesh:" << std::endl;
  for (const auto& attr : attrs)
    std::cout << attr << ", ";

  std::cout << std::endl;

  return 0;
}
