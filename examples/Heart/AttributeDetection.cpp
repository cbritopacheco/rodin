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
    attrs.insert(it->getAttribute());

  std::cout << "Cell attributes in mesh:" << std::endl;
  for (const auto& attr : attrs)
    std::cout << attr << std::endl;

  attrs.clear();

  for (auto it = mesh.getFace(); it; ++it)
    attrs.insert(it->getAttribute());

  std::cout << "\nFace attributes in mesh:" << std::endl;
  for (const auto& attr : attrs)
    std::cout << attr << std::endl;

  return 0;
}
