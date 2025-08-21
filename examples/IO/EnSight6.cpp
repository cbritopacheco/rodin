/*
 * Example usage of the EnSight6 I/O implementation
 * 
 * This example demonstrates how to use the EnSight6 format
 * for reading and writing mesh and grid function data.
 */

#include <Rodin/Geometry.h>
#include <Rodin/IO/EnSight6.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;

int main(int argc, char** argv)
{
  // Example 1: Load a mesh from EnSight6 format
  if (argc > 1)
  {
    Mesh<Context::Local> mesh;
    try
    {
      mesh.load(argv[1], FileFormat::ENSIGHT6);
      std::cout << "Successfully loaded mesh with:" << std::endl;
      std::cout << "  Vertices: " << mesh.getVertexCount() << std::endl;
      std::cout << "  Cells: " << mesh.getCellCount() << std::endl;
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error loading mesh: " << e.what() << std::endl;
      return 1;
    }
    
    // Example 2: Save the mesh back to EnSight6 format
    try
    {
      mesh.save("output.ens", FileFormat::ENSIGHT6);
      std::cout << "Successfully saved mesh to output.ens" << std::endl;
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error saving mesh: " << e.what() << std::endl;
      return 1;
    }
  }
  else
  {
    std::cout << "Usage: " << argv[0] << " <ensight6_file>" << std::endl;
    std::cout << "Example EnSight6 file format:" << std::endl;
    std::cout << R"(
EnSight6 Geometry File Format
Example mesh
node id given
element id given
coordinates
3
0  0.0  0.0  0.0
1  1.0  0.0  0.0
2  0.5  1.0  0.0
part 1
Triangle_Part
tria3
1
0 0 1 2
)" << std::endl;
  }
  
  return 0;
}