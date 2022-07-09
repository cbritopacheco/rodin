/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "MeshPrinter.h"

namespace Rodin::IO
{
   Status MeshPrinter<MeshFormat::GMSH, Traits::Serial>::print(std::ostream& os)
   {
      const auto& mfemMesh = getObject().getHandle();

      // 1. Head of the file
      //====================
      os << "$MeshFormat\n" "2.2 0 " + std::to_string(os.precision()) + "\n$EndMeshFormat\n";

      // 2. Writing nodes
      //================
      os << "$Nodes\n";
      os << mfemMesh.GetNV() << "\n";
      for (int i = 0; i < mfemMesh.GetNV(); i++)
      {
         os << i + 1
            << " " << mfemMesh.GetVertex(i)[0]
            << " " << mfemMesh.GetVertex(i)[1] << " " << mfemMesh.GetVertex(i)[2]
            << "\n"; 
      }
      os << "$EndNodes\n";

      // 3. Writing elements
      //====================
      int num_ele_bnd = mfemMesh.GetNBE(), num_ele = mfemMesh.GetNE();
      os << "$Elements\n";
      os << (num_ele_bnd + num_ele) << "\n";

      // 3.1. Writing boundary elements
      //===============================
      for (int i = 0; i < num_ele_bnd; i++)
      {
         int ele_type = 0;
         switch (mfemMesh.GetBdrElement(i)->GetType())
         {
            case mfem::Element::SEGMENT:
            {
               ele_type = 1;
               break;
            }
            case mfem::Element::TRIANGLE:
            {
               ele_type = 2;
               break;
            }
            case mfem::Element::QUADRILATERAL:
            {
               ele_type = 3;
               break;
            }
            default:
            {
               Alert::Exception("Unknown boundary element type.").raise();
            }
         }
         os << i + 1
            << " " << ele_type
            << " " << 1
            << " " << mfemMesh.GetBdrElement(i)->GetAttribute();
         const int* vertices_list = mfemMesh.GetBdrElement(i)->GetVertices();
         for (int j = 0; j < (mfemMesh.GetBdrElement(i)->GetNVertices()); j++)
            os << " " << (vertices_list[j] + 1);
         os << "\n"; 
      }

      // 3.2. Writing bulk elements
      //===========================
      for (int i = 0; i < num_ele; i++)
      {
         int ele_type = 0;
         switch (mfemMesh.GetElement(i)->GetType())
         {
            case mfem::Element::SEGMENT:
            {
               ele_type = 1;
               break;
            }
            case mfem::Element::TRIANGLE:
            {
               ele_type = 2;
               break;
            }
            case mfem::Element::QUADRILATERAL:
            {
               ele_type = 3;
               break;
            }
            case mfem::Element::TETRAHEDRON:
            {
               ele_type = 4;
               break;
            }
            case mfem::Element::HEXAHEDRON:
            {
               ele_type = 5;
               break;
            }
            case mfem::Element::PYRAMID:
            {
               ele_type = 7;
               break;
            }
            default:
            {
               Alert::Exception("Unknown element type.").raise();
            }
         }

         os << (i + 1 + num_ele_bnd)
            << " " << ele_type
            << " " << 1 << " " << mfemMesh.GetElement(i)->GetAttribute();

         const int* vertices_list = mfemMesh.GetElement(i)->GetVertices();
         for (int j = 0; j < mfemMesh.GetElement(i)->GetNVertices(); j++)
            os << " " << (vertices_list[j] + 1);
         os << "\n"; 
      }
      os << "$EndElements\n";

      return {true, {}};
   }

   IO::Status MeshPrinter<MeshFormat::MEDIT, Traits::Serial>::print(std::ostream& os)
   {
      const auto& mfemMesh = getObject().getHandle();

      int meshDim = mfemMesh.Dimension();
      int spaceDim = mfemMesh.SpaceDimension();

      os << "MeshVersionFormatted 2"
         << '\n'
         << "Dimension " << spaceDim
         << "\n\n";

      // Print vertices
      os << "Vertices"
         << '\n'
         << mfemMesh.GetNV()
         << '\n';
      for (int i = 0; i < mfemMesh.GetNV(); i++)
      {
         const double* coords = mfemMesh.GetVertex(i);
         for (int j = 0; j < mfemMesh.SpaceDimension(); j++)
            os << coords[j] << " ";
         os << 0 << '\n'; // Dummy reference for all vertices
      }

      // Print triangles
      os << '\n' << "Triangles" << '\n';
      switch (mfemMesh.Dimension())
      {
         case 2:
         {
            os << mfemMesh.GetNE() << '\n';
            mfem::Array<int> v;
            for (int i = 0; i < mfemMesh.GetNE(); i++)
            {
               mfemMesh.GetElementVertices(i, v);
               os << v[0] + 1 << " "
                  << v[1] + 1 << " "
                  << v[2] + 1 << " "
                  << mfemMesh.GetAttribute(i) << '\n';
            }
            break;
         }
         case 3:
         {
            os << mfemMesh.GetNBE() << '\n';
            mfem::Array<int> v;
            for (int i = 0; i < mfemMesh.GetNBE(); i++)
            {
               mfemMesh.GetBdrElementVertices(i, v);
               os << v[0] + 1 << " "
                  << v[1] + 1 << " "
                  << v[2] + 1 << " "
                  << mfemMesh.GetBdrAttribute(i) << '\n';
            }
            break;
         }
         default:
         {
            return {false, IO::Error{"Bad mesh dimension: " + std::to_string(meshDim)}};
         }
      }

      // Print tetrahedra
      switch (meshDim)
      {
         case 3:
         {
            os << '\n' << "Tetrahedra" << '\n'
               << mfemMesh.GetNE() << '\n';
            mfem::Array<int> v;
            for (int i = 0; i < mfemMesh.GetNE(); i++)
            {
               mfemMesh.GetElementVertices(i, v);
               os << v[0] + 1 << " "
                  << v[1] + 1 << " "
                  << v[2] + 1 << " "
                  << v[3] + 1 << " "
                  << mfemMesh.GetAttribute(i) << '\n';
            }
         }
         default:
         {
            // Do nothing
         }
      }

      // Print edges
      switch (meshDim)
      {
         case 2:
         {
            os << '\n' << "Edges" << '\n'
               << mfemMesh.GetNBE() << '\n';
            mfem::Array<int> v;
            for (int i = 0; i < mfemMesh.GetNBE(); i++)
            {
               mfemMesh.GetBdrElementVertices(i, v);
               os << v[0] + 1 << " "
                  << v[1] + 1 << " "
                  << mfemMesh.GetBdrAttribute(i) << '\n';
            }
         }
         default:
         {
            // Do nothing
         }
      }

      return {true, {}};
   }
}
