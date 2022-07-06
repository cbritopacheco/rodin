/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Printer.h"

namespace Rodin::MeshTools
{
   PrinterBase::Status Printer<MeshFormat::GMSH>::print(std::ostream& os)
   {
      // 1. Head of the file
      //====================
      os << "$MeshFormat\n" "2.2 0 " + std::to_string(os.precision()) + "\n$EndMeshFormat\n";

      // 2. Writing nodes
      //================
      os << "$Nodes\n";
      os << getMesh().GetNV() << "\n";
      for (int i = 0; i < getMesh().GetNV(); i++)
      {
         os << i + 1
            << " " << getMesh().GetVertex(i)[0]
            << " " << getMesh().GetVertex(i)[1] << " " << getMesh().GetVertex(i)[2]
            << "\n"; 
      }
      os << "$EndNodes\n";

      // 3. Writing elements
      //====================
      int num_ele_bnd = getMesh().GetNBE(), num_ele = getMesh().GetNE();
      os << "$Elements\n";
      os << (num_ele_bnd + num_ele) << "\n";

      // 3.1. Writing boundary elements
      //===============================
      for (int i = 0; i < num_ele_bnd; i++)
      {
         int ele_type = 0;
         switch (getMesh().GetBdrElement(i)->GetType())
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
            << " " << getMesh().GetBdrElement(i)->GetAttribute();
         const int* vertices_list = getMesh().GetBdrElement(i)->GetVertices();
         for (int j = 0; j < (getMesh().GetBdrElement(i)->GetNVertices()); j++)
            os << " " << (vertices_list[j] + 1);
         os << "\n"; 
      }

      // 3.2. Writing bulk elements
      //===========================
      for (int i = 0; i < num_ele; i++)
      {
         int ele_type = 0;
         switch (getMesh().GetElement(i)->GetType())
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
            << " " << 1 << " " << getMesh().GetElement(i)->GetAttribute();

         const int* vertices_list = getMesh().GetElement(i)->GetVertices();
         for (int j = 0; j < getMesh().GetElement(i)->GetNVertices(); j++)
            os << " " << (vertices_list[j] + 1);
         os << "\n"; 
      }
      os << "$EndElements\n";

      return {true, {}};
   }
}
