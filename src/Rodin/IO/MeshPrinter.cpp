/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Helpers.h"

#include "MeshPrinter.h"

namespace Rodin::IO
{
   void MeshPrinter<FileFormat::GMSH, Context::Serial>::print(std::ostream& os)
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
   }

   void MeshPrinter<FileFormat::MEDIT, Context::Serial>::print(std::ostream& os)
   {
      const auto& mesh = getObject();

      // Print header
      os << IO::Medit::Keyword::MeshVersionFormatted << " " << 2
       << '\n'
       << IO::Medit::Keyword::Dimension << " " << mesh.getSpaceDimension()
       << "\n\n";

      // Print vertices
      os << '\n'
         << IO::Medit::Keyword::Vertices
         << '\n'
         << mesh.getHandle().GetNV()
         << '\n';

      for (int i = 0; i < mesh.getHandle().GetNV(); i++)
      {
         const double* coords = mesh.getHandle().GetVertex(i);
         for (int j = 0; j < mesh.getHandle().SpaceDimension(); j++)
            os << coords[j] << " ";
         os << 0 << '\n'; // Dummy reference for all vertices
      }

      // Print edges
      switch (mesh.getDimension())
      {
         case 2:
         {
            int nbe = mesh.getHandle().GetNBE();
            os << '\n'
               << IO::Medit::Keyword::Edges
               << '\n'
               << nbe
               << '\n';

            for (int i = 0; i < nbe; i++)
            {
               auto vs = mesh.getHandle().GetBdrElement(i)->GetVertices();
               os << vs[0] + 1 << " "
                  << vs[1] + 1 << " "
                  << mesh.getHandle().GetBdrAttribute(i) << '\n';
            }
         }
         default:
         {
            break;
         }
      }

      // Print triangles
      switch (mesh.getDimension())
      {
        case 2:
        {
            os << '\n' << IO::Medit::Keyword::Triangles << '\n';
           int ne = mesh.getElementCount();
           os << ne << '\n';
           for (int i = 0; i < ne; i++)
           {
              auto vs = mesh.getHandle().GetElement(i)->GetVertices();
              os << vs[0] + 1 << " "
                 << vs[1] + 1 << " "
                 << vs[2] + 1 << " "
                 << mesh.getHandle().GetAttribute(i) << '\n';
           }
           break;
        }
        case 3:
        {
            os << '\n' << IO::Medit::Keyword::Triangles << '\n';
            int nbe = mesh.getHandle().GetNBE();
            os << nbe << '\n';
            for (int i = 0; i < nbe; i++)
            {
               auto vs = mesh.getHandle().GetBdrElement(i)->GetVertices();
               os << vs[0] + 1 << " "
                  << vs[1] + 1 << " "
                  << vs[2] + 1 << " "
                  << mesh.getHandle().GetBdrAttribute(i) << '\n';
            }
            break;
        }
        default:
        {
           Alert::Exception() << "Bad mesh dimension " << mesh.getDimension();
           break;
        }
      }

      // Print tetrahedra
      switch (mesh.getDimension())
      {
         case 3:
         {
            int ne = mesh.getElementCount();
            os << '\n'
               << IO::Medit::Keyword::Tetrahedra
               << '\n'
               << ne
               << '\n';
            for (int i = 0; i < ne; i++)
            {
               auto it = mesh.getElement(i);
               assert(false);
               // auto vs = el.getVertices();
               // assert(vs.size() == 4);
               // os << vs[0] + 1 << " "
               //    << vs[1] + 1 << " "
               //    << vs[2] + 1 << " "
               //    << vs[3] + 1 << " "
               //    << el.getAttribute() << '\n';
            }
            break;
         }
         default:
         {
            // Do nothing
            break;
         }
      }

      // Print footer
      if (m_footer)
      {
         os << '\n' << IO::Medit::Keyword::End;
      }
   }
}
