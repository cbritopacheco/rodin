/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_HPP
#define RODIN_VARIATIONAL_GRIDFUNCTIONLOADER_HPP

#include <boost/algorithm/string.hpp>

#include "Rodin/Variational/GridFunction.h"

#include "Helpers.h"

#include "GridFunctionLoader.h"

namespace Rodin::IO
{
   template <class FEC>
   void GridFunctionLoader<FileFormat::MFEM, FEC, Traits::Serial>
   ::load(std::istream& is)
   {
      auto& gf = GridFunctionLoaderBase<FEC, Traits::Serial>::getObject();
      gf.getHandle() = mfem::GridFunction(
            &gf.getFiniteElementSpace().getMesh().getHandle(), is);
      gf.getHandle().SetSpace(&gf.getFiniteElementSpace().getHandle());
   }

   template <class FEC>
   void GridFunctionLoader<FileFormat::MEDIT, FEC, Traits::Serial>
   ::load(std::istream& is)
   {
      std::string line;
      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (!line.empty())
            break;
      }

      std::istringstream iss(line);
      std::string kw;
      iss >> kw;

      int version;
      if (Medit::KeywordMap.left.count(kw))
      {
         switch (Medit::KeywordMap.left.at(kw))
         {
            case Medit::Keyword::MeshVersionFormatted:
            {
               if (!(iss >> version)) // Version is not on this line
               {
                  // Try next line
                  std::getline(is, line);
                  if (!(std::istringstream(line) >> version))
                  {
                     Alert::Exception()
                        << "Bad solution format. Encountered \""
                        << line << "\" after keyword \""
                        << Medit::Keyword::MeshVersionFormatted
                        << "\"."
                        << Alert::Raise;
                  }
               }
               break;
            }
            default:
            {
               Alert::Exception()
                  << "Bad solution format. Encountered keyword \""
                  << kw << "\" before keyword \""
                  << Medit::Keyword::MeshVersionFormatted
                  << "\"."
                  << Alert::Raise;
               break;
            }
         }
      }
      else
      {
         Alert::Exception()
            << "Bad solution format. Encountered \""
            << line << "\" before keyword \""
            << Medit::Keyword::MeshVersionFormatted
            << "\"."
            << Alert::Raise;
      }

      assert(version >= 2 && version <= 3);

      int spaceDim;
      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (line.empty())
            continue;

         std::istringstream iss(line);
         std::string kw;
         iss >> kw;
         if (kw == Medit::KeywordMap.right.at(Medit::Keyword::Dimension))
         {
            if (!(iss >> spaceDim)) // Dimension not on this line
            {
               // Try next line
               std::getline(is, line);
               if (!(std::istringstream(line) >> spaceDim))
               {
                  Alert::Exception()
                     << "Bad solution format. Encountered \""
                     << line << "\" after keyword \""
                     << Medit::Keyword::Dimension
                     << "\"."
                     << Alert::Raise;
               }
            }
            if (spaceDim < 2 || spaceDim > 3)
               Alert::Exception() << "Invalid dimension " << spaceDim << Alert::Raise;
            break;
         }
         else
         {
            Alert::Exception()
               << "Bad solution format. Encountered keyword \"" << kw
               << "\" before keyword \""
               << Medit::Keyword::Dimension
               << "\"."
               << Alert::Raise;
            return;
         }
      }

      while (std::getline(is, line))
      {
         boost::algorithm::trim(line);
         if (line.empty())
            continue;
         std::string kw;
         std::istringstream(line) >> kw;
         if (!std::isalpha(kw[0]))
            continue;
         if (Medit::KeywordMap.left.count(kw))
         {
            switch (Medit::KeywordMap.left.at(kw))
            {
               case Medit::Keyword::SolAtVertices:
               {
                  int elementCount;
                  std::getline(is, line);
                  std::istringstream(line) >> elementCount;

                  int nsols, solType;
                  std::getline(is, line);
                  std::istringstream(line) >> nsols >> solType;
                  if (nsols > 1)
                  {
                     Alert::Exception()
                        << "Solution number greater than 1 is not supported."
                        << Alert::Raise;
                  }

                  int vdim;
                  switch (solType)
                  {
                     case 1:
                     {
                        vdim = 1;
                        break;
                     }
                     case 2:
                     {
                        vdim = spaceDim;
                        break;
                     }
                     case 3:
                     {
                        Alert::Exception()
                           << "Solution type \"" << solType
                           << "\" (tensor) is not supported."
                           << Alert::Raise;
                        break;
                     }
                     default:
                     {
                        Alert::Exception()
                           << "Unknown solution type \"" << solType << "\"."
                           << Alert::Raise;
                     }
                  }

                  if (this->getObject().getFiniteElementSpace().getVectorDimension() != vdim)
                  {
                     Alert::Exception()
                        << "Mismatching dimensions (expected) "
                        << this->getObject().getFiniteElementSpace().getVectorDimension()
                        << " != "
                        << vdim
                        << " (loaded)."
                        << Alert::Raise;
                  }

                  // Read the whole data
                  switch (this->getObject().getFiniteElementSpace().getHandle().GetOrdering())
                  {
                     case mfem::Ordering::byNODES:
                     {
                        for (int i = 0; i < elementCount; i++)
                        {
                           while (std::getline(is, line))
                           {
                              boost::algorithm::trim(line);
                              if (!line.empty())
                                 break;
                           }
                           std::istringstream iss(line);
                           for (int j = 0; j < vdim; j++)
                           {
                              if(!(iss >> this->getObject().getHandle()[i + j * elementCount]))
                              {
                                 Alert::Exception()
                                    << "Bad solution format. Could not parse value in line \""
                                    << line
                                    << "\"."
                                    << Alert::Raise;
                              }
                           }
                        }
                        break;
                     }
                     case mfem::Ordering::byVDIM:
                     {
                        for (int i = 0; i < elementCount; i++)
                        {
                           while (std::getline(is, line))
                           {
                              boost::algorithm::trim(line);
                              if (!line.empty())
                                 break;
                           }
                           std::istringstream iss(line);
                           for (int j = 0; j < vdim; j++)
                           {
                              if (!(iss >> this->getObject().getHandle()[j + i * vdim]))
                              {
                                 Alert::Exception()
                                    << "Bad solution format. Could not parse value in line \""
                                    << line
                                    << "\"."
                                    << Alert::Raise;
                              }
                           }
                        }
                        break;
                     }
                  }
                  continue;
               }
               case Medit::Keyword::End:
               {
                  break;
               }
               default:
               {
                  Alert::Exception()
                     << "Encountered \"" << kw << "\". Only "
                     << Medit::Keyword::SolAtVertices
                     << " is supported."
                     << Alert::Raise;
                  continue;
               }
            }
            break;
         }
         else
         {
            Alert::Warning()
               << "Ignoring unrecognized keyword \"" << kw << "\""
               << Alert::Raise;
         }
      }
   }
}

#endif
