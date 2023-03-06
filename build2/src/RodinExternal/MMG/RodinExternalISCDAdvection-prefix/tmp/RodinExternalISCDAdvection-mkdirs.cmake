# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/carlos/Projects/rodin/third-party/ISCD/Advection"
  "/Users/carlos/Projects/rodin/build2/third-party/ISCD/Advection"
  "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix"
  "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/tmp"
  "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/src/RodinExternalISCDAdvection-stamp"
  "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/src"
  "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/src/RodinExternalISCDAdvection-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/src/RodinExternalISCDAdvection-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/carlos/Projects/rodin/build2/src/RodinExternal/MMG/RodinExternalISCDAdvection-prefix/src/RodinExternalISCDAdvection-stamp${cfgdir}") # cfgdir has leading slash
endif()
