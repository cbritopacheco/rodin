#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

namespace Rodin::Assembly
{
  template <class LinearAlgebraType, class ... Args>
  class AssemblyBase;

  template <class LinearAlgebraType, class Operand>
  class Sequential;

  template <class LinearAlgebraType, class Operand>
  class Multithreaded;

  template <class Operand>
  class OpenMP;
}

#endif