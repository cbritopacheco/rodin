#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

namespace Rodin::Assembly
{
  template <class Operand>
  class AssemblyBase;

  template <class Operand>
  class Serial;

  template <class Operand>
  class Multithreaded;

  template <class Operand>
  class OpenMP;
}

#endif
