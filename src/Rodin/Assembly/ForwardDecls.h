#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

namespace Rodin::Assembly
{
  template <class LinearAlgebraType, class Operand>
  class AssemblyBase;

  template <class LinearAlgebraType, class Operand>
  class Sequential;

  template <class LinearAlgebraType, class Operand>
  class Multithreaded;

  template <class Operand>
  class OpenMP;

  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput;

  template <class FES>
  class LinearFormAssemblyInput;

  template <class ... Ts>
  class BilinearFormTupleAssemblyInput;

  template <class Scalar, class FES, class ValueDerived>
  class DirichletBCAssemblyInput;

  template <class Mesh>
  class SequentialIteration;

  template <class Mesh>
  class MultithreadedIteration;

  template <class ... Ts>
  class Default;
}

#endif
