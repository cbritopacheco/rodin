#ifndef RODIN_ASSEMBLY_FORWARDDECLS_H
#define RODIN_ASSEMBLY_FORWARDDECLS_H

namespace Rodin::Assembly
{
  template <class LinearAlgebraType, class Operand>
  class AssemblyBase;

  template <class LinearAlgebraType, class Operand>
  class Sequential;

  template <class Mesh>
  class SequentialIteration;

  template <class TrialFES, class TestFES>
  class BilinearFormAssemblyInput;

  template <class FES>
  class LinearFormAssemblyInput;

  template <class ... Ts>
  class BilinearFormTupleAssemblyInput;

  template <class Scalar, class Solution, class FES, class ValueDerived>
  class DirichletBCAssemblyInput;

  template <class ... Ts>
  class Default;

  template <class ... Ts>
  class Generic;
}

#ifdef RODIN_USE_OPENMP
namespace Rodin::Assembly
{
  template <class LinearAlgebraType, class Operand>
  class OpenMP;

  template <class Mesh>
  class OpenMPIteration;
}
#endif

#endif
