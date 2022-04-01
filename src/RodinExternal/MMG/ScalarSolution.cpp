#include "ScalarSolution.h"

namespace Rodin::External::MMG
{
   // ---- Iterator implementation -------------------------------------------
   ScalarSolution::Iterator::Iterator(pointer ptr)
      : m_ptr(ptr)
   {}

   ScalarSolution::Iterator::reference
   ScalarSolution::Iterator::operator*()
   const
   {
      return *m_ptr;
   }

   ScalarSolution::Iterator::pointer
   ScalarSolution::Iterator::operator->()
   {
      return m_ptr;
   }

   ScalarSolution::Iterator&
   ScalarSolution::Iterator::operator++()
   {
      m_ptr++;
      return *this;
   }

   ScalarSolution::Iterator
   ScalarSolution::Iterator::operator++(int)
   {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
   }

   bool operator==(
         const ScalarSolution::Iterator& a,
         const ScalarSolution::Iterator& b)
   {
      return a.m_ptr == b.m_ptr;
   };

   bool operator!=(
         const ScalarSolution::Iterator& a,
         const ScalarSolution::Iterator& b)
   {
      return a.m_ptr != b.m_ptr;
   };

   // ---- ConstIterator implementation --------------------------------------
   ScalarSolution::ConstIterator::ConstIterator(pointer ptr)
      : m_ptr(ptr)
   {}

   ScalarSolution::ConstIterator::const_reference
   ScalarSolution::ConstIterator::operator*() const
   {
      return *m_ptr;
   }

   ScalarSolution::ConstIterator::pointer
   ScalarSolution::ConstIterator::operator->()
   {
      return m_ptr;
   }

   ScalarSolution::ConstIterator&
   ScalarSolution::ConstIterator::operator++()
   {
      m_ptr++;
      return *this;
   }

   ScalarSolution::ConstIterator
   ScalarSolution::ConstIterator::operator++(int)
   {
      ConstIterator tmp = *this;
      ++(*this);
      return tmp;
   }

   bool
   operator==(
         const ScalarSolution::ConstIterator& a,
         const ScalarSolution::ConstIterator& b)
   {
      return a.m_ptr == b.m_ptr;
   };

   bool
   operator!=(
         const ScalarSolution::ConstIterator& a,
         const ScalarSolution::ConstIterator& b)
   {
      return a.m_ptr != b.m_ptr;
   };
}
