/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXTERNAL_MMG_SCALARSOLUTION_H
#define RODIN_EXTERNAL_MMG_SCALARSOLUTION_H

#include <cassert>
#include <iterator>

#include "Solution.h"

namespace Rodin::External::MMG
{
   /**
    * @brief Scalar solution supported on a mesh.
    */
   class ScalarSolution : public SolutionBase
   {
      public:
         class Iterator
         {
            public:
               using iterator_category = std::random_access_iterator_tag;
               using difference_type   = std::ptrdiff_t;
               using value_type        = double;
               using pointer           = double*;
               using reference         = double&;

               Iterator(pointer ptr);
               reference operator*() const;
               pointer operator->();
               Iterator& operator++();
               Iterator operator++(int);
               friend bool operator==(const Iterator& a, const Iterator& b);
               friend bool operator!=(const Iterator& a, const Iterator& b);

            private:
                 pointer m_ptr;
         };

         /**
          * @internal
          */
         class ConstIterator
         {
            public:
               using iterator_category = std::forward_iterator_tag;
               using difference_type   = std::ptrdiff_t;
               using value_type        = double;
               using pointer           = double*;
               using reference         = double&;
               using const_reference   = const double&;

               ConstIterator(pointer ptr);

               const_reference operator*() const;
               pointer operator->();
               ConstIterator& operator++();
               ConstIterator operator++(int);
               friend bool operator==(const ConstIterator& a, const
                     ConstIterator& b);
               friend bool operator!=(const ConstIterator& a, const
                     ConstIterator& b);

            private:
                 pointer m_ptr;
         };

         /**
          * @returns Iterator to the beginning.
          */
         Iterator begin()
         {
            return Iterator(&getHandle()->m[1]);
         }

         /**
          * @returns Iterator to the end.
          */
         Iterator end()
         {
            return Iterator(&getHandle()->m[getHandle()->np + 1]);
         }

         /**
          * @returns ConstIterator to the beginning.
          */
         ConstIterator begin() const
         {
            return ConstIterator(&getHandle()->m[1]);
         }

         /**
          * @returns ConstIterator to the end.
          */
         ConstIterator end() const
         {
            return ConstIterator(&getHandle()->m[getHandle()->np + 1]);
         }
   };
}

#endif
