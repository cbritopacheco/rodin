#ifndef RODIN_UTILITY_BOTTOMTEMPLATE_H
#define RODIN_UTILITY_BOTTOMTEMPLATE_H

namespace Rodin::Utility
{
  template <class T>
  struct BottomTemplate
  {
    using Type = T;
  };

  template <template <class> class T, class S>
  struct BottomTemplate<T<S>>
  {
    using Type = typename BottomTemplate<S>::Type;
  };
}

#endif
