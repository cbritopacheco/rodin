// Object.cpp
#include "Object.h"
#include <cassert>
#include <petsc.h>

namespace Rodin::PETSc
{
  std::atomic<Object::Registry::Node*> Object::Registry::s_head{nullptr};

  Object::Object()
  {
    static auto registrar = []() -> bool {
      PetscBool init = PETSC_FALSE;
      PetscErrorCode ierr = PetscInitialized(&init);
      assert(ierr == PETSC_SUCCESS && init);

      ierr = PetscRegisterFinalize(&Registry::cleanup);
      assert(ierr == PETSC_SUCCESS);

      // ensure main thread’s node is pushed
      (void) Registry::local();
      return true;
    }();
    (void) registrar;

    // track *this* wrapper instance
    Registry::track(*this);
  }

  Object::~Object()
  {}

  struct Object::Registry::Node
  {
    Node* next{nullptr};
    std::vector<std::reference_wrapper<Object>> objs;

    void cleanup()
    {
      for (auto& ref : objs)
      {
        ::PetscObject& h = ref.get().getHandle();
        if (h)
          PetscObjectDestroy(&h);
      }
      objs.clear();
    }
  };

  Object::Registry::Node* Object::Registry::local()
  {
    thread_local Node node;
    static thread_local bool pushed =
      []()
      {
        Node* old = s_head.load(std::memory_order_relaxed);
        do
        {
          node.next = old;
        } while (!s_head.compare_exchange_weak(
                     old, &node,
                     std::memory_order_release,
                     std::memory_order_relaxed));
        return true;
      }();
    (void) pushed;
    return &node;
  }

  void Object::Registry::track(Object& o)
  {
    local()->objs.emplace_back(o);
  }

  PetscErrorCode Object::Registry::cleanup()
  {
    PetscFunctionBegin;
    Node* head = s_head.exchange(nullptr, std::memory_order_acquire);
    for (Node* cur = head; cur; cur = cur->next)
      cur->cleanup();
    PetscFunctionReturn(PETSC_SUCCESS);
  }
} // namespace Rodin::PETSc
