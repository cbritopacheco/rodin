#ifndef RODIN_PETSC_OBJECT_H
#define RODIN_PETSC_OBJECT_H

#include <petsc.h>
#include <vector>
#include <atomic>
#include <functional>  // for reference_wrapper

namespace Rodin::PETSc
{
  /**
   * @brief Base‐class for any wrapper around a PETSc object.
   *
   * Each derived instance is automatically tracked in a per‐thread registry,
   * and any remaining PETSc handles are destroyed at PetscFinalize().
   */
  class Object
  {
    public:
      /**
       * @brief Registry for tracking all live wrappers.
       *
       * Implements a lock‐free, singly‐linked list of per‐thread Nodes.
       */
      struct Registry
      {
        /**
         * @brief Node in the lock‐free list of per‐thread registries.
         */
        struct Node;

        static std::atomic<Node*> s_head; ///< Head of the lock‐free list of Nodes.

        /**
         * @brief Return this thread’s registry node.
         *
         * On first invocation in a thread, the node is CAS‐pushed onto s_head.
         *
         * @return Pointer to the thread‐local Node.
         */
        static Node* local();

        /**
         * @brief Register one wrapper instance for later cleanup.
         *
         * @param o  Reference to the Object to be tracked.
         */
        static void track(Object& o);

        /**
         * @brief Callback invoked by PetscFinalize().
         *
         * Iterates over all registry nodes and destroys any leftover PETSc
         * objects via PetscObjectDestroy().
         *
         * @return PETSc error code (PETSC_SUCCESS on success).
         */
        static PetscErrorCode cleanup();
      };

      /**
       * @brief Construct and register this wrapper.
       *
       * Installs the PetscFinalize hook on first construction, and tracks
       * this instance in the per‐thread registry.
       */
      Object();

      /**
       * @brief Virtual destructor.
       *
       * Declared pure to force derived classes to define their own destructor,
       * ensuring proper PETSc cleanup.
       */
      virtual ~Object() = 0;

      /**
       * @brief Access the wrapped PETSc handle.
       *
       * @return Reference to the internal PetscObject.
       */
      virtual ::PetscObject& getHandle() noexcept = 0;
  };
}

#endif
