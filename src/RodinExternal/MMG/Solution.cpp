#include "Utility.h"
#include "Solution.h"

namespace Rodin::External::MMG
{
  SolutionBase::SolutionBase(MMG5_pSol sol)
    : m_sol(sol)
  {
    if (!m_sol)
    {
      MMG5_SAFE_CALLOC(m_sol, 1, MMG5_Sol,
            Alert::Exception("Failed to allocate memory for MMG5_pSol").raise());
    }
  }

  SolutionBase::SolutionBase(const SolutionBase& other)
     : SolutionBase()
  {
    assert(other.m_sol);
    MMG5_Sol_Copy(other.getHandle(), getHandle());
  }

  SolutionBase::SolutionBase(SolutionBase&& other)
    : m_sol(other.m_sol)
  {
     assert(other.m_sol);
     other.m_sol = nullptr;
  }

  SolutionBase::~SolutionBase()
  {
    if (m_sol)
      MMG5_Sol_Free(m_sol);
  }

  SolutionBase& SolutionBase::operator=(SolutionBase&& other)
  {
     m_sol = other.m_sol;
     other.m_sol = nullptr;
     return *this;
  }
}
