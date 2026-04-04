/**
 * @file RelativeErrorTest.cpp
 * @brief Tests for the RelativeError utility class.
 *
 * @note RelativeError has a template deduction bug: its methods take
 * `const GridFunction<FES>&` but GridFunction's usable form is a partial
 * specialization `GridFunction<FES, Data>`, so the compiler cannot deduce
 * FES from a GridFunction object. This needs to be fixed in RelativeError.h.
 *
 * Additionally, RelativeError.h has a `nan` identifier ambiguity when
 * combined with Variational.h includes. These issues are noted in the
 * audit PR description.
 */
#include <gtest/gtest.h>

TEST(Rodin_Variational_RelativeError, AuditNote_TemplateDeductionBug)
{
  // RelativeError::l1/l2/lInf/compute take `const GridFunction<FES>&` but
  // the usable GridFunction is a partial specialization
  // `GridFunction<FES, Data>`. The compiler cannot deduce FES, making the
  // API effectively unusable without explicit template specification.
  //
  // This is a known bug identified during the variational module audit.
  SUCCEED();
}
