/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_H
#define RODIN_UTILITY_H

/**
 * @file
 * @brief Top level include for the Rodin::Utility module.
 *
 * The Utility module provides a comprehensive collection of template
 * metaprogramming utilities for type manipulation, compile-time computations,
 * and helper constructs used throughout the Rodin finite element library.
 *
 * ## Key Components
 *
 * - **Type Manipulation**: Extract, Wrap, Zip, Product for tuple type operations
 * - **Type Traits**: IsSpecialization, IsCompleteType, IsOneOf, HasTypeMember, HasValueMember
 * - **Parameter Pack Utilities**: ParameterPack, ForConstexpr, ForIndex
 * - **Helper Constructs**: Overloaded, Make, OptionalReference, Cast
 * - **Compile-Time Values**: False, DependentValue, IntegerSequence
 *
 * @see Rodin::Utility
 */

#include "Utility/Overloaded.h"
#include "Utility/IsSpecialization.h"
#include "Utility/False.h"

#endif
