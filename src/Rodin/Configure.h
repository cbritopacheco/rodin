/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CONFIGURE_H
#define RODIN_CONFIGURE_H

/**
 * @file
 * @brief Header containing preprocessor directives for Rodin.
 */

/**
 * @defgroup RodinDirectives Preprocessor directives
 * @brief Preprocessor directives defined by Rodin.
 */

#define RODIN_VERSION "0.0.1"

#define RODIN_VERSION_MAJOR 0

#define RODIN_VERSION_MINOR 0

#define RODIN_VERSION_PATCH 1

#define RODIN_VERSION_CODE 000001

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximal space dimension.
 */
#define RODIN_MAXIMAL_SPACE_DIMENSION 3

/**
 * @ingroup RodinDirectives
 * @brief Indicates the default polytope attribute.
 */
#define RODIN_DEFAULT_POLYTOPE_ATTRIBUTE 1

/**
 * @brief Represents the constant used for fuzzy comparison.
 */
#define RODIN_FUZZY_CONSTANT 1e-5

#endif