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
 * @brief Indicates if Rodin is built with Python support.
 *
 * # Utilization
 *
 * @code{cpp}
 *  #ifndef RODIN_WITH_PY
 *  // Code depending on Python support
 *  #endif
 * @endcode
 */
/* #undef RODIN_WITH_PY */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin::Plot is built.
 * @code{cpp}
 *  #ifndef RODIN_WITH_PLOT
 *  // Code depending on Rodin::Plot support
 *  #endif
 * @endcode
 */
/* #undef RODIN_WITH_PLOT */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with MPI support.
 *
 * # Utilization
 *
 * @code{cpp}
 *  #ifndef RODIN_USE_MPI
 *  // Code depending on MPI support
 *  #endif
 * @endcode
 */
/* #undef RODIN_USE_MPI */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with OpenMP support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifndef RODIN_USE_OPENMP
 * // Code depending on OpenMP support
 * #endif
 * @endcode
 */
#define RODIN_USE_OPENMP

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with SuiteSparse support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifndef RODIN_USE_UMFPACK
 * // Code depending on UMFPACK
 * #endif
 * @endcode
 */
/* #undef RODIN_USE_UMFPACK */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with SuiteSparse support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifndef RODIN_USE_SPQR
 * // Code depending on SPQR
 * #endif
 * @endcode
 */
/* #undef RODIN_USE_SPQR */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with SuiteSparse support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifndef RODIN_USE_CHOLMOD
 * // Code depending on CHOLMOD
 * #endif
 * @endcode
 */
/* #undef RODIN_USE_CHOLMOD */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin warnings are silenced.
 */
/* #undef RODIN_SILENCE_WARNINGS */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin exceptions are silenced.
 *
 * If defined, this directive will prevent Rodin from outputting the error
 * messages. It does not prevent the exception from being thrown.
 */
#define RODIN_SILENCE_EXCEPTIONS

/**
 * @ingroup RodinDirectives
 * @brief Indicates the Rodin resources directory.
 *
 * # Utilization
 *
 * @code{cpp}
 * std::cout << RODIN_RESOURCES_DIR << std::endl;
 * @endcode
 */
#define RODIN_RESOURCES_DIR "/home/runner/work/rodin/rodin/resources/"

/**
 * @ingroup RodinDirectives
 * @brief Indicates the Rodin third-party source directory.
 *
 * # Utilization
 *
 * @code{cpp}
 * std::cout << RODIN_THIRD_PARTY_DIR << std::endl;
 * @endcode
 */
#define RODIN_THIRD_PARTY_DIR "/home/runner/work/rodin/rodin/third-party/"

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with thread safety enabled.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifdef RODIN_THREAD_SAFE
 *   // Code that requires thread safety features.
 * #endif
 * @endcode
 */
#define RODIN_THREAD_SAFE

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin supports multithreading.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifdef RODIN_MULTITHREADED
 *   // Code that leverages multithreading capabilities.
 * #endif
 * @endcode
 */
#define RODIN_MULTITHREADED

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with Apple Accelerate support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifdef RODIN_USE_APPLE_ACCELERATE
 *   // Code that uses the Apple Accelerate framework for optimized computations.
 * #endif
 * @endcode
 */
/* #undef RODIN_USE_APPLE_ACCELERATE */

/**
 * @ingroup RodinDirectives
 * @brief Indicates if Rodin is built with Scotch support.
 *
 * # Utilization
 *
 * @code{cpp}
 * #ifdef RODIN_USE_SCOTCH
 *   // Code that uses the Scotch library for partioning.
 * #endif
 * @endcode
 */
/* #undef RODIN_USE_SCOTCH */

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
