/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_H
#define RODIN_ALERT_H

/**
 * @file
 * @brief Main header file for the Alert module.
 *
 * This file provides the main interface to the Alert module, which offers
 * comprehensive facilities for error handling, warnings, informational messages,
 * and success notifications with colored terminal output capabilities.
 *
 * The Alert module includes:
 * - Exception handling with formatted error messages
 * - Warning and informational message systems
 * - Success notifications
 * - Terminal text styling and coloring
 * - Class and namespace-specific exception handling
 *
 * @see @ref AlertModule
 */

#include "Alert/Exception.h"
#include "Alert/Success.h"
#include "Alert/Warning.h"
#include "Alert/Info.h"

#include "Alert/Text.h"
#include "Alert/Color.h"
#include "Alert/Reset.h"
#include "Alert/Stylize.h"
#include "Alert/Notation.h"
#include "Alert/Identifier.h"
#include "Alert/ClassException.h"
#include "Alert/NamespacedException.h"
#include "Alert/MemberFunctionException.h"

#include "Alert/MemberFunctionWarning.h"

#endif
