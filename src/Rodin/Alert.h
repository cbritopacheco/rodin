/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_H
#define RODIN_ALERT_H

/**
 * @file Alert.h
 * @brief Main include file for the Alert module.
 *
 * This header provides access to the complete Alert system including
 * exceptions, warnings, info messages, and success notifications with
 * colored terminal output support.
 *
 * Include this header to use the Alert module:
 * @code
 * #include <Rodin/Alert.h>
 * @endcode
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
