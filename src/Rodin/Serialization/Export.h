/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Export.h
 * @brief Polymorphic class export registrations for serialization
 * (e.g. polytope transformations), required to serialize objects
 * through base-class pointers.
 */
#ifndef RODIN_SERIALIZATION_EXPORT_H
#define RODIN_SERIALIZATION_EXPORT_H

#include <boost/serialization/export.hpp>
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/P1/P1Element.h"

// ParametricTransformation<FE> requires a scalar real-valued finite
// element (its static_assert enforces FE::RangeType == Real): geometry
// charts are built from scalar bases, with coordinates handled
// per-component. RealP1Element is therefore the only exportable P1
// instantiation; the previous Complex/Vector registrations could never
// compile (and VectorP1Element is an alias template usable only with an
// explicit scalar argument).
BOOST_CLASS_EXPORT(Rodin::Geometry::ParametricTransformation<Rodin::Variational::RealP1Element>);

#endif
