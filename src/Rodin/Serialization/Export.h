#ifndef RODIN_SERIALIZATION_EXPORT_H
#define RODIN_SERIALIZATION_EXPORT_H

#include <boost/serialization/export.hpp>
#include "Rodin/Geometry/IsoparametricTransformation.h"
#include "Rodin/Variational/P1/P1Element.h"

BOOST_CLASS_EXPORT(Rodin::Geometry::IsoparametricTransformation<Rodin::Variational::RealP1Element>);
BOOST_CLASS_EXPORT(Rodin::Geometry::IsoparametricTransformation<Rodin::Variational::ComplexP1Element>);
BOOST_CLASS_EXPORT(Rodin::Geometry::IsoparametricTransformation<Rodin::Variational::VectorP1Element>);

#endif
