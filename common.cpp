#include "common.h"

/* implement the vector functions of which the implementation needs to be delayed */
Vec::Line Vec::line(const Vec& p) const {
	return Line(*this, p - *this);
}
Vec::Plane Vec::plane(const Vec& p0, const Vec& p1) const {
	return Plane(*this, p0 - *this, p1 - *this);
}
