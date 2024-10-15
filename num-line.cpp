#include "vec.h"

num::Line::Line() {}
num::Line::Line(const Vec& _d) : d(_d) {}
num::Line::Line(const Vec& _o, const Vec& _d) : o(_o), d(_d) {}
num::Linear num::Line::fLinComb(const Line& l, size_t index, bool& parallel, float precision) const {
	/*
	*	E: o + s * d
	*	F: l.o + t * l.d
	*
	*	Solve for s and insert yields: s = l.d.y * (l.o.x - o.x) - l.d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*								   t = d.y * (l.o.x - o.x) - d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*	(same for x-z/y-z)
	*/

	/* extract the two components to work with and compute the divisor */
	const size_t _0 = (index + 1) % 3;
	const size_t _1 = (index + 2) % 3;
	const float divisor = d.c[_0] * l.d.c[_1] - d.c[_1] * l.d.c[_0];

	/* check if the two lines are parallel */
	if (parallel = (std::abs(divisor) <= precision))
		return Linear();

	/* compute the linear combination */
	const float s = (l.d.c[_1] * (l.o.c[_0] - o.c[_0]) - l.d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
	const float t = (d.c[_1] * (l.o.c[_0] - o.c[_0]) - d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
	return Linear(s, t);
}
num::Line num::Line::AxisX(float l) {
	return Line(Vec(), Vec::AxisX(l));
}
num::Line num::Line::AxisY(float l) {
	return Line(Vec(), Vec::AxisY(l));
}
num::Line num::Line::AxisZ(float l) {
	return Line(Vec(), Vec::AxisZ(l));
}
num::Line num::Line::planeX(float xPlane) const {
	return Line(o.planeX(xPlane), d.planeX(xPlane));
}
num::Line num::Line::planeY(float yPlane) const {
	return Line(o.planeY(yPlane), d.planeY(yPlane));
}
num::Line num::Line::planeZ(float zPlane) const {
	return Line(o.planeZ(zPlane), d.planeZ(zPlane));
}
num::Vec num::Line::point(float t) const {
	return o + d * t;
}
num::Line num::Line::norm() const {
	/*
	*	o - a * d = NewOrigin
	*	-> (o - a * d) * d = 0 (the line and NewOrigin should be perpendicular)
	*	Solve for a
	*/
	const float a = o.dot(d) / d.dot(d);
	return Line(o - d * a, d.norm());
}
bool num::Line::touch(const Vec& p, float precision) const {
	/*
	*	E: o + s * d
	*	Set equal to [p] and solve for s and insert
	*/

	/* find the largest component which should not be zero if the line is well defined */
	size_t index = d.comp(true);
	const float s = (p.c[index] - o.c[index]) / d.c[index];

	/* compute the point on the line where the given point is expected to be */
	const Vec t = o + d * s;

	/* ensure that the points are equal */
	return p.equal(t, precision);
}
float num::Line::find(const Vec& p) const {
	/* extract the largest component of the direction and use it to compute the scaling factor */
	size_t index = d.comp(true);
	return (p.c[index] - o.c[index]) / d.c[index];
}
bool num::Line::equal(const Line& l, float precision) const {
	return l.touch(o, precision) && l.d.parallel(d, precision);
}
bool num::Line::identical(const Line& l, float precision) const {
	return l.o.equal(o, precision) && l.d.equal(d, precision);
}
float num::Line::closestf(const Vec& p) const {
	/*
	*	o + a * d = p + v
	*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
	*	Solve for a
	*/
	return (p - o).dot(d) / d.dot(d);
}
num::Vec num::Line::closest(const Vec& p) const {
	const float a = closestf(p);
	return o + d * a - p;
}
num::Linear num::Line::closestf(const Line& l) const {
	/*
	*	E: o + s * d
	*	F: l.o + t * l.d
	*
	*	To solve: o + s * d + r * v = l.o + t * l.d
	*	both d and l.d have to be perpendicular to v
	*		=> v = d x l.d
	*
	*	Solution: r = ((l.o - o) * (d x l.d)) / v * (d x l.d)
	*	Solution: s = ((o - l.o) * (v x l.d)) / v * (d x l.d)
	*	Solution: t = ((o - l.o) * (v x d)) / v * (d x l.d)
	*/
	const Vec v = d.cross(l.d);

	/* check if the lines run in parallel */
	if (v.zero())
		return Linear(0.0f, l.closestf(o));

	/* compute the two scalars */
	const float tmp = v.dot(v);
	const Vec df = l.o - o;
	const float s = -df.dot(v.cross(l.d)) / tmp;
	const float t = -df.dot(v.cross(d)) / tmp;

	/* return the two factors */
	return Linear(s, t);
}
num::Line num::Line::closest(const Line& l) const {
	/*
	*	E: o + s * d
	*	F: l.o + t * l.d
	*
	*	To solve: o + s * d + r * v = l.o + t * l.d
	*	both d and l.d have to be perpendicular to v
	*		=> v = d x l.d
	*
	*	Solution: r = ((l.o - o) * (d x l.d)) / v * (d x l.d)
	*	Solution: s = ((o - l.o) * (v x l.d)) / v * (d x l.d)
	*	Solution: t = ((o - l.o) * (v x d)) / v * (d x l.d)
	*/
	const Vec v = d.cross(l.d);

	/* check if the lines run in parallel */
	if (v.zero())
		return Line(o, l.closest(o));

	/* compute the two scalars */
	const float tmp = v.dot(v);
	const Vec df = l.o - o;
	const float r = df.dot(v) / tmp;
	const float s = -df.dot(v.cross(l.d)) / tmp;

	/* return the final line */
	return Line(o + d * s, v * r);
}
float num::Line::intersectPlaneXf(float xPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.x) <= precision) {
		if (invalid)
			*invalid = true;
		return 0.0f;
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: x = xPlane
	*	Line: o + a * d
	*	insert into equation for x and solve for a
	*/
	return (xPlane - o.x) / d.x;
}
num::Vec num::Line::intersectPlaneX(float xPlane, bool* invalid, float precision) const {
	const float a = intersectPlaneXf(xPlane, invalid, precision);
	return o + d * a;
}
float num::Line::intersectPlaneYf(float yPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.y) <= precision) {
		if (invalid)
			*invalid = true;
		return 0.0f;
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: y = yPlane
	*	Line: o + a * d
	*	insert into equation for y and solve for a
	*/
	const float a = (yPlane - o.y) / d.y;
	return a;
}
num::Vec num::Line::intersectPlaneY(float yPlane, bool* invalid, float precision) const {
	const float a = intersectPlaneYf(yPlane, invalid, precision);
	return o + d * a;
}
float num::Line::intersectPlaneZf(float zPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.z) <= precision) {
		if (invalid)
			*invalid = true;
		return 0.0f;
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: z = zPlane
	*	Line: o + a * d
	*	insert into equation for z and solve for a
	*/
	const float a = (zPlane - o.z) / d.z;
	return a;
}
num::Vec num::Line::intersectPlaneZ(float zPlane, bool* invalid, float precision) const {
	const float a = intersectPlaneZf(zPlane, invalid, precision);
	return o + d * a;
}
num::Linear num::Line::intersectXf(const Line& l, bool* invalid, float precision) const {
	bool parallel = false;

	/* compute the linear combination */
	const Linear lin = fLinComb(l, Vec::Component::ComponentX, parallel, precision);

	/* update the invalid flag and otherwise return the result */
	if (invalid)
		*invalid = parallel;
	return lin;
}
num::Vec num::Line::intersectX(const Line& l, bool* invalid, float precision) const {
	const float s = intersectXf(l, invalid, precision).s;
	return o + d * s;
}
num::Linear num::Line::intersectYf(const Line& l, bool* invalid, float precision) const {
	bool parallel = false;

	/* compute the linear combination */
	const Linear lin = fLinComb(l, Vec::Component::ComponentY, parallel, precision);

	/* update the invalid flag and otherwise return the result */
	if (invalid)
		*invalid = parallel;
	return lin;
}
num::Vec num::Line::intersectY(const Line& l, bool* invalid, float precision) const {
	const float s = intersectYf(l, invalid, precision).s;
	return o + d * s;
}
num::Linear num::Line::intersectZf(const Line& l, bool* invalid, float precision) const {
	bool parallel = false;

	/* compute the linear combination */
	const Linear lin = fLinComb(l, Vec::Component::ComponentZ, parallel, precision);

	/* update the invalid flag and otherwise return the result */
	if (invalid)
		*invalid = parallel;
	return lin;
}
num::Vec num::Line::intersectZ(const Line& l, bool* invalid, float precision) const {
	const float s = intersectZf(l, invalid, precision).s;
	return o + d * s;
}
num::Linear num::Line::intersectf(const Line& l, bool* invalid, float precision) const {
	/*
	*	Select the axis to compute the combination for by computing the cross product between
	*	the two and then selecting the smallest component which ensures that the other two components
	*	are larger, as long as the lines are well defined
	*/
	const size_t index = d.cross(l.d).comp(false);
	bool parallel = false;

	/* compute the linear combination */
	const Linear lin = fLinComb(l, index, parallel, precision);

	/* check if the lines intersect */
	bool on = false;
	if (!parallel)
		on = num::Cmp(o.c[index] + d.c[index] * lin.s, l.o.c[index] + l.d.c[index] * lin.t, precision);

	/* update the invalid flag and return the result */
	if (invalid)
		*invalid = !on;
	return on ? lin : Linear();
}
num::Vec num::Line::intersect(const Line& l, bool* invalid, float precision) const {
	const float f = intersectf(l, invalid, precision).s;
	return o + d * f;
}

std::ostream& num::operator<<(std::ostream& out, const Line& l) {
	return (out << l.o << " -> " << l.d);
}
std::wostream& num::operator<<(std::wostream& out, const Line& l) {
	return (out << l.o << L" -> " << l.d);
}

std::istream& num::operator>>(std::istream& in, Line& l) {
	char pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != '-' || pad1 != '>')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& num::operator>>(std::wistream& in, Line& l) {
	wchar_t pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != L'-' || pad1 != L'>')
		in.setstate(std::ios::failbit);
	return in;
}
