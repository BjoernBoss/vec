#include "vec.h"

/* implement the vector object */
num::Vec::Vec() : x(0.0f), y(0.0f), z(0.0f) {}
num::Vec::Vec(float f) : x(f), y(f), z(f) {}
num::Vec::Vec(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
num::Vec num::Vec::operator+(const Vec& v) const {
	return Vec(x + v.x, y + v.y, z + v.z);
}
num::Vec num::Vec::operator-(const Vec& v) const {
	return Vec(x - v.x, y - v.y, z - v.z);
}
num::Vec num::Vec::operator-() const {
	return Vec(-x, -y, -z);
}
num::Vec num::Vec::operator*(float s) const {
	return Vec(x * s, y * s, z * s);
}
num::Vec num::Vec::operator/(float s) const {
	return Vec(x / s, y / s, z / s);
}
num::Vec& num::Vec::operator+=(const Vec& v) {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}
num::Vec& num::Vec::operator-=(const Vec& v) {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}
num::Vec& num::Vec::operator*=(float s) {
	x *= s;
	y *= s;
	z *= s;
	return *this;
}
num::Vec& num::Vec::operator/=(float s) {
	x /= s;
	y /= s;
	z /= s;
	return *this;
}
bool num::Vec::operator==(const Vec& v) const {
	return equal(v);
}
bool num::Vec::operator!=(const Vec& v) const {
	return !equal(v);
}
float num::Vec::dot(const Vec& v) const {
	return v.x * x + v.y * y + v.z * z;
}
float num::Vec::angle(const Vec& v) const {
	float dotProd = dot(v);
	float lenProd = std::sqrt(lenSquared() * v.lenSquared());
	float frac = dotProd / lenProd;

	/* check if the angle reaches the degrees where the uncertainty of floats
		might lead to values outside of the scope of arccos */
	if (frac >= 1.0f)
		return 0.0f;
	else if (frac <= -1.0f)
		return 180.0f;
	return num::ToDegree(std::acos(frac));
}
float num::Vec::len() const {
	return std::sqrt(x * x + y * y + z * z);
}
float num::Vec::lenSquared() const {
	return x * x + y * y + z * z;
}
num::Vec num::Vec::cross(const Vec& v) const {
	return Vec(
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x);
}
float num::Vec::crossX(const Vec& v) const {
	return y * v.z - z * v.y;
}
float num::Vec::crossY(const Vec& v) const {
	return z * v.x - x * v.z;
}
float num::Vec::crossZ(const Vec& v) const {
	return x * v.y - y * v.x;
}
num::Vec num::Vec::norm() const {
	float l = len();
	return Vec(x / l, y / l, z / l);
}
num::Vec num::Vec::planeX(float xPlane) const {
	return Vec(xPlane, y, z);
}
num::Vec num::Vec::planeY(float yPlane) const {
	return Vec(x, yPlane, z);
}
num::Vec num::Vec::planeZ(float zPlane) const {
	return Vec(x, y, zPlane);
}
size_t num::Vec::comp(bool largest) const {
	size_t index = 0;

	/* iterate through the components and check if one is larger */
	for (size_t i = 1; i < 3; i++) {
		if (largest ? std::abs(c[index]) >= std::abs(c[i]) : std::abs(c[index]) <= std::abs(c[i]))
			continue;
		index = i;
	}
	return index;
}
num::Vec num::Vec::rotateX(float a) const {
	a = num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x,
		y * ca - z * sa,
		y * sa + z * ca);
}
num::Vec num::Vec::rotateY(float a) const {
	a = num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x * ca + z * sa,
		y,
		z * ca - x * sa);
}
num::Vec num::Vec::rotateZ(float a) const {
	a = -num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x * ca - y * sa,
		x * sa + y * ca,
		z);
}
float num::Vec::angleX(const Vec& r) const {
	Vec ref = r.planeX();
	Vec flat = planeX();

	/* compute the angle between the two vectors and correct its sign */
	float angle = ref.angle(flat);
	return (ref.crossX(flat) > 0.0f) ? -angle : angle;
}
float num::Vec::angleY(const Vec& r) const {
	Vec ref = r.planeY();
	Vec flat = planeY();

	/* compute the angle between the two vectors and correct its sign */
	float angle = ref.angle(flat);
	return (ref.crossY(flat) > 0.0f) ? -angle : angle;
}
float num::Vec::angleZ(const Vec& r) const {
	Vec ref = r.planeZ();
	Vec flat = planeZ();

	/* compute the angle between the two vectors and correct its sign */
	float angle = ref.angle(flat);
	return (ref.crossZ(flat) > 0.0f) ? -angle : angle;
}
num::Line num::Vec::line(const Vec& p) const {
	return Line(*this, p - *this);
}
num::Plane num::Vec::plane(const Vec& p0, const Vec& p1) const {
	return Plane(*this, p0 - *this, p1 - *this);
}
num::Vec num::Vec::interpolate(const Vec& p, float t) const {
	return Vec(
		x + (p.x - x) * t,
		y + (p.y - y) * t,
		z + (p.z - z) * t
	);
}
num::Vec num::Vec::rescale(float l) const {
	float factor = std::sqrt((l * l) / lenSquared());
	return *this * factor;
}
float num::Vec::delta(const Vec& v) const {
	/* extract the largest component and use it to compute the scaling factor */
	size_t index = comp(true);
	return v.c[index] / c[index];
}
num::Vec num::Vec::scale(float f) const {
	return Vec(x * f, y * f, z * f);
}
bool num::Vec::parallel(const Vec& v, float precision) const {
	/* extract the largest components and check if the vectors are considered zero */
	const size_t largest[2] = { comp(true), v.comp(true) };
	if (std::abs(c[largest[0]]) <= precision)
		return (std::abs(v.c[largest[1]]) <= precision);
	else if (std::abs(v.c[largest[1]]) <= precision)
		return false;

	/* compute the factor with which to scale the other vector */
	const float f = c[largest[0]] / v.c[largest[1]];

	/* check if the vectors are equal when scaled */
	return equal(v * f, precision);
}
bool num::Vec::sign(const Vec& v, float precision) const {
	/* extract the largest components and check if the vectors are considered zero */
	const size_t largest[2] = { comp(true), v.comp(true) };
	if (std::abs(c[largest[0]]) <= precision)
		return (std::abs(v.c[largest[1]]) <= precision);
	else if (std::abs(v.c[largest[1]]) <= precision)
		return false;

	/* compute the factor with which to scale the other vector and ensure that the factor is positive */
	const float f = c[largest[0]] / v.c[largest[1]];
	if (f < 0.0f)
		return false;

	/* check if the vectors are equal when scaled */
	return equal(v * f, precision);
}
bool num::Vec::match(const Vec& v, float precision) const {
	return num::Cmp(dot(v), lenSquared(), precision);
}
bool num::Vec::equal(const Vec& v, float precision) const {
	/* dont subtract and then compare with zero as small errors will have a much larger effect on
	*	the result due to the canceling effects of subtraction on the information */
	return num::Cmp(x, v.x, precision) && num::Cmp(y, v.y, precision) && num::Cmp(z, v.z, precision);
}
bool num::Vec::zeroX(float precision) const {
	return num::Zero(x, precision);
}
bool num::Vec::zeroY(float precision) const {
	return num::Zero(y, precision);
}
bool num::Vec::zeroZ(float precision) const {
	return num::Zero(z, precision);
}
bool num::Vec::zero(float precision) const {
	return num::Zero(lenSquared(), precision);
}
num::Vec num::Vec::perpendicular(const Vec& v) const {
	/*
	*	(v - [this] * a) * [this] = 0
	*	Solve for a and insert
	*/
	const float a = dot(v) / dot(*this);
	return v - *this * a;
}
num::Vec num::Vec::project(const Vec& v) const {
	float factor = dot(v) / lenSquared();
	return *this * factor;
}
float num::Vec::projectf(const Vec& v) const {
	return dot(v) / lenSquared();
}

/* implement the vector multiplication from the right */
num::Vec num::operator*(float s, const Vec& v) {
	return v * s;
}

/* implement the line object */
num::Line::Linear::Linear() : s(0.0f), t(0.0f) {}
num::Line::Linear::Linear(float _s, float _t) : s(_s), t(_t) {}
num::Line::Line() {}
num::Line::Line(const Vec& _d) : d(_d) {}
num::Line::Line(const Vec& _o, const Vec& _d) : o(_o), d(_d) {}
num::Line::Linear num::Line::fLinComb(const Line& l, size_t index, bool& parallel, float precision) const {
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
num::Line::Linear num::Line::closestf(const Line& l) const {
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
num::Line::Linear num::Line::intersectXf(const Line& l, bool* invalid, float precision) const {
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
num::Line::Linear num::Line::intersectYf(const Line& l, bool* invalid, float precision) const {
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
num::Line::Linear num::Line::intersectZf(const Line& l, bool* invalid, float precision) const {
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
num::Line::Linear num::Line::intersectf(const Line& l, bool* invalid, float precision) const {
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

/* implement the plane object */
num::Plane::Linear::Linear() : s(0.0f), t(0.0f) {}
num::Plane::Linear::Linear(float _s, float _t) : s(_s), t(_t) {}
num::Plane::Plane() {}
num::Plane::Plane(const Vec& _a, const Vec& _b) : a(_a), b(_b) {}
num::Plane::Plane(const Vec& _o, const Vec& _a, const Vec& _b) : o(_o), a(_a), b(_b) {}
num::Plane::Linear num::Plane::fLinComb(const Vec& p, size_t index) const {
	/*
	*	compute the linar combination that results in the point [p] while ignoring the component passed in as index (here in X-Y plane)
	*	p = o + s * a + t * b
	*	s = ((p.x - o.x) * b.y - (p.y - o.y) * b.x) / (a.x * b.y - a.y * b.x)
	*	t = (a.x * (p.y - o.y) - a.y * (p.x - o.x)) / (a.x * b.y - a.y * b.x)
	*/
	const size_t _0 = (index + 1) % 3;
	const size_t _1 = (index + 2) % 3;

	const float divisor = a.c[_0] * b.c[_1] - a.c[_1] * b.c[_0];

	const float _v0 = p.c[_0] - o.c[_0];
	const float _v1 = p.c[_1] - o.c[_1];
	const float _s = (_v0 * b.c[_1] - _v1 * b.c[_0]) / divisor;
	const float _t = (a.c[_0] * _v1 - a.c[_1] * _v0) / divisor;
	return Linear(_s, _t);
}
num::Plane num::Plane::planeX(float xPlane) const {
	return Plane(o.planeX(xPlane), a.planeX(xPlane), b.planeX(xPlane));
}
num::Plane num::Plane::planeY(float yPlane) const {
	return Plane(o.planeY(yPlane), a.planeY(yPlane), b.planeY(yPlane));
}
num::Plane num::Plane::planeZ(float zPlane) const {
	return Plane(o.planeZ(zPlane), a.planeZ(zPlane), b.planeZ(zPlane));
}
num::Vec num::Plane::normal() const {
	return a.cross(b);
}
num::Vec num::Plane::center() const {
	return o + ((a + b) / 3);
}
num::Vec num::Plane::point(float s, float t) const {
	return o + a * s + b * t;
}
num::Vec num::Plane::point(const Linear& lin) const {
	return o + a * lin.s + b * lin.t;
}
num::Plane num::Plane::norm() const {
	/*
	*	compute this by computing the intersection point of a line with the
	*	planes normal vector as direction vector and the zero-vector as its origin
	*/
	const Vec crs = a.cross(b);

	/*
	*	E: o + s * x1 + t * x2
	*	G: [this] + f * x0
	*	Solve for f and insert
	*/
	const float f = o.dot(crs) / crs.dot(crs);
	const Vec _a = a.norm();
	return Plane(crs * f, _a, _a.perpendicular(b).norm());
}
num::Vec num::Plane::projectX(const Vec& p) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentX);
	return Vec(o.x + r.s * a.x + r.t * b.x, p.y, p.z);
}
num::Vec num::Plane::projectY(const Vec& p) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentY);
	return Vec(p.x, o.y + r.s * a.y + r.t * b.y, p.z);
}
num::Vec num::Plane::projectZ(const Vec& p) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentZ);
	return Vec(p.x, p.y, o.z + r.s * a.z + r.t * b.z);
}
num::Vec num::Plane::project(const Vec& v) const {
	/*
	*	compute the projection onto the normal of the plane and subtract it from p
	*	as this will result in only the part on the vector within the plane
	*/
	const Vec crs = a.cross(b);
	return v - crs.project(v);
}
bool num::Plane::inTriangleX(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentX);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool num::Plane::inTriangleY(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentY);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool num::Plane::inTriangleZ(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentZ);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool num::Plane::inTriangle(const Vec& p, bool* touching, float precision) const {
	/*
	*	find the smallest component of the cross product which ensures that
	*	the other two components are larger, as long as the plane is well defined
	*/
	size_t index = a.cross(b).comp(false);

	/* compute the linear combination across the other two axes */
	Linear r = fLinComb(p, index);

	/* check if the touching property should be validated */
	if (touching != 0)
		*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool num::Plane::inConeX(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentX);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool num::Plane::inConeY(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentY);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool num::Plane::inConeZ(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Vec::Component::ComponentZ);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool num::Plane::inCone(const Vec& p, bool* touching, float precision) const {
	/*
	*	find the smallest component of the cross product which ensures that
	*	the other two components are larger, as long as the plane is well defined
	*/
	size_t index = a.cross(b).comp(false);

	/* compute the linear combination across the other two axes */
	Linear r = fLinComb(p, index);

	/* check if the touching property should be validated */
	if (touching != 0)
		*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool num::Plane::touch(const Vec& p, float precision) const {
	/*
	*	find the smallest component of the cross product which ensures that
	*	the other two components are larger, as long as the plane is well defined
	*/
	size_t index = a.cross(b).comp(false);

	/* compute the linear combination across the other two axes */
	Linear r = fLinComb(p, index);

	/* compute the point on the plane where the given point is expected to be */
	const Vec t = o + a * r.s + b * r.t;

	/* compare the point on the plane with the given point */
	return p.equal(t, precision);
}
bool num::Plane::equal(const Plane& p, float precision) const {
	return p.touch(o, precision) && a.cross(b).parallel(p.normal(), precision);
}
bool num::Plane::identical(const Plane& p, float precision) const {
	return p.o.equal(o, precision) && p.a.equal(a, precision) && p.b.equal(b, precision);
}
num::Vec num::Plane::closest(const Vec& p) const {
	/*
	*	compute this by computing the intersection point of a line with the
	*	planes normal vector as direction vector and the point as origin
	*/
	const Vec crs = a.cross(b);

	/*
	*	E: o + s * x1 + t * x2
	*	G: [this] + f * x0
	*	Solve for f and insert
	*/
	const float f = (o - p).dot(crs) / crs.dot(crs);
	return crs * f;
}
num::Vec num::Plane::steepestX() const {
	/*
	*	If a and b are perpendicular to each other and have the same length,
	*	the vector of steepest ascent can be computed as follows:
	*
	*	steepest = a * a.x + b * b.x
	*/
	Vec t = a.perpendicular(b);

	/* stretch the vectors to the same lengths */
	const float lF = a.lenSquared() / t.lenSquared();
	t = t * std::sqrt(lF);

	/* compute the vector of steepest ascent */
	return a * a.x + t * t.x;
}
num::Vec num::Plane::steepestY() const {
	/*
	*	If a and b are perpendicular to each other and have the same length,
	*	the vector of steepest ascent can be computed as follows:
	*
	*	steepest = a * a.y + b * b.y
	*/
	Vec t = a.perpendicular(b);

	/* stretch the vectors to the same lengths */
	const float lF = a.lenSquared() / t.lenSquared();
	t = t * std::sqrt(lF);

	/* compute the vector of steepest ascent */
	return a * a.y + t * t.y;
}
num::Vec num::Plane::steepestZ() const {
	/*
	*	If a and b are perpendicular to each other and have the same length,
	*	the vector of steepest ascent can be computed as follows:
	*
	*	steepest = a * a.z + b * b.z
	*/
	Vec t = a.perpendicular(b);

	/* stretch the vectors to the same lengths */
	const float lF = a.lenSquared() / t.lenSquared();
	t = t * std::sqrt(lF);

	/* compute the vector of steepest ascent */
	return a * a.z + t * t.z;
}
num::Line num::Plane::intersectPlaneX(float xPlane, bool* invalid, float precision) const {
	/*
	*	order the extent-vectors in order to have the one with the larger x component
	*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
	*/
	const Vec& _x0 = std::abs(a.x) < std::abs(b.x) ? b : a;
	const Vec& _x1 = &_x0 == &a ? b : a;
	if (std::abs(_x0.x) <= precision) {
		if (invalid)
			*invalid = true;
		return Line();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: x = xPlane
	*	F: o + s * x0 + t * x1
	*	Solve for s and insert to have t left over, which will be the free variable in the Line
	*/
	return Line(
		o + Vec(1.0f, _x0.y / _x0.x, _x0.z / _x0.x) * (xPlane - o.x),
		Vec(0.0f, _x1.y - (_x0.y * _x1.x / _x0.x), _x1.z - (_x0.z * _x1.x / _x0.x))
	);
}
num::Line num::Plane::intersectPlaneY(float yPlane, bool* invalid, float precision) const {
	/*
	*	order the extent-vectors in order to have the one with the larger y component
	*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
	*/
	const Vec& _x0 = std::abs(a.y) < std::abs(b.y) ? b : a;
	const Vec& _x1 = &_x0 == &a ? b : a;
	if (std::abs(_x0.y) <= precision) {
		if (invalid)
			*invalid = true;
		return Line();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: y = yPlane
	*	F: o + s * x0 + t * x1
	*	Solve for s and insert to have t left over, which will be the free variable in the Line
	*/
	return Line(
		o + Vec(_x0.x / _x0.y, 1.0f, _x0.z / _x0.y) * (yPlane - o.y),
		Vec(_x1.x - (_x0.x * _x1.y / _x0.y), 0.0f, _x1.z - (_x0.z * _x1.y / _x0.y))
	);
}
num::Line num::Plane::intersectPlaneZ(float zPlane, bool* invalid, float precision) const {
	/*
	*	order the extent-vectors in order to have the one with the larger z component
	*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
	*/
	const Vec& _x0 = std::abs(a.z) < std::abs(b.z) ? b : a;
	const Vec& _x1 = &_x0 == &a ? b : a;
	if (std::abs(_x0.z) <= precision) {
		if (invalid)
			*invalid = true;
		return Line();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: z = zPlane
	*	F: o + s * x0 + t * x1
	*	Solve for s and insert to have t left over, which will be the free variable in the Line
	*/
	return Line(
		o + Vec(_x0.x / _x0.z, _x0.y / _x0.z, 1.0f) * (zPlane - o.z),
		Vec(_x1.x - (_x0.x * _x1.z / _x0.z), _x1.y - (_x0.y * _x1.z / _x0.z), 0.0f)
	);
}
num::Line num::Plane::intersect(const Plane& p, bool* invalid, float precision) const {
	const Vec crs = a.cross(b);

	/* select the extent-vector which is less parallel to the plane and check if the planes are parallel */
	const float dtValue[2] = { crs.dot(p.a), crs.dot(p.b) };
	const float aDtValue[2] = { std::abs(dtValue[0]), std::abs(dtValue[1]) };
	if (aDtValue[0] <= precision && aDtValue[1] <= precision) {
		if (invalid)
			*invalid = true;
		return Line();
	}
	else if (invalid)
		*invalid = false;
	const Vec& _x2 = aDtValue[0] >= aDtValue[1] ? p.a : p.b;
	const Vec& _x3 = aDtValue[0] >= aDtValue[1] ? p.b : p.a;
	const float dt = 1.0f / (aDtValue[0] >= aDtValue[1] ? dtValue[0] : dtValue[1]);

	/*
	*	E: o + s * a + t * b
	*	F: p.o + m * p.a + n * p.b
	*	Solving for m/n yields: m * p.a * (a x b) + n * p.b * (a x b) = (o - p.o) * (a x b)
	*	The more p.a/p.b are parallel to the plane E, the more it gets damped, thus they are sorted by how parallel they are.
	*	By inserting m into the plane F, the following formula is yielded.
	*/
	return Line(
		p.o + _x2 * (((o - p.o).dot(crs)) * dt),
		_x3 - _x2 * (_x3.dot(crs) * dt)
	);
}
num::Plane::Linear num::Plane::intersectf(const Line& l, bool* invalid, float precision) const {
	const Vec crs = a.cross(b);

	/* check if the line and the plane are parallel */
	if (std::abs(crs.dot(l.d)) <= precision) {
		if (invalid)
			*invalid = true;
		return Linear();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: o + s * a + t * b
	*	G: l.o + f * l.d
	*	Solve for f and insert
	*
	*	f = ((o - l.o) * (a x b)) / (l.d * (a x b))
	*	s = ((o - l.o) * (l.d x b)) / (l.d * (a x b))
	*	t = ((o - l.o) * (a x l.d)) / (l.d * (a x b))
	*/
	const float divisor = l.d.dot(crs);
	const float s = (o - l.o).dot(l.d.cross(b)) / divisor;
	const float t = (o - l.o).dot(a.cross(l.d)) / divisor;
	return Linear(s, t);
}
num::Vec num::Plane::intersect(const Line& l, bool* invalid, float precision) const {
	const Vec crs = a.cross(b);

	/* check if the line and the plane are parallel */
	if (std::abs(crs.dot(l.d)) <= precision) {
		if (invalid)
			*invalid = true;
		return Vec();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: o + s * a + t * b
	*	G: l.o + f * l.d
	*	Solve for f and insert
	*
	*	f = ((o - l.o) * (a x b)) / (l.d * (a x b))
	*	s = ((o - l.o) * (l.d x b)) / (l.d * (a x b))
	*	t = ((o - l.o) * (a x l.d)) / (l.d * (a x b))
	*/
	const float a = (o - l.o).dot(crs) / l.d.dot(crs);
	return l.o + l.d * a;
}
num::Plane::Linear num::Plane::linearX(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentX);
}
num::Plane::Linear num::Plane::linearY(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentY);
}
num::Plane::Linear num::Plane::linearZ(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentZ);
}
num::Plane::Linear num::Plane::linear(const Vec& p, bool* touching, float precision) const {
	/*
	*	find the smallest component of the cross product which ensures that
	*	the other two components are larger, as long as the plane is well defined
	*/
	size_t index = a.cross(b).comp(false);

	/* compute the linear combination across the other two axes */
	Linear r = fLinComb(p, index);

	/* check if the touching property should be validated */
	if (touching != 0)
		*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
	return r;
}
