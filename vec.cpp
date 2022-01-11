#include "vec.h"

/* implement the vector object */
Vec::Vec() : x(0.0f), y(0.0f), z(0.0f) {}
Vec::Vec(float f) : x(f), y(f), z(f) {}
Vec::Vec(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
Vec Vec::operator+(const Vec& v) const {
	return Vec(x + v.x, y + v.y, z + v.z);
}
Vec Vec::operator-(const Vec& v) const {
	return Vec(x - v.x, y - v.y, z - v.z);
}
Vec Vec::operator-() const {
	return Vec(-x, -y, -z);
}
Vec Vec::operator*(float s) const {
	return Vec(x * s, y * s, z * s);
}
Vec Vec::operator/(float s) const {
	return Vec(x / s, y / s, z / s);
}
Vec& Vec::operator+=(const Vec& v) {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}
Vec& Vec::operator-=(const Vec& v) {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}
Vec& Vec::operator*=(float s) {
	x *= s;
	y *= s;
	z *= s;
	return *this;
}
Vec& Vec::operator/=(float s) {
	x /= s;
	y /= s;
	z /= s;
	return *this;
}
float Vec::dot(const Vec& v) const {
	return v.x * x + v.y * y + v.z * z;
}
float Vec::angle(const Vec& v) const {
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
float Vec::len() const {
	return std::sqrt(x * x + y * y + z * z);
}
float Vec::lenSquared() const {
	return x * x + y * y + z * z;
}
Vec Vec::cross(const Vec& v) const {
	return Vec(
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x);
}
float Vec::crossX(const Vec& v) const {
	return y * v.z - z * v.y;
}
float Vec::crossY(const Vec& v) const {
	return z * v.x - x * v.z;
}
float Vec::crossZ(const Vec& v) const {
	return x * v.y - y * v.x;
}
Vec Vec::norm() const {
	float l = len();
	return Vec(x / l, y / l, z / l);
}
Vec Vec::planeX() const {
	return Vec(0.0f, y, z);
}
Vec Vec::planeY() const {
	return Vec(x, 0.0f, z);
}
Vec Vec::planeZ() const {
	return Vec(x, y, 0.0f);
}
Vec Vec::rotateX(float a) const {
	a = num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x,
		y * ca - z * sa,
		y * sa + z * ca);
}
Vec Vec::rotateY(float a) const {
	a = num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x * ca + z * sa,
		y,
		z * ca - x * sa);
}
Vec Vec::rotateZ(float a) const {
	a = -num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x * ca - y * sa,
		x * sa + y * ca,
		z);
}
size_t Vec::comp(bool largest) const {
	size_t index = 0;

	/* iterate through the components and check if one is larger */
	for (size_t i = 1; i < 3; i++) {
		if (largest ? std::abs(c[index]) >= std::abs(c[i]) : std::abs(c[index]) <= std::abs(c[i]))
			continue;
		index = i;
	}
	return index;
}
Vec::Line Vec::line(const Vec& p) const {
	return Line(*this, p - *this);
}
Vec::Plane Vec::plane(const Vec& p0, const Vec& p1) const {
	return Plane(*this, p0 - *this, p1 - *this);
}
Vec Vec::scale(float l) const {
	float factor = std::sqrt((l * l) / lenSquared());
	return *this * factor;
}
float Vec::scale(const Vec& v) const {
	/* extract the largest component and use it to compute the scaling factor */
	size_t index = v.comp(true);
	return v.c[index] / c[index];
}
bool Vec::parallel(const Vec& v, float precision) const {
	/* check if the vectors are considered zero */
	const float lens[2] = { lenSquared(), v.lenSquared() };
	if (lens[0] <= precision)
		return (lens[1] <= precision);
	else if (lens[1] <= precision)
		return false;

	/* check if the vectors point in the same direction, when scaled and transformed by their sign */
	return num::Cmp(std::abs(dot(v)), std::sqrt(lens[0] * lens[1]), precision);
}
bool Vec::sign(const Vec& v, float precision) const {
	/* check if the vectors are considered zero */
	const float lens[2] = { lenSquared(), v.lenSquared() };
	if (lens[0] <= precision)
		return (lens[1] <= precision);
	else if (lens[1] <= precision)
		return false;

	/* check if the vectors point in the same direction, when scaled */
	return num::Cmp(dot(v), std::sqrt(lens[0] * lens[1]), precision);
}
bool Vec::same(const Vec& v, float precision) const {
	return num::Cmp(dot(v), lenSquared(), precision);
}
bool Vec::identical(const Vec& v, float precision) const {
	return num::Cmp((v - *this).lenSquared(), 0.0f, precision);
}
bool Vec::zeroX(float precision) const {
	return num::Cmp(dot(planeX()), lenSquared(), precision);
}
bool Vec::zeroY(float precision) const {
	return num::Cmp(dot(planeY()), lenSquared(), precision);
}
bool Vec::zeroZ(float precision) const {
	return num::Cmp(dot(planeZ()), lenSquared(), precision);
}
bool Vec::zero(float precision) const {
	return num::Cmp(0.0f, lenSquared(), precision);
}
Vec Vec::perpendicular(const Vec& v) const {
	/*
	*	(v - [this] * a) * [this] = 0
	*	Solve for a and insert
	*/
	const float a = dot(v) / dot(*this);
	return v - *this * a;
}
Vec Vec::project(const Vec& v) const {
	float factor = dot(v) / lenSquared();
	return *this * factor;
}
float Vec::projectFactor(const Vec& v) const {
	return dot(v) / lenSquared();
}

/* implement the line object */
Vec::Line::Linear::Linear() : s(0.0f), t(0.0f) {}
Vec::Line::Linear::Linear(float _s, float _t) : s(_s), t(_t) {}
Vec::Line::Line() {}
Vec::Line::Line(const Vec& _d) : d(_d) {}
Vec::Line::Line(const Vec& _o, const Vec& _d) : o(_o), d(_d) {}
Vec::Line Vec::Line::planeX() const {
	return Line(o.planeX(), d.planeX());
}
Vec::Line Vec::Line::planeY() const {
	return Line(o.planeY(), d.planeY());
}
Vec::Line Vec::Line::planeZ() const {
	return Line(o.planeZ(), d.planeZ());
}
Vec Vec::Line::point(float f) const {
	return o + d * f;
}
Vec::Line Vec::Line::norm() const {
	/*
	*	o - a * d = NewOrigin
	*	-> (o - a * d) * d = 0 (the line and NewOrigin should be perpendicular)
	*	Solve for a
	*/
	const float a = o.dot(d) / d.dot(d);
	return Line(o - d * a, d.norm());
}
bool Vec::Line::touch(const Vec& p, float precision) const {
	/*
	*	E: o + s * d
	*	Set equal to [p] and solve for s and insert
	*/

	/* find the largest component which should not be zero if the line is well defined */
	size_t index = d.comp(true);
	const float s = (p.c[index] - o.c[index]) / d.c[index];

	/* compute the point on the line where the given point is expected to be */
	const Vec t = o + d * s;

	/* compare the point on the line with the given point */
	const float lens[2] = { p.lenSquared(), t.lenSquared() };
	return num::Cmp(p.dot(t), std::sqrt(lens[0] * lens[1]), precision);
}
bool Vec::Line::same(const Line& l, float precision) const {
	return l.touch(o, precision) && l.d.parallel(d, precision);
}
bool Vec::Line::identical(const Line& l, float precision) const {
	return l.o.identical(o, precision) && l.d.identical(d, precision);
}
Vec Vec::Line::closest(const Vec& p) const {
	/*
	*	o + a * d = p + v
	*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
	*	Solve for a
	*/
	const float a = (p - o).dot(d) / d.dot(d);
	return o + d * a - p;
}
float Vec::Line::closestFactor(const Vec& p) const {
	/*
	*	o + a * d = p + v
	*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
	*	Solve for a
	*/
	return (p - o).dot(d) / d.dot(d);
}
Vec::Line Vec::Line::closest(const Line& l) const {
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
	if (v.same(Vec()))
		return Line(o, l.closest(o));

	/* compute the two scalars */
	const float tmp = v.dot(v);
	const Vec df = l.o - o;
	const float r = df.dot(v) / tmp;
	const float s = -df.dot(v.cross(l.d)) / tmp;

	/* compute the final result */
	return Line(o + d * s, v * r);
}
Vec::Line::Linear Vec::Line::closestFactor(const Line& l) const {
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
	if (v.same(Vec()))
		return Linear(0.0f, l.closestFactor(o));

	/* compute the two scalars */
	const float tmp = v.dot(v);
	const Vec df = l.o - o;
	const float s = -df.dot(v.cross(l.d)) / tmp;
	const float t = -df.dot(v.cross(d)) / tmp;

	/* return the two factors */
	return Linear(s, t);
}
Vec Vec::Line::intersectX(float xPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.x) <= precision) {
		if (invalid)
			*invalid = true;
		return Vec();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: x = xPlane
	*	Line: o + a * d
	*	Solve for a and insert
	*/
	const float a = (xPlane - o.x) / d.x;
	return o + d * a;
}
Vec Vec::Line::intersectY(float yPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.y) <= precision) {
		if (invalid)
			*invalid = true;
		return Vec();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: y = yPlane
	*	Line: o + a * d
	*	Solve for a and insert
	*/
	const float a = (yPlane - o.y) / d.y;
	return o + d * a;
}
Vec Vec::Line::intersectZ(float zPlane, bool* invalid, float precision) const {
	/* check if the line and the plane are parallel */
	if (std::abs(d.z) <= precision) {
		if (invalid)
			*invalid = true;
		return Vec();
	}
	else if (invalid)
		*invalid = false;

	/*
	*	E: z = zPlane
	*	Line: o + a * d
	*	Solve for a and insert
	*/
	const float a = (zPlane - o.z) / d.z;
	return o + d * a;
}
float Vec::Line::intersectXFactor(float xPlane, bool* invalid, float precision) const {
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
	*	Solve for a and insert
	*/
	const float a = (xPlane - o.x) / d.x;
	return a;
}
float Vec::Line::intersectYFactor(float yPlane, bool* invalid, float precision) const {
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
	*	Solve for a and insert
	*/
	const float a = (yPlane - o.y) / d.y;
	return a;
}
float Vec::Line::intersectZFactor(float zPlane, bool* invalid, float precision) const {
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
	*	Solve for a and insert
	*/
	const float a = (zPlane - o.z) / d.z;
	return a;
}
Vec Vec::Line::intersect(const Line& l, bool* invalid, float precision) const {
	/*
	*	E: o + s * d
	*	F: l.o + t * l.d
	*
	*	Solve for s and insert yields: s = l.d.y * (l.o.x - o.x) - l.d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*								   t = d.y * (l.o.x - o.x) - d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*	(same for x-z/y-z)
	*/

	/*
	*	Select the axis to compute the combination for and compute the divisor.
	*	Select the axis by computing the cross product between the two and then selecting the smallest two components of it.
	*/
	const size_t i = d.cross(l.d).comp(true);
	const size_t _0 = (i + 1) % 3;
	const size_t _1 = (i + 2) % 3;
	const float divisor = d.c[_0] * l.d.c[_1] - d.c[_1] * l.d.c[_0];

	/* check if the two lines are parallel */
	Vec pt;
	bool on = false;
	if (std::abs(divisor) <= precision) {
		pt = o;
		on = l.touch(pt, precision);
	}

	/* compute the points and check if they are equal */
	else {
		const float s = (l.d.c[_1] * (l.o.c[_0] - o.c[_0]) - l.d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
		const float t = (d.c[_1] * (l.o.c[_0] - o.c[_0]) - d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
		pt = o + d * s;
		on = pt.same(l.o + l.d * t);
	}

	/* return the point if it is on the line */
	if (invalid)
		*invalid = !on;
	return on ? pt : Vec();
}
Vec::Line::Linear Vec::Line::intersectFactor(const Line& l, bool* invalid, float precision) const {
	/*
	*	E: o + s * d
	*	F: l.o + t * l.d
	*
	*	Solve for s and insert yields: s = l.d.y * (l.o.x - o.x) - l.d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*								   t = d.y * (l.o.x - o.x) - d.x * (l.o.y - o.y) / (d.x * l.d.y - d.y * l.d.x)
	*	(same for x-z/y-z)
	*/

	/*
	*	Select the axis to compute the combination for and compute the divisor.
	*	Select the axis by computing the cross product between the two and then selecting the smallest two components of it.
	*/
	const size_t i = d.cross(l.d).comp(true);
	const size_t _0 = (i + 1) % 3;
	const size_t _1 = (i + 2) % 3;
	const float divisor = d.c[_0] * l.d.c[_1] - d.c[_1] * l.d.c[_0];

	/* check if the two lines are parallel */
	bool on = false;
	Linear lin;
	if (std::abs(divisor) <= precision)
		on = l.touch(o, precision);

	/* compute the points and check if they are equal */
	else {
		lin.s = (l.d.c[_1] * (l.o.c[_0] - o.c[_0]) - l.d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
		lin.t = (d.c[_1] * (l.o.c[_0] - o.c[_0]) - d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
		Vec pt = o + d * lin.s;
		on = pt.same(l.o + l.d * lin.t);
	}

	/* return the point if it is on the line */
	if (invalid)
		*invalid = !on;
	return on ? lin : Linear();
}

/* implement the plane object */
Vec::Plane::Linear::Linear() : s(0.0f), t(0.0f) {}
Vec::Plane::Linear::Linear(float _s, float _t) : s(_s), t(_t) {}
Vec::Plane::Plane() {}
Vec::Plane::Plane(const Vec& _a, const Vec& _b) : a(_a), b(_b) {}
Vec::Plane::Plane(const Vec& _o, const Vec& _a, const Vec& _b) : o(_o), a(_a), b(_b) {}
Vec::Plane::Linear Vec::Plane::fLinComb(const Vec& p, size_t index) const {
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
Vec::Plane Vec::Plane::planeX() const {
	return Plane(o.planeX(), a.planeX(), b.planeX());
}
Vec::Plane Vec::Plane::planeY() const {
	return Plane(o.planeY(), a.planeY(), b.planeY());
}
Vec::Plane Vec::Plane::planeZ() const {
	return Plane(o.planeZ(), a.planeZ(), b.planeZ());
}
Vec Vec::Plane::normal() const {
	return a.cross(b);
}
Vec Vec::Plane::center() const {
	return o + ((a + b) / 3);
}
Vec Vec::Plane::point(float s, float t) const {
	return o + a * s + b * t;
}
Vec Vec::Plane::point(const Linear& lin) const {
	return o + a * lin.s + b * lin.t;
}
Vec::Plane Vec::Plane::norm() const {
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
Vec Vec::Plane::projectX(const Vec& p) const {
	const Linear r = fLinComb(p, Component::ComponentX);
	return Vec(o.x + r.s * a.x + r.t * b.x, p.y, p.z);
}
Vec Vec::Plane::projectY(const Vec& p) const {
	const Linear r = fLinComb(p, Component::ComponentY);
	return Vec(p.x, o.y + r.s * a.y + r.t * b.y, p.z);
}
Vec Vec::Plane::projectZ(const Vec& p) const {
	const Linear r = fLinComb(p, Component::ComponentZ);
	return Vec(p.x, p.y, o.z + r.s * a.z + r.t * b.z);
}
Vec Vec::Plane::project(const Vec& v) const {
	/*
	*	compute the projection onto the normal of the plane and subtract it from p
	*	as this will result in only the part on the vector within the plane
	*/
	const Vec crs = a.cross(b);
	return v - crs.project(v);
}
bool Vec::Plane::inTriangleX(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentX);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool Vec::Plane::inTriangleY(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentY);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool Vec::Plane::inTriangleZ(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentZ);
	return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
}
bool Vec::Plane::inTriangle(const Vec& p, bool* touching, float precision) const {
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
bool Vec::Plane::inConeX(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentX);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool Vec::Plane::inConeY(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentY);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool Vec::Plane::inConeZ(const Vec& p, float precision) const {
	const Linear r = fLinComb(p, Component::ComponentZ);
	return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
}
bool Vec::Plane::inCone(const Vec& p, bool* touching, float precision) const {
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
bool Vec::Plane::touch(const Vec& p, float precision) const {
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
	const float lens[2] = { p.lenSquared(), t.lenSquared() };
	return num::Cmp(p.dot(t), std::sqrt(lens[0] * lens[1]), precision);
}
bool Vec::Plane::same(const Plane& p, float precision) const {
	return p.touch(o, precision) && a.cross(b).parallel(p.normal(), precision);
}
bool Vec::Plane::identical(const Plane& p, float precision) const {
	return p.o.identical(o, precision) && p.a.identical(a, precision) && p.b.identical(b, precision);
}
Vec Vec::Plane::closest(const Vec& p) const {
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
Vec Vec::Plane::steepestX() const {
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
Vec Vec::Plane::steepestY() const {
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
Vec Vec::Plane::steepestZ() const {
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
Vec::Line Vec::Plane::intersectX(float xPlane, bool* invalid, float precision) const {
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
Vec::Line Vec::Plane::intersectY(float yPlane, bool* invalid, float precision) const {
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
Vec::Line Vec::Plane::intersectZ(float zPlane, bool* invalid, float precision) const {
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
Vec::Line Vec::Plane::intersect(const Plane& p, bool* invalid, float precision) const {
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
Vec Vec::Plane::intersect(const Line& l, bool* invalid, float precision) const {
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
Vec::Plane::Linear Vec::Plane::linearX(const Vec& p) const {
	return fLinComb(p, Component::ComponentX);
}
Vec::Plane::Linear Vec::Plane::linearY(const Vec& p) const {
	return fLinComb(p, Component::ComponentY);
}
Vec::Plane::Linear Vec::Plane::linearZ(const Vec& p) const {
	return fLinComb(p, Component::ComponentZ);
}
Vec::Plane::Linear Vec::Plane::linear(const Vec& p, bool* touching, float precision) const {
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
