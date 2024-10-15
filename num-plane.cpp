#include "vec.h"

num::Plane::Plane() {}
num::Plane::Plane(const Vec& _a, const Vec& _b) : a(_a), b(_b) {}
num::Plane::Plane(const Vec& _o, const Vec& _a, const Vec& _b) : o(_o), a(_a), b(_b) {}
num::Linear num::Plane::fLinComb(const Vec& p, size_t index) const {
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
num::Plane num::Plane::AxisX(float d) {
	return Plane(Vec::AxisX(d), Vec::AxisY(), Vec::AxisZ());
}
num::Plane num::Plane::AxisY(float d) {
	return Plane(Vec::AxisY(d), Vec::AxisX(), Vec::AxisZ());
}
num::Plane num::Plane::AxisZ(float d) {
	return Plane(Vec::AxisZ(d), Vec::AxisX(), Vec::AxisY());
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
float num::Plane::area() const {
	/* the magnitude of the vector of the cross product is equivalent to the area of the
	*	parallelogram created by the two component vectors of the cross product */
	return a.cross(b).len() / 2.0f;
}
float num::Plane::areaX() const {
	/* the magnitude of the vector of the cross product is equivalent to the area of the
	*	parallelogram created by the two component vectors of the cross product */
	return a.crossX(b) / 2.0f;
}
float num::Plane::areaY() const {
	/* the magnitude of the vector of the cross product is equivalent to the area of the
	*	parallelogram created by the two component vectors of the cross product */
	return a.crossY(b) / 2.0f;
}
float num::Plane::areaZ() const {
	/* the magnitude of the vector of the cross product is equivalent to the area of the
	*	parallelogram created by the two component vectors of the cross product */
	return a.crossZ(b) / 2.0f;
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
num::Linear num::Plane::intersectf(const Line& l, bool* invalid, float precision) const {
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
num::Linear num::Plane::linearX(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentX);
}
num::Linear num::Plane::linearY(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentY);
}
num::Linear num::Plane::linearZ(const Vec& p) const {
	return fLinComb(p, Vec::Component::ComponentZ);
}
num::Linear num::Plane::linear(const Vec& p, bool* touching, float precision) const {
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

std::ostream& num::operator<<(std::ostream& out, const Plane& p) {
	return (out << p.o << " -> " << p.a << " | " << p.b);
}
std::wostream& num::operator<<(std::wostream& out, const Plane& p) {
	return (out << p.o << L" -> " << p.a << L" | " << p.b);
}

std::istream& num::operator>>(std::istream& in, Plane& p) {
	char pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != '-' || pad1 != '>' || pad2 != '|')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& num::operator>>(std::wistream& in, Plane& p) {
	wchar_t pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != L'-' || pad1 != L'>' || pad2 != L'|')
		in.setstate(std::ios::failbit);
	return in;

}
