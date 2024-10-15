#include "vec.h"

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
num::Vec num::Vec::AxisX(float l) {
	return Vec(l, 0.0f, 0.0f);
}
num::Vec num::Vec::AxisY(float l) {
	return Vec(0.0f, l, 0.0f);
}
num::Vec num::Vec::AxisZ(float l) {
	return Vec(0.0f, 0.0f, l);
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
	a = num::ToRadian(a);
	const float sa = std::sin(a);
	const float ca = std::cos(a);
	return Vec(
		x * ca - y * sa,
		x * sa + y * ca,
		z);
}
float num::Vec::angleX(const Vec& v) const {
	Vec flat = planeX();
	Vec target = v.planeX();

	/* compute the angle between the two vectors and correct its sign */
	float angle = flat.angle(target);
	return (flat.crossX(target) < 0.0f) ? -angle : angle;
}
float num::Vec::angleY(const Vec& v) const {
	Vec flat = planeY();
	Vec target = v.planeY();

	/* compute the angle between the two vectors and correct its sign */
	float angle = flat.angle(target);
	return (flat.crossY(target) < 0.0f) ? -angle : angle;
}
float num::Vec::angleZ(const Vec& v) const {
	Vec flat = planeZ();
	Vec target = v.planeZ();

	/* compute the angle between the two vectors and correct its sign */
	float angle = flat.angle(target);
	return (flat.crossZ(target) < 0.0f) ? -angle : angle;
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
bool num::Vec::match(const Vec& v, float precision) const {
	return num::Cmp(dot(v), lenSquared(), precision);
}
bool num::Vec::negligibleX(float precision) const {
	return num::Cmp(lenSquared(), planeX().lenSquared(), precision);
}
bool num::Vec::negligibleY(float precision) const {
	return num::Cmp(lenSquared(), planeY().lenSquared(), precision);
}
bool num::Vec::negligibleZ(float precision) const {
	return num::Cmp(lenSquared(), planeZ().lenSquared(), precision);
}
bool num::Vec::isPerpendicular(const Vec& v, float precision) const {
	return num::Zero(dot(v), precision);
}
bool num::Vec::isAcuteAngle(const Vec& v, float precision) const {
	return dot(v) >= -precision;
}
bool num::Vec::isObtuseAngle(const Vec& v, float precision) const {
	return dot(v) <= precision;
}
float num::Vec::projectf(const Vec& v) const {
	return dot(v) / lenSquared();
}
num::Vec num::Vec::project(const Vec& v) const {
	return *this * projectf(v);
}
num::Vec num::Vec::perpendicular(const Vec& v) const {
	return v - project(v);
}
float num::Vec::reachf(const Vec& v) const {
	return v.lenSquared() / dot(v);
}
num::Vec num::Vec::reach(const Vec& v) const {
	return *this * reachf(v);
}
num::Vec num::Vec::passing(const Vec& v) const {
	return v.reach(*this) - *this;
}
float num::Vec::passPointf(const Vec& v) const {
	/* check if the vectors are perpendicular in which case the point will always be considered passed */
	if (num::Zero(dot(v)))
		return 1.0f;

	/* compute the factor required to let this vector reach v and return it if its greater than 1 */
	return std::max(1.0f, reachf(v));
}
num::Vec num::Vec::passPoint(const Vec& v) const {
	return *this * passPointf(v);
}

num::Vec num::operator*(float s, const Vec& v) {
	return v * s;
}

std::ostream& num::operator<<(std::ostream& out, const Vec& v) {
	return (out << v.x << ", " << v.y << ", " << v.z);
}
std::wostream& num::operator<<(std::wostream& out, const Vec& v) {
	return (out << v.x << L", " << v.y << L", " << v.z);
}
std::istream& num::operator>>(std::istream& in, Vec& v) {
	char pad0 = 0, pad1 = 0;
	in >> v.x >> pad0 >> v.y >> pad1 >> v.z;
	if (pad0 != ',' || pad1 != ',')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& num::operator>>(std::wistream& in, Vec& v) {
	wchar_t pad0 = 0, pad1 = 0;
	in >> v.x >> pad0 >> v.y >> pad1 >> v.z;
	if (pad0 != L',' || pad1 != L',')
		in.setstate(std::ios::failbit);
	return in;
}
