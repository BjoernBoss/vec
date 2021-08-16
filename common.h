#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>

/* topview: (BO3)
*	       +x
*	       ^
*	       |
*	+z <---+---
*	       |
*
*/

/* define the float comparison data */
static constexpr float FloatCmpPrecision = 0.00001f;
static bool FloatCompare(float a, float b, float p = FloatCmpPrecision) {
	if (std::isnan(a) || std::isnan(b))
		return false;
	if (a == 0.0f)
		return std::abs(b) <= p;
	if (b == 0.0f)
		return std::abs(a) <= p;
	const float _a = std::abs(a);
	const float _b = std::abs(b);
	return std::abs(a - b) <= std::min(_a, _b) * p;
}

/* define the vector class */
struct Vec {
public:
	enum Component : uint8_t {
		ComponentX = 0,
		ComponentY = 2,
		ComponentZ = 1
	};

private:
	static constexpr float ConstPi = 3.1415926536f;

public:
	static constexpr float toRadian(float deg) {
		return (deg * ConstPi) / 180.0f;
	}
	static constexpr float toDegree(float deg) {
		return (deg * 180.0f) / ConstPi;
	}

public:
	union {
		struct {
			float x;
			float z;
			float y;
		};
		float c[3];
	};

public:
	Vec() : x(0.0f), y(0.0f), z(0.0f) {}
	Vec(Vec&& v) noexcept : x(v.x), y(v.y), z(v.z) {}
	Vec(const Vec& v) : x(v.x), y(v.y), z(v.z) {}
	Vec(float f) : x(f), y(f), z(f) {}
	Vec(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	Vec& operator=(Vec&& v) noexcept {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	Vec& operator=(const Vec& v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	Vec operator+(const Vec& v) const {
		return Vec(x + v.x, y + v.y, z + v.z);
	}
	Vec operator-(const Vec& v) const {
		return Vec(x - v.x, y - v.y, z - v.z);
	}
	Vec operator-() const {
		return Vec(-x, -y, -z);
	}
	Vec operator*(float s) const {
		return Vec(x * s, y * s, z * s);
	}
	Vec operator/(float s) const {
		return Vec(x / s, y / s, z / s);
	}
	Vec& operator+=(const Vec& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vec& operator-=(const Vec& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	Vec& operator*=(float s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	Vec& operator/=(float s) {
		x /= s;
		y /= s;
		z /= s;
		return *this;
	}

public:
	float dot(const Vec& v) const {
		return v.x * x + v.y * y + v.z * z;
	}
	static float angle(float x, float y) {
		float deg = toDegree(std::atan2(x, y));
		if (deg < 0)
			deg += 360.0f;
		return deg;
	}
	float angle(const Vec& v) const {
		float dotProd = dot(v);
		float lenProd = std::sqrt(lenSquared() * v.lenSquared());
		float frac = dotProd / lenProd;

		/* check if the angle reaches the deges where the uncertainty of floats
			might lead to values outside of the scope of arccos */
		if (frac >= 1.0f)
			return 0.0f;
		else if (frac <= -1.0f)
			return 180.0f;
		return  toDegree(std::acos(frac));
	}
	float len() const {
		return std::sqrt(x * x + y * y + z * z);
	}
	float lenSquared() const {
		return x * x + y * y + z * z;
	}
	Vec cross(const Vec& v) const {
		return Vec(
			y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x);
	}
	float crossX(const Vec& v) const {
		return y * v.z - z * v.y;
	}
	float crossY(const Vec& v) const {
		return z * v.x - x * v.z;
	}
	float crossZ(const Vec& v) const {
		return x * v.y - y * v.x;
	}
	Vec norm() const {
		float l = len();
		return Vec(x / l, y / l, z / l);
	}
	Vec planeX() const {
		return Vec(0.0f, y, z);
	}
	Vec planeY() const {
		return Vec(x, 0.0f, z);
	}
	Vec planeZ() const {
		return Vec(x, y, 0.0f);
	}
	Vec rotateX(float a) const {
		a = toRadian(a);
		const float sa = std::sin(a);
		const float ca = std::cos(a);
		return Vec(
			x,
			y * ca - z * sa,
			y * sa + z * ca);
	}
	Vec rotateY(float a) const {
		a = toRadian(a);
		const float sa = std::sin(a);
		const float ca = std::cos(a);
		return Vec(
			x * ca + z * sa,
			y,
			z * ca - x * sa);
	}
	Vec rotateZ(float a) const {
		a = -toRadian(a);
		const float sa = std::sin(a);
		const float ca = std::cos(a);
		return Vec(
			x * ca - y * sa,
			x * sa + y * ca,
			z);
	}
	size_t comp(bool largest) const {
		size_t index = 0;

		/* iterate through the components and check if one is larger */
		for (size_t i = 1; i < 3; i++) {
			if (largest ? std::abs(c[index]) >= std::abs(c[i]) : std::abs(c[index]) <= std::abs(c[i]))
				continue;
			index = i;
		}
		return index;
	}

public:
	struct Line;
	struct Plane;

public:
	/* construct the line [this:(p-this)] */
	Line line(const Vec& p) const;

	/* construct the plane [this:(p0-this):(p1-this)] */
	Plane plane(const Vec& p0, const Vec& p1) const;

	/* check if [this] and [v] are linear combinations of each other */
	bool parallel(const Vec& v, float precision = FloatCmpPrecision) const {
		/* check if the vectors are considered zero */
		const float lens[2] = { lenSquared(), v.lenSquared() };
		if (lens[0] <= precision)
			return (lens[1] <= precision);
		else if (lens[1] <= precision)
			return false;

		/* check if the vectors point in the same direction, when scaled and transformed by their sign */
		return FloatCompare(std::abs(dot(v)), std::sqrt(lens[0] * lens[1]), precision);
	}

	/* check if [this] and [v] are equal */
	bool equal(const Vec& v, float precision = FloatCmpPrecision) const {
		return FloatCompare(dot(v), lenSquared(), precision);
	}

	/* check if the x-component of this vector is negligible */
	bool zeroX(float precision = FloatCmpPrecision) const {
		return FloatCompare(dot(planeX()), lenSquared(), precision);
	}

	/* check if the y-component of this vector is negligible */
	bool zeroY(float precision = FloatCmpPrecision) const {
		return FloatCompare(dot(planeY()), lenSquared(), precision);
	}

	/* check if the z-component of this vector is negligible */
	bool zeroZ(float precision = FloatCmpPrecision) const {
		return FloatCompare(dot(planeZ()), lenSquared(), precision);
	}

	/* compute a vector which is a linear combination of [this] and [v] but is perpendicular to [this] */
	Vec perpendicular(const Vec& v) const {
		/*
		*	(v - [this] * a) * [this] = 0
		*	Solve for a and insert
		*/
		const float a = dot(v) / dot(*this);
		return v - *this * a;
	}
};
struct Vec::Line {
	Vec o;
	Vec d;

public:
	Line() {}
	Line(Line&& l) noexcept : o(l.o), d(l.d) {}
	Line(const Line& l) : o(l.o), d(l.d) {}
	Line(const Vec& _d) : d(_d) {}
	Line(const Vec& _o, const Vec& _d) : o(_o), d(_d) {}

public:
	struct Linear {
		float s;
		float t;

	public:
		Linear() : s(0.0f), t(0.0f) {}
		Linear(float _s, float _t) : s(_s), t(_t) {}
	};

public:
	Line& operator=(Line&& l) noexcept {
		o = l.o;
		d = l.d;
		return *this;
	}
	Line& operator=(const Line& l) {
		o = l.o;
		d = l.d;
		return *this;
	}

public:
	Line planeX() const {
		return Line(o.planeX(), d.planeX());
	}
	Line planeY() const {
		return Line(o.planeY(), d.planeY());
	}
	Line planeZ() const {
		return Line(o.planeZ(), d.planeZ());
	}

public:
	/* compute a point on this line */
	Vec point(float f) const {
		return o + d * f;
	}

	/* compute a normalized origin that is close to zero normalize the direction */
	Line norm() const {
		/*
		*	o - a * d = NewOrigin
		*	-> (o - a * d) * d = 0 (the line and NewOrigin should be perpendicular)
		*	Solve for a
		*/
		const float a = o.dot(d) / d.dot(d);
		return Line(o - d * a, d.norm());
	}

	/* check if [p] lies on the line */
	bool touch(const Vec& p, float precision = FloatCmpPrecision) const {
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
		return FloatCompare(p.dot(t), std::sqrt(lens[0] * lens[1]), precision);
	}

	/* check if the line [l] and this line describe the same line */
	bool same(const Line& l, float precision = FloatCmpPrecision) const {
		return l.touch(o, precision) && l.d.parallel(d, precision);
	}

	/* compute the shortest vector which connects [p] to a point on the line (automatically perpendicular) */
	Vec closest(const Vec& p) const {
		/*
		*	o + a * d = p + v
		*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
		*	Solve for a
		*/
		const float a = (p - o).dot(d) / d.dot(d);
		return o + d * a - p;
	}

	/* compute the factor for which this line reaches the point closest to p (automatically perpendicular) */
	float closestFactor(const Vec& p) const {
		/*
		*	o + a * d = p + v
		*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
		*	Solve for a
		*/
		return (p - o).dot(d) / d.dot(d);
	}

	/* compute the shortest line which intersects this line and the line [l] */
	Line closest(const Line& l) const {
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
		if (v.equal(Vec()))
			return Line(o, l.closest(o));

		/* compute the two scalars */
		const float tmp = v.dot(v);
		const Vec df = l.o - o;
		const float r = df.dot(v) / tmp;
		const float s = -df.dot(v.cross(l.d)) / tmp;

		/* compute the final result */
		return Line(o + d * s, v * r);
	}

	/* compute the shortest line which intersects this line and the line [l] */
	Linear closestFactor(const Line& l) const {
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
		if (v.equal(Vec()))
			return Linear(0.0f, l.closestFactor(o));

		/* compute the two scalars */
		const float tmp = v.dot(v);
		const Vec df = l.o - o;
		const float s = -df.dot(v.cross(l.d)) / tmp;
		const float t = -df.dot(v.cross(d)) / tmp;

		/* return the two factors */
		return Linear(s, t);
	}

	/* compute the intersection point of this line and the Y-Z plane at [xPlane] */
	Vec intersectX(float xPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersection point of the line [this:x0] and the X-Z plane at [yPlane] */
	Vec intersectY(float yPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersection point of the line [this:x0] and the X-Y plane at [zPlane] */
	Vec intersectZ(float zPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersection point of this line and the line [l] */
	Vec intersect(const Line& l, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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
			Vec pt = o + d * s;
			on = pt.equal(l.o + l.d * t);
		}

		/* return the point if it is on the line */
		if (invalid)
			*invalid = !on;
		return on ? pt : Vec();
	}

	/* compute the factor to scale this lines direction and [l]'s direction with to intersect the lines */
	Linear intersectFactor(const Line& l, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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
			on = pt.equal(l.o + l.d * lin.t);
		}

		/* return the point if it is on the line */
		if (invalid)
			*invalid = !on;
		return on ? lin : Linear();
	}
};
struct Vec::Plane {
	Vec o;
	Vec a;
	Vec b;

public:
	Plane() {}
	Plane(Plane&& p) noexcept : o(p.o), a(p.a), b(p.b) {}
	Plane(const Plane& p) : o(p.o), a(p.a), b(p.b) {}
	Plane(const Vec& _a, const Vec& _b) : a(_a), b(_b) {}
	Plane(const Vec& _o, const Vec& _a, const Vec& _b) : o(_o), a(_a), b(_b) {}

public:
	Plane& operator=(Plane&& p) noexcept {
		o = p.o;
		a = p.a;
		b = p.b;
		return *this;
	}
	Plane& operator=(const Plane& p) {
		o = p.o;
		a = p.a;
		b = p.b;
		return *this;
	}

public:
	/*
	*	compute the linear combination of the two extent vectors to the point [this] in a plane (here in X-Y plane)
	*
	*	p = o + s * a + t * b
	*	s = ((p.x - o.x) * b.y - (p.y - o.y) * b.x) / (a.x * b.y - a.y * b.x)
	*	t = (a.x * (p.y - o.y) - a.y * (p.x - o.x)) / (a.x * b.y - a.y * b.x)
	*/
	struct Linear {
		float s;
		float t;

	public:
		Linear() : s(0.0f), t(0.0f) {}
		Linear(float _s, float _t) : s(_s), t(_t) {}
	};

private:
	Linear fLinCombX(const Vec& p) const {
		const float divisor = a.y * b.z - a.z * b.y;

		const float _y = p.y - o.y;
		const float _z = p.z - o.z;
		const float _s = (_y * b.z - _z * b.y) / divisor;
		const float _t = (a.y * _z - a.z * _y) / divisor;
		return Linear(_s, _t);
	}
	Linear fLinCombY(const Vec& p) const {
		const float divisor = a.x * b.z - a.z * b.x;

		const float _x = p.x - o.x;
		const float _z = p.z - o.z;
		const float _s = (_x * b.z - _z * b.x) / divisor;
		const float _t = (a.x * _z - a.z * _x) / divisor;
		return Linear(_s, _t);
	}
	Linear fLinCombZ(const Vec& p) const {
		const float divisor = a.x * b.y - a.y * b.x;

		const float _x = p.x - o.x;
		const float _y = p.y - o.y;
		const float _s = (_x * b.y - _y * b.x) / divisor;
		const float _t = (a.x * _y - a.y * _x) / divisor;
		return Linear(_s, _t);
	}

public:
	Plane planeX() const {
		return Plane(o.planeX(), a.planeX(), b.planeX());
	}
	Plane planeY() const {
		return Plane(o.planeY(), a.planeY(), b.planeY());
	}
	Plane planeZ() const {
		return Plane(o.planeZ(), a.planeZ(), b.planeZ());
	}

public:
	/* compute the normal vector of the plane */
	Vec normal() const {
		return a.cross(b);
	}

	/* compute the center point of the triangle produced by the two extent vectors and the origin */
	Vec center() const {
		return o + ((a + b) / 3);
	}

	/* compute a point on this plane */
	Vec point(float s, float t) const {
		return o + a * s + b * t;
	}

	/* compute a point on this plane */
	Vec point(const Linear& lin) const {
		return o + a * lin.s + b * lin.t;
	}

	/* compute a normalized origin that is close to zero and normalize the directions as well as orient them to be perpendicular */
	Plane norm() const {
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

	/* compute the vector of if [p] is being projected onto the plane viewed orthogonally from the Y-Z plane */
	Vec projectX(const Vec& p) const {
		const Linear r = fLinCombX(p);
		return Vec(o.x + r.s * a.x + r.t * b.x, p.y, p.z);
	}

	/* compute the vector of if [p] is being projected onto the plane viewed orthogonally from the X-Z plane */
	Vec projectY(const Vec& p) const {
		const Linear r = fLinCombY(p);
		return Vec(p.x, o.y + r.s * a.y + r.t * b.y, p.z);
	}

	/* compute the vector of if [p] is being projected onto the plane viewed orthogonally from the X-Y plane */
	Vec projectZ(const Vec& p) const {
		const Linear r = fLinCombZ(p);
		return Vec(p.x, p.y, o.z + r.s * a.z + r.t * b.z);
	}

	/* check if [p] lies within the triangle when projected orthogonally onto the Y-Z plane */
	bool inTriangleX(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombX(p);
		return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
	}

	/* check if [p] lies within the triangle when projected orthogonally onto the X-Z plane */
	bool inTriangleY(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombY(p);
		return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
	}

	/* check if [p] lies within the triangle when projected orthogonally onto the X-Y plane */
	bool inTriangleZ(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombZ(p);
		return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
	}

	/* check if [p] lies within the cone of a and b when projected orthogonally onto the Y-Z plane */
	float inConeX(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombX(p);
		return (r.s >= -precision && r.t >= -precision) ? r.s + r.t : -1.0f;
	}

	/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Z plane */
	float inConeY(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombY(p);
		return (r.s >= -precision && r.t >= -precision) ? r.s + r.t : -1.0f;
	}

	/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Y plane */
	float inConeZ(const Vec& p, float precision = FloatCmpPrecision) const {
		const Linear r = fLinCombZ(p);
		return (r.s >= -precision && r.t >= -precision) ? r.s + r.t : -1.0f;
	}

	/* check if [p] lies on the plane */
	bool touch(const Vec& p, float precision = FloatCmpPrecision) const {
		/*
		*	find the smallest component of the cross product which ensures that
		*	the other two components are larger, as long as the plane is well defined
		*/
		size_t index = a.cross(b).comp(false);

		/* compute the linear combination across the other two axes */
		Linear r;
		if (index == Component::ComponentX)
			r = fLinCombX(p);
		else if (index == Component::ComponentY)
			r = fLinCombY(p);
		else
			r = fLinCombZ(p);

		/* compute the point on the plane where the given point is expected to be */
		const Vec t = o + a * r.s + b * r.t;

		/* compare the point on the plane with the given point */
		const float lens[2] = { p.lenSquared(), t.lenSquared() };
		return FloatCompare(p.dot(t), std::sqrt(lens[0] * lens[1]), precision);
	}

	/* check if the plane [p] and this plane describe the same plane */
	bool same(const Plane& p, float precision = FloatCmpPrecision) const {
		return p.touch(o, precision) && a.cross(b).parallel(p.normal(), precision);
	}

	/* compute the shortest vector which connects [p] to a point on the plane (automatically perpendicular) */
	Vec closest(const Vec& p) const {
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

	/* comopute the vector of steepest ascent in the plane for the X axis */
	Vec steepestX() const {
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

	/* comopute the vector of steepest ascent in the plane for the Y axis */
	Vec steepestY() const {
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

	/* comopute the vector of steepest ascent in the plane for the Z axis */
	Vec steepestZ() const {
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

	/* compute the intersecting line between the Y-Z plane at [xPlane] and this plane */
	Line intersectX(float xPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersecting line between the X-Z plane at [yPlane] and this plane */
	Line intersectY(float yPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersecting line between the X-Y plane at [zPlane] and this plane */
	Line intersectZ(float zPlane, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersection line of the this and the plane [p] */
	Line intersect(const Plane& p, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the intersection point of the line [l] and this plane */
	Vec intersect(const Line& l, bool* invalid = 0, float precision = FloatCmpPrecision) const {
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

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the Y-Z plane */
	Linear linearX(const Vec& p) const {
		return fLinCombX(p);
	}

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the X-Z plane */
	Linear linearY(const Vec& p) const {
		return fLinCombY(p);
	}

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the X-Y plane */
	Linear linearZ(const Vec& p) const {
		return fLinCombZ(p);
	}
};
