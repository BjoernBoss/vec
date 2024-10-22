#pragma once

#include "num-common.h"
#include "num-vec.h"

namespace num {
	struct Plane {
		num::Vecf o;
		num::Vecf a;
		num::Vecf b;

	public:
		constexpr Plane() = default;
		constexpr Plane(const num::Vecf& a, const num::Vecf& b) : a{ a }, b{ b } {}
		constexpr Plane(const num::Vecf& o, const num::Vecf& a, const num::Vecf& b) : o{ o }, a{ a }, b{ b } {}

	private:
		/* compute the linear combination of the two extent vectors to the point in a plane based on the index axis */
		constexpr num::Linear fLinComb(const num::Vecf& p, size_t index) const {
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
			return num::Linear{ _s, _t };
		}

	public:
		/* create a plane parallel to the Y-Z plane at distance [d] to the origin */
		static constexpr num::Plane AxisX(float d = 0.0f) {
			return num::Plane{ num::Vecf::AxisX(d), num::Vecf::AxisY(), num::Vecf::AxisZ() };
		}

		/* create a plane parallel to the X-Z plane at distance [d] to the origin */
		static constexpr num::Plane AxisY(float d = 1.0f) {
			return num::Plane{ num::Vecf::AxisY(d), num::Vecf::AxisX(), num::Vecf::AxisZ() };
		}

		/* create a plane parallel to the X-Y plane at distance [d] to the origin */
		static constexpr num::Plane AxisZ(float d = 1.0f) {
			return num::Plane{ num::Vecf::AxisZ(d), num::Vecf::AxisX(), num::Vecf::AxisY() };
		}

	public:
		/* compute the plane [this] projected onto the Y-Z plane */
		constexpr num::Plane planeX(float xPlane = 0.0f) const {
			return num::Plane{ o.planeX(xPlane), a.planeX(xPlane), b.planeX(xPlane) };
		}

		/* compute the plane [this] projected onto the X-Z plane */
		constexpr num::Plane planeY(float yPlane = 0.0f) const {
			return num::Plane{ o.planeY(yPlane), a.planeY(yPlane), b.planeY(yPlane) };
		}

		/* compute the plane [this] projected onto the X-Y plane */
		constexpr num::Plane planeZ(float zPlane = 0.0f) const {
			return num::Plane{ o.planeZ(zPlane), a.planeZ(zPlane), b.planeZ(zPlane) };
		}

		/* compute the normal vector of the plane */
		constexpr num::Vecf normal() const {
			return a.cross(b);
		}

		/* compute the area of the triangle created by the plane [this] */
		float area() const {
			/* the magnitude of the vector of the cross product is equivalent to the area of the
			*	parallelogram created by the two component vectors of the cross product */
			return a.cross(b).len() / 2.0f;
		}

		/* compute the area of the triangle created by the plane [this] when projected onto the Y-Z plane */
		constexpr float areaX() const {
			/* the magnitude of the vector of the cross product is equivalent to the area of the
			*	parallelogram created by the two component vectors of the cross product */
			return a.crossX(b) / 2.0f;
		}

		/* compute the area of the triangle created by the plane [this] when projected onto the X-Z plane */
		constexpr float areaY() const {
			/* the magnitude of the vector of the cross product is equivalent to the area of the
			*	parallelogram created by the two component vectors of the cross product */
			return a.crossY(b) / 2.0f;
		}

		/* compute the area of the triangle created by the plane [this] when projected onto the X-Y plane */
		constexpr float areaZ() const {
			/* the magnitude of the vector of the cross product is equivalent to the area of the
			*	parallelogram created by the two component vectors of the cross product */
			return a.crossZ(b) / 2.0f;
		}

		/* compute the center point of the triangle produced by the two extent vectors and the origin */
		constexpr num::Vecf center() const {
			return o + ((a + b) / 3);
		}

		/* compute a point on plane [this] */
		constexpr num::Vecf point(float s, float t) const {
			return o + a * s + b * t;
		}

		/* compute a point on plane [this] */
		constexpr num::Vecf point(const num::Linear& lin) const {
			return o + a * lin.s + b * lin.t;
		}

		/* compute a normalized origin that is close to zero and normalize the directions as well as orient them to be perpendicular */
		num::Plane norm() const {
			/*
			*	compute this by computing the intersection point of a line with the
			*	planes normal vector as direction vector and the zero-vector as its origin
			*/
			const num::Vecf crs = a.cross(b);

			/*
			*	E: o + s * x1 + t * x2
			*	G: [this] + f * x0
			*	Solve for f and insert
			*/
			const float f = o.dot(crs) / crs.dot(crs);
			const num::Vecf _a = a.norm();
			return num::Plane{ crs * f, _a, _a.perpendicular(b).norm() };
		}

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the Y-Z plane */
		constexpr num::Vecf projectX(const num::Vecf& p) const {
			const num::Linear r = fLinComb(p, num::ComponentX);
			return num::Vecf{ o.x + r.s * a.x + r.t * b.x, p.y, p.z };
		}

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the X-Z plane */
		constexpr num::Vecf projectY(const num::Vecf& p) const {
			const num::Linear r = fLinComb(p, num::ComponentY);
			return num::Vecf{ p.x, o.y + r.s * a.y + r.t * b.y, p.z };
		}

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the X-Y plane */
		constexpr num::Vecf projectZ(const num::Vecf& p) const {
			const num::Linear r = fLinComb(p, num::ComponentZ);
			return num::Vecf{ p.x, p.y, o.z + r.s * a.z + r.t * b.z };
		}

		/* compute the vector if [v] is being projected onto the plane */
		constexpr num::Vecf project(const num::Vecf& v) const {
			/*
			*	compute the projection onto the normal of the plane and subtract it from p
			*	as this will result in only the part on the vector within the plane
			*/
			const num::Vecf crs = a.cross(b);
			return v - crs.project(v);
		}

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the Y-Z plane */
		constexpr bool inTriangleX(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentX);
			return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
		}

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the X-Z plane */
		constexpr bool inTriangleY(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentY);
			return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
		}

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the X-Y plane */
		constexpr bool inTriangleZ(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentZ);
			return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
		}

		/* check if [p] lies within the triangle of a and b */
		constexpr bool inTriangle(const num::Vecf& p, bool* touching = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	find the smallest component of the cross product which ensures that
			*	the other two components are larger, as long as the plane is well defined
			*/
			size_t index = a.cross(b).comp(false);

			/* compute the linear combination across the other two axes */
			num::Linear r = fLinComb(p, index);

			/* check if the touching property should be validated */
			if (touching != 0)
				*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
			return r.s >= -precision && r.t >= -precision && (r.s + r.t) <= (1.0f + precision);
		}

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the Y-Z plane */
		constexpr bool inConeX(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentX);
			return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
		}

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Z plane */
		constexpr bool inConeY(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentY);
			return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
		}

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Y plane */
		constexpr bool inConeZ(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			const num::Linear r = fLinComb(p, num::ComponentZ);
			return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
		}

		/* check if [p] lies within the cone of a and b */
		constexpr bool inCone(const num::Vecf& p, bool* touching = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	find the smallest component of the cross product which ensures that
			*	the other two components are larger, as long as the plane is well defined
			*/
			size_t index = a.cross(b).comp(false);

			/* compute the linear combination across the other two axes */
			num::Linear r = fLinComb(p, index);

			/* check if the touching property should be validated */
			if (touching != 0)
				*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
			return r.s >= -precision && r.t >= -precision && (r.s <= 1.0f + precision) && (r.t <= 1.0f + precision);
		}

		/* check if [p] lies on the plane */
		constexpr bool touch(const num::Vecf& p, float precision = num::Precision<float>::Def) const {
			/*
			*	find the smallest component of the cross product which ensures that
			*	the other two components are larger, as long as the plane is well defined
			*/
			size_t index = a.cross(b).comp(false);

			/* compute the linear combination across the other two axes */
			num::Linear r = fLinComb(p, index);

			/* compute the point on the plane where the given point is expected to be */
			const Vecf t = point(r);

			/* compare the point on the plane with the given point */
			return p.equal(t, precision);
		}

		/* check if the plane [p] and plane [this] describe the same plane */
		constexpr bool equal(const num::Plane& p, float precision = num::Precision<float>::Def) const {
			return p.touch(o, precision) && a.cross(b).parallel(p.normal(), precision);
		}

		/* check if the plane [p] and plane [this] describe the identically same plane */
		constexpr bool identical(const num::Plane& p, float precision = num::Precision<float>::Def) const {
			return p.o.equal(o, precision) && p.a.equal(a, precision) && p.b.equal(b, precision);
		}

		/* compute the shortest vector which connects [p] to a point on the plane (automatically perpendicular) */
		constexpr num::Vecf closest(const num::Vecf& p) const {
			/*
			*	compute this by computing the intersection point of a line with the
			*	planes normal vector as direction vector and the point as origin
			*/
			const num::Vecf crs = a.cross(b);

			/*
			*	E: o + s * x1 + t * x2
			*	G: [this] + f * x0
			*	Solve for f and insert
			*/
			const float f = (o - p).dot(crs) / crs.dot(crs);
			return crs * f;
		}

		/* compute the vector of steepest ascent in the plane for the X axis */
		constexpr num::Vecf steepestX() const {
			/*
			*	If a and b are perpendicular to each other and have the same length,
			*	the vector of steepest ascent can be computed as follows:
			*
			*	steepest = a * a.x + b * b.x
			*/
			num::Vecf t = a.perpendicular(b);

			/* stretch the vectors to the same lengths */
			const float lF = a.lenSquared() / t.lenSquared();
			t = t * std::sqrt(lF);

			/* compute the vector of steepest ascent */
			return a * a.x + t * t.x;
		}

		/* compute the vector of steepest ascent in the plane for the Y axis */
		constexpr num::Vecf steepestY() const {
			/*
			*	If a and b are perpendicular to each other and have the same length,
			*	the vector of steepest ascent can be computed as follows:
			*
			*	steepest = a * a.y + b * b.y
			*/
			num::Vecf t = a.perpendicular(b);

			/* stretch the vectors to the same lengths */
			const float lF = a.lenSquared() / t.lenSquared();
			t = t * std::sqrt(lF);

			/* compute the vector of steepest ascent */
			return a * a.y + t * t.y;
		}

		/* compute the vector of steepest ascent in the plane for the Z axis */
		constexpr num::Vecf steepestZ() const {
			/*
			*	If a and b are perpendicular to each other and have the same length,
			*	the vector of steepest ascent can be computed as follows:
			*
			*	steepest = a * a.z + b * b.z
			*/
			num::Vecf t = a.perpendicular(b);

			/* stretch the vectors to the same lengths */
			const float lF = a.lenSquared() / t.lenSquared();
			t = t * std::sqrt(lF);

			/* compute the vector of steepest ascent */
			return a * a.z + t * t.z;
		}

		/* compute the intersecting line between the Y-Z plane at [xPlane] and plane [this] (invalid if parallel: returns null line) */
		constexpr num::Line intersectPlaneX(float xPlane, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	order the extent-vectors in order to have the one with the larger x component
			*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
			*/
			const num::Vecf& _x0 = num::Abs(a.x) < num::Abs(b.x) ? b : a;
			const num::Vecf& _x1 = &_x0 == &a ? b : a;
			if (num::Abs(_x0.x) <= precision) {
				if (invalid)
					*invalid = true;
				return num::Line{};
			}
			else if (invalid)
				*invalid = false;

			/*
			*	E: x = xPlane
			*	F: o + s * x0 + t * x1
			*	Solve for s and insert to have t left over, which will be the free variable in the Line
			*/
			return num::Line{
				o + num::Vecf(1.0f, _x0.y / _x0.x, _x0.z / _x0.x) * (xPlane - o.x),
				num::Vecf(0.0f, _x1.y - (_x0.y * _x1.x / _x0.x), _x1.z - (_x0.z * _x1.x / _x0.x))
			};
		}

		/* compute the intersecting line between the X-Z plane at [yPlane] and plane [this] (invalid if parallel: returns null line) */
		constexpr num::Line intersectPlaneY(float yPlane, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	order the extent-vectors in order to have the one with the larger y component
			*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
			*/
			const num::Vecf& _x0 = num::Abs(a.y) < num::Abs(b.y) ? b : a;
			const num::Vecf& _x1 = &_x0 == &a ? b : a;
			if (num::Abs(_x0.y) <= precision) {
				if (invalid)
					*invalid = true;
				return num::Line{};
			}
			else if (invalid)
				*invalid = false;

			/*
			*	E: y = yPlane
			*	F: o + s * x0 + t * x1
			*	Solve for s and insert to have t left over, which will be the free variable in the Line
			*/
			return num::Line{
				o + num::Vecf(_x0.x / _x0.y, 1.0f, _x0.z / _x0.y) * (yPlane - o.y),
				num::Vecf(_x1.x - (_x0.x * _x1.y / _x0.y), 0.0f, _x1.z - (_x0.z * _x1.y / _x0.y))
			};
		}

		/* compute the intersecting line between the X-Y plane at [zPlane] and plane [this] (invalid if parallel: returns null line) */
		constexpr num::Line intersectPlaneZ(float zPlane, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	order the extent-vectors in order to have the one with the larger z component
			*	at the front and check if the value is not equal to zero (which would else imply that the planes are parallel)
			*/
			const num::Vecf& _x0 = num::Abs(a.z) < num::Abs(b.z) ? b : a;
			const num::Vecf& _x1 = &_x0 == &a ? b : a;
			if (num::Abs(_x0.z) <= precision) {
				if (invalid)
					*invalid = true;
				return num::Line{};
			}
			else if (invalid)
				*invalid = false;

			/*
			*	E: z = zPlane
			*	F: o + s * x0 + t * x1
			*	Solve for s and insert to have t left over, which will be the free variable in the Line
			*/
			return num::Line(
				o + num::Vecf(_x0.x / _x0.z, _x0.y / _x0.z, 1.0f) * (zPlane - o.z),
				num::Vecf(_x1.x - (_x0.x * _x1.z / _x0.z), _x1.y - (_x0.y * _x1.z / _x0.z), 0.0f)
			);
		}

		/* compute the intersection line of the plane [this] and the plane [p] (invalid if parallel: returns null line) */
		constexpr num::Line intersect(const num::Plane& p, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			const num::Vecf crs = a.cross(b);

			/* select the extent-vector which is less parallel to the plane and check if the planes are parallel */
			const float dtValue[2] = { crs.dot(p.a), crs.dot(p.b) };
			const float aDtValue[2] = { num::Abs(dtValue[0]), num::Abs(dtValue[1]) };
			if (aDtValue[0] <= precision && aDtValue[1] <= precision) {
				if (invalid)
					*invalid = true;
				return Line();
			}
			else if (invalid)
				*invalid = false;
			const num::Vecf& _x2 = aDtValue[0] >= aDtValue[1] ? p.a : p.b;
			const num::Vecf& _x3 = aDtValue[0] >= aDtValue[1] ? p.b : p.a;
			const float dt = 1.0f / (aDtValue[0] >= aDtValue[1] ? dtValue[0] : dtValue[1]);

			/*
			*	E: o + s * a + t * b
			*	F: p.o + m * p.a + n * p.b
			*	Solving for m/n yields: m * p.a * (a x b) + n * p.b * (a x b) = (o - p.o) * (a x b)
			*	The more p.a/p.b are parallel to the plane E, the more it gets damped, thus they are sorted by how parallel they are.
			*	By inserting m into the plane F, the following formula is yielded.
			*/
			return num::Line{
				p.o + _x2 * (((o - p.o).dot(crs)) * dt),
				_x3 - _x2 * (_x3.dot(crs) * dt)
			};
		}

		/* compute the intersection factors of the plane [this] and the line [l] (invalid if parallel: returns 0, 0) */
		constexpr num::Linear intersectf(const num::Line& l, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			const num::Vecf crs = a.cross(b);

			/* check if the line and the plane are parallel */
			if (num::Abs(crs.dot(l.d)) <= precision) {
				if (invalid)
					*invalid = true;
				return num::Linear{};
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
			return num::Linear{ s, t };
		}

		/* compute the intersection point of the plane [this] and the line [l] (invalid if parallel: returns null vector) */
		constexpr num::Vecf intersect(const num::Line& l, bool* invalid = 0, float precision = num::Precision<float>::Def) const {
			const num::Vecf crs = a.cross(b);

			/* check if the line and the plane are parallel */
			if (num::Abs(crs.dot(l.d)) <= precision) {
				if (invalid)
					*invalid = true;
				return num::Vecf{};
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

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the Y-Z plane */
		constexpr num::Linear linearX(const num::Vecf& p) const {
			return fLinComb(p, num::ComponentX);
		}

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the X-Z plane */
		constexpr num::Linear linearY(const num::Vecf& p) const {
			return fLinComb(p, num::ComponentY);
		}

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the X-Y plane */
		constexpr num::Linear linearZ(const num::Vecf& p) const {
			return fLinComb(p, num::ComponentZ);
		}

		/* compute the linear combination to reach the point [p] when the point lies on the plane */
		constexpr num::Linear linear(const num::Vecf& p, bool* touching = 0, float precision = num::Precision<float>::Def) const {
			/*
			*	find the smallest component of the cross product which ensures that
			*	the other two components are larger, as long as the plane is well defined
			*/
			size_t index = a.cross(b).comp(false);

			/* compute the linear combination across the other two axes */
			num::Linear r = fLinComb(p, index);

			/* check if the touching property should be validated */
			if (touching != 0)
				*touching = num::Cmp(p.c[index] - o.c[index], r.s * a.c[index] + r.t * b.c[index], precision);
			return r;
		}
	};
}
