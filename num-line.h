#pragma once

#include "num-common.h"
#include "num-vec.h"

namespace num {
	template <class Type>
	struct Line {
		num::Vec<Type> o;
		num::Vec<Type> d;

	public:
		constexpr Line() = default;
		constexpr Line(const num::Vec<Type>& d) : d{ d } {}
		constexpr Line(const num::Vec<Type>& o, const num::Vec<Type>& d) : o{ o }, d{ d } {}

	private:
		/* compute the linear combination of the line [this] and line [l] to intersect based on the two other axes than the index axis */
		constexpr num::Linear<Type> fLinComb(const num::Line<Type>& l, size_t index, bool& parallel, Type precision) const {
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
			const Type divisor = d.c[_0] * l.d.c[_1] - d.c[_1] * l.d.c[_0];

			/* check if the two lines are parallel */
			parallel = (num::Abs(divisor) <= precision);
			if (parallel)
				return num::Linear<Type>{};

			/* compute the linear combination */
			const Type s = (l.d.c[_1] * (l.o.c[_0] - o.c[_0]) - l.d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
			const Type t = (d.c[_1] * (l.o.c[_0] - o.c[_0]) - d.c[_0] * (l.o.c[_1] - o.c[_1])) / divisor;
			return num::Linear<Type>{ s, t };
		}

	public:
		/* create a line along the x axis with the length [l] and its origin the the origin */
		static constexpr num::Line<Type> AxisX(Type l = 1) {
			return num::Line<Type>{ num::Vec<Type>::AxisX(l) };
		}

		/* create a line along the y axis with the length [l] and its origin the the origin */
		static constexpr num::Line<Type> AxisY(Type l = 1) {
			return num::Line<Type>{ num::Vec<Type>::AxisY(l) };
		}

		/* create a line along the z axis with the length [l] and its origin the the origin */
		static constexpr num::Line<Type> AxisZ(Type l = 1) {
			return num::Line<Type>{ num::Vec<Type>::AxisZ(l) };
		}

	public:
		/* compute the line [this] projected onto the Y-Z plane */
		constexpr num::Line<Type> planeX(Type xPlane = 0) const {
			return num::Line<Type>{ o.planeX(xPlane), d.planeX(xPlane) };
		}

		/* compute the line [this] projected onto the X-Z plane */
		constexpr num::Line<Type> planeY(Type yPlane = 0) const {
			return num::Line<Type>{ o.planeY(yPlane), d.planeY(yPlane) };
		}

		/* compute the line [this] projected onto the X-Y plane */
		constexpr num::Line<Type> planeZ(Type zPlane = 0) const {
			return num::Line<Type>{ o.planeZ(zPlane), d.planeZ(zPlane) };
		}

		/* compute a point on line [this] */
		constexpr num::Vec<Type> point(Type t) const {
			return o + d * t;
		}

		/* compute a normalized origin that is close to zero normalize the direction */
		constexpr num::Line<Type> norm() const {
			/*
			*	o - a * d = NewOrigin
			*	-> (o - a * d) * d = 0 (the line and NewOrigin should be perpendicular)
			*	Solve for a
			*/
			const Type a = o.dot(d) / d.dot(d);
			return num::Line<Type>{ o - d * a, d.norm() };
		}

		/* check if [p] lies on line [this] */
		constexpr bool touch(const num::Vec<Type>& p, Type precision = num::Const<Type>::Precision) const {
			/*
			*	E: o + s * d
			*	Set equal to [p] and solve for s and insert
			*/

			/* find the largest component which should not be zero if the line is well defined */
			size_t index = d.comp(true);
			const Type s = (p.c[index] - o.c[index]) / d.c[index];

			/* compute the point on the line where the given point is expected to be */
			const num::Vec<Type> t = point(s);

			/* ensure that the points are equal */
			return p.equal(t, precision);
		}

		/* returns a position along the [this] line where the point [p] lies (result only valid if it lies on the line) */
		constexpr Type find(const num::Vec<Type>& p) const {
			/* extract the largest component of the direction and use it to compute the scaling factor */
			size_t index = d.comp(true);
			return (p.c[index] - o.c[index]) / d.c[index];
		}

		/* check if the line [l] and line [this] describe the same line */
		constexpr bool equal(const num::Line<Type>& l, Type precision = num::Const<Type>::Precision) const {
			return l.touch(o, precision) && l.d.parallel(d, precision);
		}

		/* check if the line [l] and line [this] describe the identically same line */
		constexpr bool identical(const num::Line<Type>& l, Type precision = num::Const<Type>::Precision) const {
			return l.o.equal(o, precision) && l.d.equal(d, precision);
		}

		/* compute the factor for which line [this] reaches the point closest to p (automatically perpendicular) */
		constexpr Type closestf(const num::Vec<Type>& p) const {
			/*
			*	o + a * d = p + v
			*	-> (o + a * d - p) * d = 0 (the line and [p:v] should be perpendicular)
			*	Solve for a
			*/
			return (p - o).dot(d) / d.dot(d);
		}

		/* compute the shortest vector which connects [p] to a point on the line (automatically perpendicular) */
		constexpr num::Vec<Type> closest(const num::Vec<Type>& p) const {
			const Type a = closestf(p);
			return point(a) - p;
		}

		/* compute the shortest line which intersects line [this] and the line [l] */
		constexpr num::Linear<Type> closestf(const num::Line<Type>& l) const {
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
			const num::Vec<Type> v = d.cross(l.d);

			/* check if the lines run in parallel */
			if (v.zero())
				return num::Linear<Type>{ 0, l.closestf(o) };

			/* compute the two scalars */
			const Type tmp = v.dot(v);
			const num::Vec<Type> df = l.o - o;
			const Type s = -df.dot(v.cross(l.d)) / tmp;
			const Type t = -df.dot(v.cross(d)) / tmp;

			/* return the two factors */
			return num::Linear<Type>{ s, t };
		}

		/* compute the shortest line which intersects line [this] and the line [l] */
		constexpr num::Line<Type> closest(const num::Line<Type>& l) const {
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
			const num::Vec<Type> v = d.cross(l.d);

			/* check if the lines run in parallel */
			if (v.zero())
				return num::Line<Type>{ o, l.closest(o) };

			/* compute the two scalars */
			const Type tmp = v.dot(v);
			const num::Vec<Type> df = l.o - o;
			const Type r = df.dot(v) / tmp;
			const Type s = -df.dot(v.cross(l.d)) / tmp;

			/* return the final line */
			return num::Line<Type>{ point(s), v* r };
		}

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the Y-Z plane at [xPlane] (invalid if parallel: returns [this] origin) */
		constexpr Type intersectPlaneXf(Type xPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			/* check if the line and the plane are parallel */
			if (num::Abs(d.x) <= precision) {
				if (invalid)
					*invalid = true;
				return 0;
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

		/* compute the intersection point of line [this] and the Y-Z plane at [xPlane] (invalid if parallel: returns null vector) */
		constexpr num::Vec<Type> intersectPlaneX(Type xPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type a = intersectPlaneXf(xPlane, invalid, precision);
			return point(a);
		}

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the X-Z plane at [yPlane] (invalid if parallel: returns [this] origin) */
		constexpr Type intersectPlaneYf(Type yPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			/* check if the line and the plane are parallel */
			if (num::Abs(d.y) <= precision) {
				if (invalid)
					*invalid = true;
				return 0;
			}
			else if (invalid)
				*invalid = false;

			/*
			*	E: y = yPlane
			*	Line: o + a * d
			*	insert into equation for y and solve for a
			*/
			const Type a = (yPlane - o.y) / d.y;
			return a;
		}

		/* compute the intersection point of line [this] and the X-Z plane at [yPlane] (invalid if parallel: returns null vector) */
		constexpr num::Vec<Type> intersectPlaneY(Type yPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type a = intersectPlaneYf(yPlane, invalid, precision);
			return point(a);
		}

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the X-Y plane at [zPlane] (invalid if parallel: returns [this] origin) */
		constexpr Type intersectPlaneZf(Type zPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			/* check if the line and the plane are parallel */
			if (num::Abs(d.z) <= precision) {
				if (invalid)
					*invalid = true;
				return 0;
			}
			else if (invalid)
				*invalid = false;

			/*
			*	E: z = zPlane
			*	Line: o + a * d
			*	insert into equation for z and solve for a
			*/
			const Type a = (zPlane - o.z) / d.z;
			return a;
		}

		/* compute the intersection point of line [this] and the X-Y plane at [zPlane] (invalid if parallel: returns null vector) */
		constexpr num::Vec<Type> intersectPlaneZ(Type zPlane, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type a = intersectPlaneZf(zPlane, invalid, precision);
			return point(a);
		}

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the X-Y plane (invalid if no intersection point: returns 0, 0) */
		constexpr num::Linear<Type> intersectXf(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			bool parallel = false;

			/* compute the linear combination */
			const num::Linear<Type> lin = fLinComb(l, num::ComponentX, parallel, precision);

			/* update the invalid flag and otherwise return the result */
			if (invalid)
				*invalid = parallel;
			return lin;
		}

		/* compute the intersection point of line [this] and line [l] when viewed in the X-Y plane (invalid if no intersection point: returns [this] origin) */
		constexpr num::Vec<Type> intersectX(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type s = intersectXf(l, invalid, precision).s;
			return point(s);
		}

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the X-Z plane (invalid if no intersection point: returns 0, 0) */
		constexpr num::Linear<Type> intersectYf(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			bool parallel = false;

			/* compute the linear combination */
			const num::Linear<Type> lin = fLinComb(l, num::ComponentY, parallel, precision);

			/* update the invalid flag and otherwise return the result */
			if (invalid)
				*invalid = parallel;
			return lin;
		}

		/* compute the intersection point of line [this] and line [l] when viewed in the X-Z plane (invalid if no intersection point: returns [this] origin) */
		constexpr num::Vec<Type> intersectY(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type s = intersectYf(l, invalid, precision).s;
			return point(s);
		}

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the Y-Z plane (invalid if no intersection point: returns 0, 0) */
		constexpr num::Linear<Type> intersectZf(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			bool parallel = false;

			/* compute the linear combination */
			const num::Linear<Type> lin = fLinComb(l, num::ComponentZ, parallel, precision);

			/* update the invalid flag and otherwise return the result */
			if (invalid)
				*invalid = parallel;
			return lin;
		}

		/* compute the intersection point of line [this] and line [l] when viewed in the Y-Z plane (invalid if no intersection point: returns [this] origin) */
		constexpr num::Vec<Type> intersectZ(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type s = intersectZf(l, invalid, precision).s;
			return point(s);
		}

		/* compute the factor to scale line [this] and line [l] with to intersect the lines (invalid if no intersection point: returns 0, 0) */
		constexpr num::Linear<Type> intersectf(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			/*
			*	Select the axis to compute the combination for by computing the cross product between
			*	the two and then selecting the smallest component which ensures that the other two components
			*	are larger, as long as the lines are well defined
			*/
			const size_t index = d.cross(l.d).comp(false);
			bool parallel = false;

			/* compute the linear combination */
			const num::Linear<Type> lin = fLinComb(l, index, parallel, precision);

			/* check if the lines intersect */
			bool on = false;
			if (!parallel)
				on = num::Cmp(o.c[index] + d.c[index] * lin.s, l.o.c[index] + l.d.c[index] * lin.t, precision);

			/* update the invalid flag and return the result */
			if (invalid)
				*invalid = !on;
			return on ? lin : num::Linear<Type>{};
		}

		/* compute the intersection point of line [this] and the line [l] (invalid if no intersection point: returns [this] origin) */
		constexpr num::Vec<Type> intersect(const num::Line<Type>& l, bool* invalid = 0, Type precision = num::Const<Type>::Precision) const {
			const Type f = intersectf(l, invalid, precision).s;
			return point(f);
		}
	};
}
