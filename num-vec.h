/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright (c) 2024 Bjoern Boss Henrichsen */
#pragma once

#include "num-common.h"

namespace num {
	template <std::floating_point> struct Line;
	template <std::floating_point> struct Plane;

	/* defines the component layout in memory */
	enum Component : uint8_t {
		ComponentX = 0,
		ComponentY = 1,
		ComponentZ = 2,
	};

	template <std::floating_point Type>
	struct Vec {
	public:
		union {
			struct {
				Type x;
				Type y;
				Type z;
			};
			Type c[3];
		};

	public:
		constexpr Vec() : x{ 0 }, y{ 0 }, z{ 0 } {}
		constexpr Vec(Type f) : x{ f }, y{ f }, z{ f } {}
		constexpr Vec(Type x, Type y, Type z) : x{ x }, y{ y }, z{ z } {}

	public:
		constexpr num::Vec<Type> operator+(const num::Vec<Type>& v) const {
			return num::Vec<Type>{ x + v.x, y + v.y, z + v.z };
		}
		constexpr num::Vec<Type> operator-(const num::Vec<Type>& v) const {
			return num::Vec<Type>{ x - v.x, y - v.y, z - v.z };
		}
		constexpr num::Vec<Type> operator-() const {
			return num::Vec<Type>{ -x, -y, -z };
		}
		constexpr num::Vec<Type> operator*(Type s) const {
			return num::Vec<Type>{ x* s, y* s, z* s };
		}
		constexpr num::Vec<Type> operator/(Type s) const {
			return num::Vec<Type>{ x / s, y / s, z / s };
		}
		constexpr num::Vec<Type>& operator+=(const num::Vec<Type>& v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		constexpr num::Vec<Type>& operator-=(const num::Vec<Type>& v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		constexpr num::Vec<Type>& operator*=(Type s) {
			x *= s;
			y *= s;
			z *= s;
			return *this;
		}
		constexpr num::Vec<Type>& operator/=(Type s) {
			x /= s;
			y /= s;
			z /= s;
			return *this;
		}
		constexpr bool operator==(const num::Vec<Type>& v) const {
			return identical(v);
		}
		constexpr bool operator!=(const num::Vec<Type>& v) const {
			return !(*this == v);
		}

	public:
		/* create a vector on the x axis with the length [l] */
		static constexpr num::Vec<Type> AxisX(Type l = 1) {
			return num::Vec<Type>{ l, 0, 0 };
		}

		/* create a vector on the y axis with the length [l] */
		static constexpr num::Vec<Type> AxisY(Type l = 1) {
			return num::Vec<Type>{ 0, l, 0 };
		}

		/* create a vector on the z axis with the length [l] */
		static constexpr num::Vec<Type> AxisZ(Type l = 1) {
			return num::Vec<Type>{ 0, 0, l };
		}

	public:
		/* compute the dot product [this] * [v] */
		constexpr Type dot(const num::Vec<Type>& v) const {
			return v.x * x + v.y * y + v.z * z;
		}

		/* compute the angle between the vector [this] and [v] [0; 180] */
		constexpr Type angle(const num::Vec<Type>& v) const {
			Type dotProd = dot(v);
			Type lenProd = std::sqrt(lenSquared() * v.lenSquared());
			Type frac = dotProd / lenProd;

			/* check if the angle reaches the degrees where the uncertainty of floats
				might lead to values outside of the scope of arccos */
			if (frac >= 1)
				return 0;
			else if (frac <= -1)
				return 180;
			return num::ToDegree(std::acos(frac));
		}

		/* compute the squared length of the vector */
		constexpr Type lenSquared() const {
			return dot(*this);
		}

		/* compute the length of the vector */
		constexpr Type len() const {
			return std::sqrt(lenSquared());
		}

		/* compute the the cross product [this] x [v] */
		constexpr num::Vec<Type> cross(const num::Vec<Type>& v) const {
			return Vec<Type>{
				y* v.z - z * v.y,
					z* v.x - x * v.z,
					x* v.y - y * v.x
			};
		}

		/* compute the x component of the cross product [this] x [v] */
		constexpr Type crossX(const num::Vec<Type>& v) const {
			return y * v.z - z * v.y;
		}

		/* compute the y component of the cross product [this] x [v] */
		constexpr Type crossY(const num::Vec<Type>& v) const {
			return z * v.x - x * v.z;
		}

		/* compute the z component of the cross product [this] x [v] */
		constexpr Type crossZ(const num::Vec<Type>& v) const {
			return x * v.y - y * v.x;
		}

		/* compute the vector [this] normalized with the length 1 */
		constexpr num::Vec<Type> norm() const {
			return *this / len();
		}

		/* compute the vector [this] projected onto the Y-Z plane */
		constexpr num::Vec<Type> planeX(Type xPlane = 0) const {
			return num::Vec<Type>{ xPlane, y, z };
		}

		/* compute the vector [this] projected onto the X-Z plane */
		constexpr num::Vec<Type> planeY(Type yPlane = 0) const {
			return num::Vec<Type>{ x, yPlane, z };
		}

		/* compute the vector [this] projected onto the X-Y plane */
		constexpr num::Vec<Type> planeZ(Type zPlane = 0) const {
			return num::Vec<Type>{ x, y, zPlane };
		}

		/* compute the index of the largest or smallest component of the vector */
		constexpr size_t comp(bool largest) const {
			size_t index = 0;

			/* iterate through the components and check if one is larger */
			for (size_t i = 1; i < 3; i++) {
				if (largest ? num::Abs(c[index]) < num::Abs(c[i]) : num::Abs(c[index]) > num::Abs(c[i]))
					index = i;
			}
			return index;
		}

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the x axis when it points towards the observer */
		constexpr num::Vec<Type> rotateX(Type a) const {
			a = num::ToRadian(a);
			const Type sa = std::sin(a);
			const Type ca = std::cos(a);
			return num::Vec<Type>{
				x,
					y* ca - z * sa,
					y* sa + z * ca
			};
		}

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the y axis when it points towards the observer */
		constexpr num::Vec<Type> rotateY(Type a) const {
			a = num::ToRadian(a);
			const Type sa = std::sin(a);
			const Type ca = std::cos(a);
			return num::Vec<Type>{
				x* ca + z * sa,
					y,
					z* ca - x * sa
			};
		}

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the z axis when it points towards the observer */
		constexpr num::Vec<Type> rotateZ(Type a) const {
			a = num::ToRadian(a);
			const Type sa = std::sin(a);
			const Type ca = std::cos(a);
			return num::Vec<Type>{
				x* ca - y * sa,
					x* sa + y * ca,
					z
			};
		}

		/* compute the angle to rotate [this] by on the x axis to match the vector [v] when projected onto the y/z plane [-180; 180] */
		constexpr Type angleX(const num::Vec<Type>& v) const {
			num::Vec<Type> flat = planeX();
			num::Vec<Type> target = v.planeX();

			/* compute the angle between the two vectors and correct its sign */
			Type angle = flat.angle(target);
			return (flat.crossX(target) < 0) ? -angle : angle;
		}

		/* compute the angle to rotate [this] by on the y axis to match the vector [v] when projected onto the x/z plane [-180; 180] */
		constexpr Type angleY(const num::Vec<Type>& v) const {
			num::Vec<Type> flat = planeY();
			num::Vec<Type> target = v.planeY();

			/* compute the angle between the two vectors and correct its sign */
			Type angle = flat.angle(target);
			return (flat.crossY(target) < 0) ? -angle : angle;
		}

		/* compute the angle to rotate [this] by on the z axis to match the vector [v] when projected onto the x/y plane [-180; 180] */
		constexpr Type angleZ(const num::Vec<Type>& v) const {
			num::Vec<Type> flat = planeZ();
			num::Vec<Type> target = v.planeZ();

			/* compute the angle between the two vectors and correct its sign */
			Type angle = flat.angle(target);
			return (flat.crossZ(target) < 0) ? -angle : angle;
		}

		/* construct the line [this:(p-this)] */
		constexpr num::Line<Type> line(const num::Vec<Type>& p) const;

		/* construct the plane [this:(p0-this):(p1-this)] */
		constexpr num::Plane<Type> plane(const num::Vec<Type>& p0, const num::Vec<Type>& p1) const;

		/* construct a vector interpolated between [this] and the vector [v] at [t] */
		constexpr num::Vec<Type> interpolate(const num::Vec<Type>& p, Type t) const {
			return num::Vec<Type>{
				x + (p.x - x) * t,
					y + (p.y - y) * t,
					z + (p.z - z) * t
			};
		}

		/* compute the factor which multiplied with [this] will result in a vector which is parallel to [this] but has length [l] */
		constexpr Type rescalef(Type l) const {
			return std::sqrt((l * l) / lenSquared());
		}

		/* construct a vector which is parallel to [this] but has length [l] */
		constexpr num::Vec<Type> rescale(Type l) const {
			return *this * rescalef(l);
		}

		/* compute the factor with which to scale [this] to be equal to vector [v] (result only valid if the vectors are parallel) */
		constexpr Type delta(const num::Vec<Type>& v) const {
			/* extract the largest component and use it to compute the scaling factor */
			size_t index = comp(true);
			return v.c[index] / c[index];
		}

		/* construct a vector parallel to [this] but scaled by [f] */
		constexpr num::Vec<Type> scale(Type f) const {
			return num::Vec<Type>{ x* f, y* f, z* f };
		}

		/* check if [this] and [v] describe the same vector but scaled by any factor */
		constexpr bool parallel(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			/* extract the largest components and check if the vectors are considered zero */
			const size_t largest[2] = { comp(true), v.comp(true) };
			if (num::Abs(c[largest[0]]) <= precision)
				return (num::Abs(v.c[largest[1]]) <= precision);
			else if (num::Abs(v.c[largest[1]]) <= precision)
				return false;

			/* compute the factor with which to scale the other vector */
			const Type f = c[largest[0]] / v.c[largest[1]];

			/* check if the vectors are equal when scaled */
			return match(v * f, precision);
		}

		/* check if [this] and [v] describe the same vector but scaled by a positive factor */
		constexpr bool sign(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			/* extract the largest components and check if the vectors are considered zero */
			const size_t largest[2] = { comp(true), v.comp(true) };
			if (num::Abs(c[largest[0]]) <= precision)
				return (num::Abs(v.c[largest[1]]) <= precision);
			else if (num::Abs(v.c[largest[1]]) <= precision)
				return false;

			/* compute the factor with which to scale the other vector and ensure that the factor is positive */
			const Type f = c[largest[0]] / v.c[largest[1]];
			if (f < 0)
				return false;

			/* check if the vectors are equal when scaled */
			return match(v * f, precision);
		}

		/* check if [this] and [v] are identical (all components are weighted the same) */
		constexpr bool identical(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			/* dont subtract and then compare with zero as small errors will have a much larger effect on
			*	the result due to the canceling effects of subtraction on the information */
			return num::Cmp(x, v.x, precision) && num::Cmp(y, v.y, precision) && num::Cmp(z, v.z, precision);
		}

		/* check if the x-component of this vector is zero */
		constexpr bool zeroX(Type precision = num::Const<Type>::Precision) const {
			return num::Zero(x, precision);
		}

		/* check if the y-component of this vector is zero */
		constexpr bool zeroY(Type precision = num::Const<Type>::Precision) const {
			return num::Zero(y, precision);
		}

		/* check if the z-component of this vector is zero */
		constexpr bool zeroZ(Type precision = num::Const<Type>::Precision) const {
			return num::Zero(z, precision);
		}

		/* check if the vector is zero */
		constexpr bool zero(Type precision = num::Const<Type>::Precision) const {
			return num::Zero(lenSquared(), precision);
		}

		/* check if [this] and [v] describe the same vector (point into the same direction with the same length in regard to the precision and the magnitude of the separate components) */
		constexpr bool match(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			return num::Cmp(dot(v), lenSquared(), precision);
		}

		/* check if the x-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleX(Type precision = num::Const<Type>::Precision) const {
			return num::Cmp(lenSquared(), planeX().lenSquared(), precision);
		}

		/* check if the y-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleY(Type precision = num::Const<Type>::Precision) const {
			return num::Cmp(lenSquared(), planeY().lenSquared(), precision);
		}

		/* check if the z-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleZ(Type precision = num::Const<Type>::Precision) const {
			return num::Cmp(lenSquared(), planeZ().lenSquared(), precision);
		}

		/* check if [this] and [v] are perpendicular to each other */
		constexpr bool isPerpendicular(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			return num::Zero(dot(v), precision);
		}

		/* check if [this] and [v] form an acute angle to each other (includes perpendicular) */
		constexpr bool isAcuteAngle(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			return dot(v) >= -precision;
		}

		/* check if [this] and [v] form an obtuse angle to each other (includes perpendicular) */
		constexpr bool isObtuseAngle(const num::Vec<Type>& v, Type precision = num::Const<Type>::Precision) const {
			return dot(v) <= precision;
		}

		/* compute the factor which multiplied with [this] will result in a projection of [v] onto [this] and thereby parallel to [this] */
		constexpr Type projectf(const num::Vec<Type>& v) const {
			return dot(v) / lenSquared();
		}

		/* construct a vector which is a projection of [v] onto [this] and thereby parallel to [this] */
		constexpr num::Vec<Type> project(const num::Vec<Type>& v) const {
			return *this * projectf(v);
		}

		/* construct a vector which is perpendicular to [this] and contains only the components of [v] perpendicular to [this] (lies in the plane defined by [this] and [v]) */
		constexpr num::Vec<Type> perpendicular(const num::Vec<Type>& v) const {
			return v - project(v);
		}

		/* compute the factor which multiplied with [this] will result in a vector parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr Type reachf(const num::Vec<Type>& v) const {
			return v.lenSquared() / dot(v);
		}

		/* construct a vector which is parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr num::Vec<Type> reach(const num::Vec<Type>& v) const {
			return *this * reachf(v);
		}

		/* construct a vector which is perpendicular to [this] but when added to [this] will be parallel to [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr num::Vec<Type> passing(const num::Vec<Type>& v) const {
			return v.reach(*this) - *this;
		}

		/* compute the factor which multiplied with [this] will at least pass the vector [v] if it has not already been passed */
		constexpr Type passPointf(const num::Vec<Type>& v) const {
			/* check if the vectors are perpendicular in which case the point will always be considered passed */
			if (num::Zero(dot(v)))
				return 1;

			/* compute the factor required to let this vector reach v and return it if its greater than 1 */
			return std::max(1, reachf(v));
		}

		/* construct the vector which parallel to [this] and will at least pass the vector [v] if it has not already been passed */
		constexpr num::Vec<Type> passPoint(const num::Vec<Type>& v) const {
			return *this * passPointf(v);
		}
	};

	template <std::floating_point Type>
	constexpr num::Vec<Type> operator*(Type s, const num::Vec<Type>& v) {
		return v * s;
	}
}
