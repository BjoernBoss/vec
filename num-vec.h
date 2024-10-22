#pragma once

#include "num-common.h"

namespace num {
	/* defines the component layout in memory */
	enum Component : uint8_t {
		ComponentX = 0,
		ComponentY = 1,
		ComponentZ = 2,
	};

	struct Vec {
	public:
		union {
			struct {
				float x;
				float y;
				float z;
			};
			float c[3];
		};

	public:
		constexpr Vec() : x{ 0 }, y{ 0 }, z{ 0 } {}
		constexpr Vec(float f) : x{ f }, y{ f }, z{ f } {}
		constexpr Vec(float x, float y, float z) : x{ x }, y{ y }, z{ z } {}

	public:
		constexpr num::Vec operator+(const num::Vec& v) const {
			return num::Vec{ x + v.x, y + v.y, z + v.z };
		}
		constexpr num::Vec operator-(const num::Vec& v) const {
			return num::Vec{ x - v.x, y - v.y, z - v.z };
		}
		constexpr num::Vec operator-() const {
			return num::Vec{ -x, -y, -z };
		}
		constexpr num::Vec operator*(float s) const {
			return num::Vec{ x * s, y * s, z * s };
		}
		constexpr num::Vec operator/(float s) const {
			return num::Vec{ x / s, y / s, z / s };
		}
		constexpr num::Vec& operator+=(const num::Vec& v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		constexpr num::Vec& operator-=(const num::Vec& v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		constexpr num::Vec& operator*=(float s) {
			x *= s;
			y *= s;
			z *= s;
			return *this;
		}
		constexpr num::Vec& operator/=(float s) {
			x /= s;
			y /= s;
			z /= s;
			return *this;
		}
		constexpr bool operator==(const num::Vec& v) const {
			return equal(v);
		}
		constexpr bool operator!=(const num::Vec& v) const {
			return !(*this == v);
		}

	public:
		/* create a vector on the x axis with the length [l] */
		static constexpr num::Vec AxisX(float l = 1.0f) {
			return num::Vec{ l, 0.0f, 0.0f };
		}

		/* create a vector on the y axis with the length [l] */
		static constexpr num::Vec AxisY(float l = 1.0f) {
			return num::Vec{ 0.0f, l, 0.0f };
		}

		/* create a vector on the z axis with the length [l] */
		static constexpr num::Vec AxisZ(float l = 1.0f) {
			return num::Vec{ 0.0f, 0.0f, l };
		}

	public:
		/* compute the dot product [this] * [v] */
		constexpr float dot(const num::Vec& v) const {
			return v.x * x + v.y * y + v.z * z;
		}

		/* compute the angle between the vector [this] and [v] [0; 180] */
		float angle(const Vec& v) const {
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

		/* compute the squared length of the vector */
		constexpr float lenSquared() const {
			return dot(*this);
		}

		/* compute the length of the vector */
		float len() const {
			return std::sqrt(lenSquared());
		}

		/* compute the the cross product [this] x [v] */
		constexpr num::Vec cross(const num::Vec& v) const {
			return Vec{
				y * v.z - z * v.y,
				z * v.x - x * v.z,
				x * v.y - y * v.x
			};
		}

		/* compute the x component of the cross product [this] x [v] */
		constexpr float crossX(const num::Vec& v) const {
			return y * v.z - z * v.y;
		}

		/* compute the y component of the cross product [this] x [v] */
		constexpr float crossY(const num::Vec& v) const {
			return z * v.x - x * v.z;
		}

		/* compute the z component of the cross product [this] x [v] */
		constexpr float crossZ(const num::Vec& v) const {
			return x * v.y - y * v.x;
		}

		/* compute the vector [this] normalized with the length 1 */
		constexpr num::Vec norm() const {
			return *this / len();
		}

		/* compute the vector [this] projected onto the Y-Z plane */
		constexpr num::Vec planeX(float xPlane = 0.0f) const {
			return num::Vec{ xPlane, y, z };
		}

		/* compute the vector [this] projected onto the X-Z plane */
		constexpr num::Vec planeY(float yPlane = 0.0f) const {
			return num::Vec{ x, yPlane, z };
		}

		/* compute the vector [this] projected onto the X-Y plane */
		constexpr num::Vec planeZ(float zPlane = 0.0f) const {
			return num::Vec{ x, y, zPlane };
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
		num::Vec rotateX(float a) const {
			a = num::ToRadian(a);
			const float sa = std::sin(a);
			const float ca = std::cos(a);
			return num::Vec{
				x,
				y * ca - z * sa,
				y * sa + z * ca
			};
		}

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the y axis when it points towards the observer */
		num::Vec rotateY(float a) const {
			a = num::ToRadian(a);
			const float sa = std::sin(a);
			const float ca = std::cos(a);
			return num::Vec{
				x * ca + z * sa,
				y,
				z * ca - x * sa
			};
		}

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the z axis when it points towards the observer */
		num::Vec rotateZ(float a) const {
			a = num::ToRadian(a);
			const float sa = std::sin(a);
			const float ca = std::cos(a);
			return num::Vec{
				x * ca - y * sa,
				x * sa + y * ca,
				z
			};
		}

		/* compute the angle to rotate [this] by on the x axis to match the vector [v] when projected onto the y/z plane [-180; 180] */
		float angleX(const num::Vec& v) const {
			num::Vec flat = planeX();
			num::Vec target = v.planeX();

			/* compute the angle between the two vectors and correct its sign */
			float angle = flat.angle(target);
			return (flat.crossX(target) < 0.0f) ? -angle : angle;
		}

		/* compute the angle to rotate [this] by on the y axis to match the vector [v] when projected onto the x/z plane [-180; 180] */
		float angleY(const num::Vec& v) const {
			num::Vec flat = planeY();
			num::Vec target = v.planeY();

			/* compute the angle between the two vectors and correct its sign */
			float angle = flat.angle(target);
			return (flat.crossY(target) < 0.0f) ? -angle : angle;
		}

		/* compute the angle to rotate [this] by on the z axis to match the vector [v] when projected onto the x/y plane [-180; 180] */
		float angleZ(const num::Vec& v) const {
			num::Vec flat = planeZ();
			num::Vec target = v.planeZ();

			/* compute the angle between the two vectors and correct its sign */
			float angle = flat.angle(target);
			return (flat.crossZ(target) < 0.0f) ? -angle : angle;
		}

		/* construct the line [this:(p-this)] */
		constexpr num::Line line(const num::Vec& p) const;

		/* construct the plane [this:(p0-this):(p1-this)] */
		constexpr num::Plane plane(const num::Vec& p0, const num::Vec& p1) const;

		/* construct a vector interpolated between [this] and the vector [v] at [t] */
		constexpr num::Vec interpolate(const num::Vec& p, float t) const {
			return num::Vec{
				x + (p.x - x) * t,
				y + (p.y - y) * t,
				z + (p.z - z) * t
			};
		}

		/* compute the factor which multiplied with [this] will result in a vector which is parallel to [this] but has length [l] */
		float rescalef(float l) const {
			return std::sqrt((l * l) / lenSquared());
		}

		/* construct a vector which is parallel to [this] but has length [l] */
		constexpr num::Vec rescale(float l) const {
			return *this * rescalef(l);
		}

		/* compute the factor with which to scale [this] to be equal to vector [v] (result only valid if the vectors are parallel) */
		constexpr float delta(const num::Vec& v) const {
			/* extract the largest component and use it to compute the scaling factor */
			size_t index = comp(true);
			return v.c[index] / c[index];
		}

		/* construct a vector parallel to [this] but scaled by [f] */
		constexpr num::Vec scale(float f) const {
			return num::Vec{ x * f, y * f, z * f };
		}

		/* check if [this] and [v] describe the same vector but scaled by any factor */
		constexpr bool parallel(const num::Vec& v, float precision = num::Precision) const {
			/* extract the largest components and check if the vectors are considered zero */
			const size_t largest[2] = { comp(true), v.comp(true) };
			if (num::Abs(c[largest[0]]) <= precision)
				return (num::Abs(v.c[largest[1]]) <= precision);
			else if (num::Abs(v.c[largest[1]]) <= precision)
				return false;

			/* compute the factor with which to scale the other vector */
			const float f = c[largest[0]] / v.c[largest[1]];

			/* check if the vectors are equal when scaled */
			return equal(v * f, precision);
		}

		/* check if [this] and [v] describe the same vector but scaled by a positive factor */
		constexpr bool sign(const num::Vec& v, float precision = num::Precision) const {
			/* extract the largest components and check if the vectors are considered zero */
			const size_t largest[2] = { comp(true), v.comp(true) };
			if (num::Abs(c[largest[0]]) <= precision)
				return (num::Abs(v.c[largest[1]]) <= precision);
			else if (num::Abs(v.c[largest[1]]) <= precision)
				return false;

			/* compute the factor with which to scale the other vector and ensure that the factor is positive */
			const float f = c[largest[0]] / v.c[largest[1]];
			if (f < 0.0f)
				return false;

			/* check if the vectors are equal when scaled */
			return equal(v * f, precision);
		}

		/* check if [this] and [v] are identical (all components are weighted the same) */
		constexpr bool equal(const num::Vec& v, float precision = num::Precision) const {
			/* dont subtract and then compare with zero as small errors will have a much larger effect on
			*	the result due to the canceling effects of subtraction on the information */
			return num::Cmp(x, v.x, precision) && num::Cmp(y, v.y, precision) && num::Cmp(z, v.z, precision);
		}

		/* check if the x-component of this vector is zero */
		constexpr bool zeroX(float precision = num::Precision) const {
			return num::Zero(x, precision);
		}

		/* check if the y-component of this vector is zero */
		constexpr bool zeroY(float precision = num::Precision) const {
			return num::Zero(y, precision);
		}

		/* check if the z-component of this vector is zero */
		constexpr bool zeroZ(float precision = num::Precision) const {
			return num::Zero(z, precision);
		}

		/* check if the vector is zero */
		constexpr bool zero(float precision = num::Precision) const {
			return num::Zero(lenSquared(), precision);
		}

		/* check if [this] and [v] describe the same vector (point into the same direction with the same length in regard to the precision and the magnitude of the separate components) */
		constexpr bool match(const num::Vec& v, float precision = num::Precision) const {
			return num::Cmp(dot(v), lenSquared(), precision);
		}

		/* check if the x-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleX(float precision = num::Precision) const {
			return num::Cmp(lenSquared(), planeX().lenSquared(), precision);
		}

		/* check if the y-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleY(float precision = num::Precision) const {
			return num::Cmp(lenSquared(), planeY().lenSquared(), precision);
		}

		/* check if the z-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		constexpr bool negligibleZ(float precision = num::Precision) const {
			return num::Cmp(lenSquared(), planeZ().lenSquared(), precision);
		}

		/* check if [this] and [v] are perpendicular to each other */
		constexpr bool isPerpendicular(const num::Vec& v, float precision = num::Precision) const {
			return num::Zero(dot(v), precision);
		}

		/* check if [this] and [v] form an acute angle to each other (includes perpendicular) */
		constexpr bool isAcuteAngle(const num::Vec& v, float precision = num::Precision) const {
			return dot(v) >= -precision;
		}

		/* check if [this] and [v] form an obtuse angle to each other (includes perpendicular) */
		constexpr bool isObtuseAngle(const num::Vec& v, float precision = num::Precision) const {
			return dot(v) <= precision;
		}

		/* compute the factor which multiplied with [this] will result in a projection of [v] onto [this] and thereby parallel to [this] */
		constexpr float projectf(const num::Vec& v) const {
			return dot(v) / lenSquared();
		}

		/* construct a vector which is a projection of [v] onto [this] and thereby parallel to [this] */
		constexpr num::Vec project(const num::Vec& v) const {
			return *this * projectf(v);
		}

		/* construct a vector which is perpendicular to [this] and contains only the components of [v] perpendicular to [this] (lies in the plane defined by [this] and [v]) */
		constexpr num::Vec perpendicular(const num::Vec& v) const {
			return v - project(v);
		}

		/* compute the factor which multiplied with [this] will result in a vector parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr float reachf(const num::Vec& v) const {
			return v.lenSquared() / dot(v);
		}

		/* construct a vector which is parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr num::Vec reach(const num::Vec& v) const {
			return *this * reachf(v);
		}

		/* construct a vector which is perpendicular to [this] but when added to [this] will be parallel to [v] (invalid result if [this] and [v] are perpendicular) */
		constexpr num::Vec passing(const num::Vec& v) const {
			return v.reach(*this) - *this;
		}

		/* compute the factor which multiplied with [this] will at least pass the vector [v] if it has not already been passed */
		constexpr float passPointf(const num::Vec& v) const {
			/* check if the vectors are perpendicular in which case the point will always be considered passed */
			if (num::Zero(dot(v)))
				return 1.0f;

			/* compute the factor required to let this vector reach v and return it if its greater than 1 */
			return std::max(1.0f, reachf(v));
		}

		/* construct the vector which parallel to [this] and will at least pass the vector [v] if it has not already been passed */
		constexpr num::Vec passPoint(const num::Vec& v) const {
			return *this * passPointf(v);
		}
	};

	constexpr num::Vec operator*(float s, const num::Vec& v) {
		return v * s;
	}
}
