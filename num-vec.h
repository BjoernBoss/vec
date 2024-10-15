#pragma once

#include "num-common.h"

namespace num {
	struct Vec {
	public:
		/* defines the component layout in memory */
		enum Component : uint8_t {
			ComponentX = 0,
			ComponentY = 1,
			ComponentZ = 2,
		};

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
		Vec();
		Vec(float f);
		Vec(float x, float y, float z);
		Vec operator+(const Vec& v) const;
		Vec operator-(const Vec& v) const;
		Vec operator-() const;
		Vec operator*(float s) const;
		Vec operator/(float s) const;
		Vec& operator+=(const Vec& v);
		Vec& operator-=(const Vec& v);
		Vec& operator*=(float s);
		Vec& operator/=(float s);
		bool operator==(const Vec& v) const;
		bool operator!=(const Vec& v) const;

	public:
		/* create a vector on the x axis with the length [l] */
		static Vec AxisX(float l = 1.0f);

		/* create a vector on the y axis with the length [l] */
		static Vec AxisY(float l = 1.0f);

		/* create a vector on the z axis with the length [l] */
		static Vec AxisZ(float l = 1.0f);

	public:
		/* compute the dot product [this] * [v] */
		float dot(const Vec& v) const;

		/* compute the angle between the vector [this] and [v] [0; 180] */
		float angle(const Vec& v) const;

		/* compute the length of the vector */
		float len() const;

		/* compute the squared length of the vector */
		float lenSquared() const;

		/* compute the the cross product [this] x [v] */
		Vec cross(const Vec& v) const;

		/* compute the x component of the cross product [this] x [v] */
		float crossX(const Vec& v) const;

		/* compute the y component of the cross product [this] x [v] */
		float crossY(const Vec& v) const;

		/* compute the z component of the cross product [this] x [v] */
		float crossZ(const Vec& v) const;

		/* compute the vector [this] normalized with the length 1 */
		Vec norm() const;

		/* compute the vector [this] projected onto the Y-Z plane */
		Vec planeX(float xPlane = 0.0f) const;

		/* compute the vector [this] projected onto the X-Z plane */
		Vec planeY(float yPlane = 0.0f) const;

		/* compute the vector [this] projected onto the X-Y plane */
		Vec planeZ(float zPlane = 0.0f) const;

		/* compute the index of the largest or smallest component of the vector */
		size_t comp(bool largest) const;

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the x axis when it points towards the observer */
		Vec rotateX(float a) const;

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the y axis when it points towards the observer */
		Vec rotateY(float a) const;

		/* compute the vector when rotating [this] by [a] degrees counterclockwise along the z axis when it points towards the observer */
		Vec rotateZ(float a) const;

		/* compute the angle to rotate [this] by on the x axis to match the vector [v] when projected onto the y/z plane [-180; 180] */
		float angleX(const Vec& v) const;

		/* compute the angle to rotate [this] by on the y axis to match the vector [v] when projected onto the x/z plane [-180; 180] */
		float angleY(const Vec& v) const;

		/* compute the angle to rotate [this] by on the z axis to match the vector [v] when projected onto the x/y plane [-180; 180] */
		float angleZ(const Vec& v) const;

		/* construct the line [this:(p-this)] */
		Line line(const Vec& p) const;

		/* construct the plane [this:(p0-this):(p1-this)] */
		Plane plane(const Vec& p0, const Vec& p1) const;

		/* construct a vector interpolated between [this] and the vector [v] at [t] */
		Vec interpolate(const Vec& p, float t) const;

		/* construct a vector which is parallel to [this] but has length [l] */
		Vec rescale(float l) const;

		/* compute the factor with which to scale [this] to be equal to vector [v] (result only valid if the vectors are parallel) */
		float delta(const Vec& v) const;

		/* construct a vector parallel to [this] but scaled by [f] */
		Vec scale(float f) const;

		/* check if [this] and [v] describe the same vector but scaled by any factor */
		bool parallel(const Vec& v, float precision = num::Precision) const;

		/* check if [this] and [v] describe the same vector but scaled by a positive factor */
		bool sign(const Vec& v, float precision = num::Precision) const;

		/* check if [this] and [v] are identical (all components are weighted the same) */
		bool equal(const Vec& v, float precision = num::Precision) const;

		/* check if the x-component of this vector is zero */
		bool zeroX(float precision = num::Precision) const;

		/* check if the y-component of this vector is zero */
		bool zeroY(float precision = num::Precision) const;

		/* check if the z-component of this vector is zero */
		bool zeroZ(float precision = num::Precision) const;

		/* check if the vector is zero */
		bool zero(float precision = num::Precision) const;

		/* check if [this] and [v] describe the same vector (point into the same direction with the same length in regard to the precision and the magnitude of the separate components) */
		bool match(const Vec& v, float precision = num::Precision) const;

		/* check if the x-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		bool negligibleX(float precision = num::Precision) const;

		/* check if the y-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		bool negligibleY(float precision = num::Precision) const;

		/* check if the z-component of the vector is negligible relative to the other components (in regard to the precision and the magnitude of the separate components) */
		bool negligibleZ(float precision = num::Precision) const;

		/* check if [this] and [v] are perpendicular to each other */
		bool isPerpendicular(const Vec& v, float precision = num::Precision) const;

		/* check if [this] and [v] form an acute angle to each other (includes perpendicular) */
		bool isAcuteAngle(const Vec& v, float precision = num::Precision) const;

		/* check if [this] and [v] form an obtuse angle to each other (includes perpendicular) */
		bool isObtuseAngle(const Vec& v, float precision = num::Precision) const;

		/* compute the factor which multiplied with [this] will result in a projection of [v] onto [this] and thereby parallel to [this] */
		float projectf(const Vec& v) const;

		/* compute a vector which is a projection of [v] onto [this] and thereby parallel to [this] */
		Vec project(const Vec& v) const;

		/* compute a vector which is perpendicular to [this] and contains only the components of [v] perpendicular to [this] (lies in the plane defined by [this] and [v]) */
		Vec perpendicular(const Vec& v) const;

		/* compute the factor which multiplied with [this] will result in a vector parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		float reachf(const Vec& v) const;

		/* compute a vector which is parallel to [this] which reaches the point [v] i.e. it is perpendicular to [v] when subtracted from [v] (invalid result if [this] and [v] are perpendicular) */
		Vec reach(const Vec& v) const;

		/* compute a vector which is perpendicular to [this] but when added to [this] will be parallel to [v] (invalid result if [this] and [v] are perpendicular) */
		Vec passing(const Vec& v) const;

		/* compute the factor which multiplied with [this] will at least pass the vector [v] if it has not already been passed */
		float passPointf(const Vec& v) const;

		/* compute the vector which parallel to [this] and will at least pass the vector [v] if it has not already been passed */
		Vec passPoint(const Vec& v) const;
	};

	Vec operator*(float s, const Vec& v);

	std::ostream& operator<<(std::ostream& out, const Vec& v);
	std::wostream& operator<<(std::wostream& out, const Vec& v);
	std::istream& operator>>(std::istream& in, Vec& v);
	std::wistream& operator>>(std::wistream& in, Vec& v);
}
