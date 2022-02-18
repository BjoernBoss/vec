#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <utility>

#include "float-base.h"

#include <istream>
#include <ostream>

namespace num {
	/*
	*	- right-handed system
	*	- counterclockwise rotations when the corresponding axis points towards the observer
	*	- angles are calculated in degrees
	*/
	struct Line;
	struct Plane;

	/* define the linear helper object */
	struct Linear {
		float s;
		float t;

	public:
		Linear();
		Linear(float s, float t);
	};

	/* define the vector object */
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

		/* check if [this] and [v] describe the same vector (point into the same direction with the same length in regard to the precision and the magnitude of the separate components) */
		bool match(const Vec& v, float precision = num::Precision) const;

		/* check if [this] and [v] are identical (all components are weighted the same) */
		bool equal(const Vec& v, float precision = num::Precision) const;

		/* check if the x-component of this vector is negligible */
		bool zeroX(float precision = num::Precision) const;

		/* check if the y-component of this vector is negligible */
		bool zeroY(float precision = num::Precision) const;

		/* check if the z-component of this vector is negligible */
		bool zeroZ(float precision = num::Precision) const;

		/* check if the vector is negligible */
		bool zero(float precision = num::Precision) const;

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

	/* define the line object */
	struct Line {
		Vec o;
		Vec d;

	public:
		Line();
		Line(const Vec& d);
		Line(const Vec& o, const Vec& d);

	private:
		/* compute the linear combination of the line [this] and line [l] to intersect based on the two other axes than the index axis */
		Linear fLinComb(const Line& l, size_t index, bool& parallel, float precision) const;

	public:
		/* create a line along the x axis with the length [l] and its origin the the origin */
		static Line AxisX(float l = 1.0f);

		/* create a line along the y axis with the length [l] and its origin the the origin */
		static Line AxisY(float l = 1.0f);

		/* create a line along the z axis with the length [l] and its origin the the origin */
		static Line AxisZ(float l = 1.0f);

	public:
		/* compute the line [this] projected onto the Y-Z plane */
		Line planeX(float xPlane = 0.0f) const;

		/* compute the line [this] projected onto the X-Z plane */
		Line planeY(float yPlane = 0.0f) const;

		/* compute the line [this] projected onto the X-Y plane */
		Line planeZ(float zPlane = 0.0f) const;

		/* compute a point on line [this] */
		Vec point(float t) const;

		/* compute a normalized origin that is close to zero normalize the direction */
		Line norm() const;

		/* check if [p] lies on line [this] */
		bool touch(const Vec& p, float precision = num::Precision) const;

		/* returns a position along the [this] line where the point [p] lies (result only valid if it lies on the line) */
		float find(const Vec& p) const;

		/* check if the line [l] and line [this] describe the same line */
		bool equal(const Line& l, float precision = num::Precision) const;

		/* check if the line [l] and line [this] describe the identically same line */
		bool identical(const Line& l, float precision = num::Precision) const;

		/* compute the factor for which line [this] reaches the point closest to p (automatically perpendicular) */
		float closestf(const Vec& p) const;

		/* compute the shortest vector which connects [p] to a point on the line (automatically perpendicular) */
		Vec closest(const Vec& p) const;

		/* compute the shortest line which intersects line [this] and the line [l] */
		Linear closestf(const Line& l) const;

		/* compute the shortest line which intersects line [this] and the line [l] */
		Line closest(const Line& l) const;

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the Y-Z plane at [xPlane] (invalid if parallel: returns [this] origin) */
		float intersectPlaneXf(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and the Y-Z plane at [xPlane] (invalid if parallel: returns null vector) */
		Vec intersectPlaneX(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the X-Z plane at [yPlane] (invalid if parallel: returns [this] origin) */
		float intersectPlaneYf(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and the X-Z plane at [yPlane] (invalid if parallel: returns null vector) */
		Vec intersectPlaneY(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] with to reach the intersection point of this line and the X-Y plane at [zPlane] (invalid if parallel: returns [this] origin) */
		float intersectPlaneZf(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and the X-Y plane at [zPlane] (invalid if parallel: returns null vector) */
		Vec intersectPlaneZ(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the X-Y plane (invalid if no intersection point: returns 0, 0) */
		Linear intersectXf(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and line [l] when viewed in the X-Y plane (invalid if no intersection point: returns [this] origin) */
		Vec intersectX(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the X-Z plane (invalid if no intersection point: returns 0, 0) */
		Linear intersectYf(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and line [l] when viewed in the X-Z plane (invalid if no intersection point: returns [this] origin) */
		Vec intersectY(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] and line [l] with to intersect the lines when viewed in the Y-Z plane (invalid if no intersection point: returns 0, 0) */
		Linear intersectZf(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and line [l] when viewed in the Y-Z plane (invalid if no intersection point: returns [this] origin) */
		Vec intersectZ(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the factor to scale line [this] and line [l] with to intersect the lines (invalid if no intersection point: returns 0, 0) */
		Linear intersectf(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of line [this] and the line [l] (invalid if no intersection point: returns [this] origin) */
		Vec intersect(const Line& l, bool* invalid = 0, float precision = num::Precision) const;
	};

	/* define the plane object */
	struct Plane {
		Vec o;
		Vec a;
		Vec b;

	public:
		Plane();
		Plane(const Vec& a, const Vec& b);
		Plane(const Vec& o, const Vec& a, const Vec& b);

	private:
		/* compute the linear combination of the two extent vectors to the point in a plane based on the index axis */
		Linear fLinComb(const Vec& p, size_t index) const;

	public:
		/* create a plane parallel to the Y-Z plane at distance [d] to the origin */
		static Plane AxisX(float d = 0.0f);

		/* create a plane parallel to the X-Z plane at distance [d] to the origin */
		static Plane AxisY(float d = 1.0f);

		/* create a plane parallel to the X-Y plane at distance [d] to the origin */
		static Plane AxisZ(float d = 1.0f);

	public:
		/* compute the plane [this] projected onto the Y-Z plane */
		Plane planeX(float xPlane = 0.0f) const;

		/* compute the plane [this] projected onto the X-Z plane */
		Plane planeY(float yPlane = 0.0f) const;

		/* compute the plane [this] projected onto the X-Y plane */
		Plane planeZ(float zPlane = 0.0f) const;

		/* compute the normal vector of the plane */
		Vec normal() const;

		/* compute the area of the triangle created by the plane [this] */
		float area() const;

		/* compute the center point of the triangle produced by the two extent vectors and the origin */
		Vec center() const;

		/* compute a point on plane [this] */
		Vec point(float s, float t) const;

		/* compute a point on plane [this] */
		Vec point(const Linear& lin) const;

		/* compute a normalized origin that is close to zero and normalize the directions as well as orient them to be perpendicular */
		Plane norm() const;

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the Y-Z plane */
		Vec projectX(const Vec& p) const;

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the X-Z plane */
		Vec projectY(const Vec& p) const;

		/* compute the vector if [p] is being projected onto the plane viewed orthogonally from the X-Y plane */
		Vec projectZ(const Vec& p) const;

		/* compute the vector if [v] is being projected onto the plane */
		Vec project(const Vec& v) const;

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the Y-Z plane */
		bool inTriangleX(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the X-Z plane */
		bool inTriangleY(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the triangle of a and b when projected orthogonally onto the X-Y plane */
		bool inTriangleZ(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the triangle of a and b */
		bool inTriangle(const Vec& p, bool* touching = 0, float precision = num::Precision) const;

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the Y-Z plane */
		bool inConeX(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Z plane */
		bool inConeY(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the cone of a and b when projected orthogonally onto the X-Y plane */
		bool inConeZ(const Vec& p, float precision = num::Precision) const;

		/* check if [p] lies within the cone of a and b */
		bool inCone(const Vec& p, bool* touching = 0, float precision = num::Precision) const;

		/* check if [p] lies on the plane */
		bool touch(const Vec& p, float precision = num::Precision) const;

		/* check if the plane [p] and plane [this] describe the same plane */
		bool equal(const Plane& p, float precision = num::Precision) const;

		/* check if the plane [p] and plane [this] describe the identically same plane */
		bool identical(const Plane& p, float precision = num::Precision) const;

		/* compute the shortest vector which connects [p] to a point on the plane (automatically perpendicular) */
		Vec closest(const Vec& p) const;

		/* compute the vector of steepest ascent in the plane for the X axis */
		Vec steepestX() const;

		/* compute the vector of steepest ascent in the plane for the Y axis */
		Vec steepestY() const;

		/* compute the vector of steepest ascent in the plane for the Z axis */
		Vec steepestZ() const;

		/* compute the intersecting line between the Y-Z plane at [xPlane] and plane [this] (invalid if parallel: returns null line) */
		Line intersectPlaneX(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersecting line between the X-Z plane at [yPlane] and plane [this] (invalid if parallel: returns null line) */
		Line intersectPlaneY(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersecting line between the X-Y plane at [zPlane] and plane [this] (invalid if parallel: returns null line) */
		Line intersectPlaneZ(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection line of the plane [this] and the plane [p] (invalid if parallel: returns null line) */
		Line intersect(const Plane& p, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection factors of the plane [this] and the line [l] (invalid if parallel: returns 0, 0) */
		Linear intersectf(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the intersection point of the plane [this] and the line [l] (invalid if parallel: returns null vector) */
		Vec intersect(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the Y-Z plane */
		Linear linearX(const Vec& p) const;

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the X-Z plane */
		Linear linearY(const Vec& p) const;

		/* compute the linear combination to reach the point [p] on plane [this] when projected orthogonally onto the X-Y plane */
		Linear linearZ(const Vec& p) const;

		/* compute the linear combination to reach the point [p] when the point lies on the plane */
		Linear linear(const Vec& p, bool* touching = 0, float precision = num::Precision) const;
	};

	/* define the vector multiplication from the right */
	Vec operator*(float s, const Vec& v);

	/* define the out streaming operator */
	std::ostream& operator<<(std::ostream& out, const Linear& l);
	std::wostream& operator<<(std::wostream& out, const Linear& l);
	std::ostream& operator<<(std::ostream& out, const Vec& v);
	std::wostream& operator<<(std::wostream& out, const Vec& v);
	std::ostream& operator<<(std::ostream& out, const Line& l);
	std::wostream& operator<<(std::wostream& out, const Line& l);
	std::ostream& operator<<(std::ostream& out, const Plane& p);
	std::wostream& operator<<(std::wostream& out, const Plane& p);

	/* define the in streaming operator */
	std::istream& operator>>(std::istream& in, Linear& l);
	std::wistream& operator>>(std::wistream& in, Linear& l);
	std::istream& operator>>(std::istream& in, Vec& v);
	std::wistream& operator>>(std::wistream& in, Vec& v);
	std::istream& operator>>(std::istream& in, Line& l);
	std::wistream& operator>>(std::wistream& in, Line& l);
	std::istream& operator>>(std::istream& in, Plane& p);
	std::wistream& operator>>(std::wistream& in, Plane& p);
}
