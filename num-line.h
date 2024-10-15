#pragma once

#include "num-common.h"

namespace num {
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

	std::ostream& operator<<(std::ostream& out, const Line& l);
	std::wostream& operator<<(std::wostream& out, const Line& l);

	std::istream& operator>>(std::istream& in, Line& l);
	std::wistream& operator>>(std::wistream& in, Line& l);
}
