#pragma once

#include "num-common.h"

namespace num {
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

		/* compute the area of the triangle created by the plane [this] when projected onto the Y-Z plane */
		float areaX() const;

		/* compute the area of the triangle created by the plane [this] when projected onto the X-Z plane */
		float areaY() const;

		/* compute the area of the triangle created by the plane [this] when projected onto the X-Y plane */
		float areaZ() const;

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

	std::ostream& operator<<(std::ostream& out, const Plane& p);
	std::wostream& operator<<(std::wostream& out, const Plane& p);

	std::istream& operator>>(std::istream& in, Plane& p);
	std::wistream& operator>>(std::wistream& in, Plane& p);
}
