#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>

/* define the float number helper */
namespace num {
	static constexpr float ConstPi = 3.1415926536f;

	/* define the float comparison function */
	static constexpr float Precision = 0.00001f;
	static bool Cmp(float a, float b, float p = Precision) {
		if (std::isnan(a) || std::isnan(b))
			return false;
		if (a == 0.0f)
			return std::abs(b) <= 0.01f * p;
		if (b == 0.0f)
			return std::abs(a) <= 0.01f * p;
		const float _a = std::abs(a);
		const float _b = std::abs(b);
		return std::abs(a - b) <= std::min(_a, _b) * p;
	}

	/* define the angle conversion functions */
	static constexpr float ToRadian(float deg) {
		return (deg * ConstPi) / 180.0f;
	}
	static constexpr float ToDegree(float deg) {
		return (deg * 180.0f) / ConstPi;
	}
	static float ToAngle(float x, float y) {
		float deg = ToDegree(std::atan2(x, y));
		if (deg < 0)
			deg += 360.0f;
		return deg;
	}
}

/* define the vector object */
struct Vec {
public:
	/* defines the component layout in memory */
	enum Component : uint8_t {
		ComponentX = 0,
		ComponentY = 2,
		ComponentZ = 1
	};

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

public:
	float dot(const Vec& v) const;
	float angle(const Vec& v) const;
	float len() const;
	float lenSquared() const;
	Vec cross(const Vec& v) const;
	float crossX(const Vec& v) const;
	float crossY(const Vec& v) const;
	float crossZ(const Vec& v) const;
	Vec norm() const;
	Vec planeX() const;
	Vec planeY() const;
	Vec planeZ() const;
	Vec rotateX(float a) const;
	Vec rotateY(float a) const;
	Vec rotateZ(float a) const;
	size_t comp(bool largest) const;

public:
	struct Line;
	struct Plane;

public:
	/* construct the line [this:(p-this)] */
	Line line(const Vec& p) const;

	/* construct the plane [this:(p0-this):(p1-this)] */
	Plane plane(const Vec& p0, const Vec& p1) const;

	/* construct a vector which is parallel to [this] but has length l */
	Vec scale(float l) const;

	/* construct the vector interpolated between [this] and the vector [v] at t */
	Vec interpolate(const Vec& p, float t) const;

	/* compute the factor with which to scale [this] to be equal to vector [v] (result only valid if the vectors are parallel) */
	float scale(const Vec& v) const;

	/* check if [this] and [v] describe the same vector but scaled by any factor */
	bool parallel(const Vec& v, float precision = num::Precision) const;

	/* check if [this] and [v] describe the same vector but scaled by a positive factor */
	bool sign(const Vec& v, float precision = num::Precision) const;

	/* check if [this] and [v] describe the same vector (point into the same direction with a the same length in regard to the precision) */
	bool same(const Vec& v, float precision = num::Precision) const;

	/* check if [this] and [v] are identical (all components are weighted the same) */
	bool identical(const Vec& v, float precision = num::Precision) const;

	/* check if the x-component of this vector is negligible */
	bool zeroX(float precision = num::Precision) const;

	/* check if the y-component of this vector is negligible */
	bool zeroY(float precision = num::Precision) const;

	/* check if the z-component of this vector is negligible */
	bool zeroZ(float precision = num::Precision) const;

	/* check if the vector is negligible */
	bool zero(float precision = num::Precision) const;

	/* compute a vector which is a linear combination of [this] and [v] but is perpendicular to [this] */
	Vec perpendicular(const Vec& v) const;

	/* compute a vector which is a projection of [v] onto [this] and thereby parallel to [this] */
	Vec project(const Vec& v) const;

	/* compute the factor which multiplied with [this] will result in a projection of [v] onto [this] and thereby parallel to [this] */
	float projectFactor(const Vec& v) const;
};

/* define the line object */
struct Vec::Line {
	Vec o;
	Vec d;

public:
	struct Linear {
		float s;
		float t;

	public:
		Linear();
		Linear(float s, float t);
	};

public:
	Line();
	Line(const Vec& d);
	Line(const Vec& o, const Vec& d);

public:
	Line planeX() const;
	Line planeY() const;
	Line planeZ() const;

public:
	/* compute a point on this line */
	Vec point(float f) const;

	/* compute a normalized origin that is close to zero normalize the direction */
	Line norm() const;

	/* check if [p] lies on the line */
	bool touch(const Vec& p, float precision = num::Precision) const;

	/* check if the line [l] and this line describe the same line */
	bool same(const Line& l, float precision = num::Precision) const;

	/* check if the line [l] and this line describe the identically same line */
	bool identical(const Line& l, float precision = num::Precision) const;

	/* compute the shortest vector which connects [p] to a point on the line (automatically perpendicular) */
	Vec closest(const Vec& p) const;

	/* compute the factor for which this line reaches the point closest to p (automatically perpendicular) */
	float closestFactor(const Vec& p) const;

	/* compute the shortest line which intersects this line and the line [l] */
	Line closest(const Line& l) const;

	/* compute the shortest line which intersects this line and the line [l] */
	Linear closestFactor(const Line& l) const;

	/* compute the intersection point of this line and the Y-Z plane at [xPlane] */
	Vec intersectX(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersection point of this line and the X-Z plane at [yPlane] */
	Vec intersectY(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersection point of this line and the X-Y plane at [zPlane] */
	Vec intersectZ(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the factor to scale this line with to reach the intersection point of this line and the Y-Z plane at [xPlane] */
	float intersectXFactor(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the factor to scale this line with to reach the intersection point of this line and the X-Z plane at [yPlane] */
	float intersectYFactor(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the factor to scale this line with to reach the intersection point of this line and the X-Y plane at [zPlane] */
	float intersectZFactor(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersection point of this line and the line [l] */
	Vec intersect(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the factor to scale this lines direction and [l]'s direction with to intersect the lines */
	Linear intersectFactor(const Line& l, bool* invalid = 0, float precision = num::Precision) const;
};

/* define the plane object */
struct Vec::Plane {
	Vec o;
	Vec a;
	Vec b;

public:
	struct Linear {
		float s;
		float t;

	public:
		Linear();
		Linear(float _s, float _t);
	};

public:
	Plane();
	Plane(const Vec& a, const Vec& b);
	Plane(const Vec& o, const Vec& a, const Vec& b);

private:
	/* compute the linear combination of the two extent vectors to the point [this] in a plane based on the index */
	Linear fLinComb(const Vec& p, size_t index) const;

public:
	Plane planeX() const;
	Plane planeY() const;
	Plane planeZ() const;

public:
	/* compute the normal vector of the plane */
	Vec normal() const;

	/* compute the center point of the triangle produced by the two extent vectors and the origin */
	Vec center() const;

	/* compute a point on this plane */
	Vec point(float s, float t) const;

	/* compute a point on this plane */
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

	/* check if the plane [p] and this plane describe the same plane */
	bool same(const Plane& p, float precision = num::Precision) const;

	/* check if the plane [p] and this plane describe the identically same plane */
	bool identical(const Plane& p, float precision = num::Precision) const;

	/* compute the shortest vector which connects [p] to a point on the plane (automatically perpendicular) */
	Vec closest(const Vec& p) const;

	/* compute the vector of steepest ascent in the plane for the X axis */
	Vec steepestX() const;

	/* compute the vector of steepest ascent in the plane for the Y axis */
	Vec steepestY() const;

	/* compute the vector of steepest ascent in the plane for the Z axis */
	Vec steepestZ() const;

	/* compute the intersecting line between the Y-Z plane at [xPlane] and this plane */
	Line intersectX(float xPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersecting line between the X-Z plane at [yPlane] and this plane */
	Line intersectY(float yPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersecting line between the X-Y plane at [zPlane] and this plane */
	Line intersectZ(float zPlane, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersection line of the this and the plane [p] */
	Line intersect(const Plane& p, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the intersection point of the line [l] and this plane */
	Vec intersect(const Line& l, bool* invalid = 0, float precision = num::Precision) const;

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the Y-Z plane */
	Linear linearX(const Vec& p) const;

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the X-Z plane */
	Linear linearY(const Vec& p) const;

	/* compute the linear combination to reach the point [p] on this plane when projected orthogonally onto the X-Y plane */
	Linear linearZ(const Vec& p) const;

	/* compute the linear combination to reach the point [p] when the point lies on the plane */
	Linear linear(const Vec& p, bool* touching = 0, float precision = num::Precision) const;
};
