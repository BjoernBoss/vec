#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <utility>
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

	static constexpr float Pi = 3.1415926536f;
	static constexpr float SqrtTwo = 1.4142135624f;
	static constexpr float Precision = 0.00001f;
	static constexpr float ZeroPrecisionFactor = 0.01f;

	/* define the float zero comparison function */
	bool Zero(float a, float p = Precision);

	/* define the float comparison function */
	bool Cmp(float a, float b, float p = Precision);

	/* define the angle conversion functions */
	static constexpr float ToRadian(float deg) {
		return (deg * Pi) / 180.0f;
	}
	static constexpr float ToDegree(float deg) {
		return (deg * 180.0f) / Pi;
	}
	float ToAngle(float x, float y);

	/* define the angle comparison functions which computes the angle to add to [base] to reach [test] in degrees */
	float AngleDiff(float base, float test);

	/* define the angle comparison functions which computes the absolute difference between [base] and [test] in degrees */
	float AngleAbs(float base, float test);

	struct Linear {
	public:
		float s = 0.0f;
		float t = 0.0f;

	public:
		Linear();
		Linear(float s, float t);
	};

	std::ostream& operator<<(std::ostream& out, const Linear& l);
	std::wostream& operator<<(std::wostream& out, const Linear& l);
	std::istream& operator>>(std::istream& in, Linear& l);
	std::wistream& operator>>(std::wistream& in, Linear& l);
}
