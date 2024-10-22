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
	static constexpr float Precision = 0.00001f;
	static constexpr float ZeroPrecisionFactor = 0.01f;

	/* float abs-function (not using std implementation to allow for constexpr) */
	constexpr float Abs(float v) {
		return (v < 0 ? -v : v);
	}

	/* check if number can be considered zero */
	constexpr bool Zero(float a, float p = num::Precision) {
		/* dont check for nan as nan will fail this check and thereby return false by default */
		return num::Abs(a) <= num::ZeroPrecisionFactor * p;
	}

	/* compare the values for equality, given the corresponding precision */
	constexpr bool Cmp(float a, float b, float p = num::Precision) {
		if (std::isnan(a) || std::isnan(b))
			return false;
		if (a == 0.0f)
			return num::Zero(b);
		if (b == 0.0f)
			return num::Zero(a);
		const float _a = num::Abs(a);
		const float _b = num::Abs(b);
		return num::Abs(a - b) <= std::min(_a, _b) * p;
	}

	constexpr float ToRadian(float deg) {
		return (deg * num::Pi) / 180.0f;
	}

	constexpr float ToDegree(float deg) {
		return (deg * 180.0f) / num::Pi;
	}

	constexpr float ToAngle(float x, float y) {
		float deg = num::ToDegree(std::atan2(x, y));
		if (deg < 0)
			deg += 360.0f;
		return deg;
	}

	/* compute the angle to add to [base] to reach [test] in degrees */
	constexpr float AngleDiff(float base, float test) {
		float diff = test - base;
		if (diff <= -180.0f)
			diff += 360.0f;
		else if (diff > 180.0f)
			diff -= 360.0f;
		return diff;
	}

	/* compute the absolute difference between [base] and [test] in degrees */
	constexpr float AngleAbs(float base, float test) {
		float diff = num::Abs(test - base);
		if (diff > 180.0f)
			diff = 360.0f - diff;
		return diff;
	}

	struct Linear {
	public:
		float s = 0.0f;
		float t = 0.0f;

	public:
		constexpr Linear() : s{ 0 }, t{ 0 } {}
		constexpr Linear(float s, float t) : s{ s }, t{ t } {}
	};

	std::ostream& operator<<(std::ostream& out, const Linear& l);
	std::wostream& operator<<(std::wostream& out, const Linear& l);
	std::istream& operator>>(std::istream& in, Linear& l);
	std::wistream& operator>>(std::wistream& in, Linear& l);
}
