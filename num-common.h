#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <utility>

namespace num {
	/*
	*	- right-handed system
	*	- counterclockwise rotations when the corresponding axis points towards the observer
	*	- angles are calculated in degrees
	*/

	static constexpr double Pi = 3.1415926536;
	static constexpr double ZeroPrecisionFactor = 0.01;

	template <class Type> struct Precision;
	template <> struct Precision<float> {
		static constexpr float Def = 0.00001f;
	};
	template <> struct Precision<double> {
		static constexpr double Def = 0.00000001;
	};

	/* float abs-function (not using std implementation to allow for constexpr) */
	template <class Type>
	constexpr Type Abs(Type v) {
		return (v < 0 ? -v : v);
	}

	/* check if number can be considered zero */
	template <class Type>
	constexpr bool Zero(Type a, Type p = num::Precision<Type>::Def) {
		/* dont check for nan as nan will fail this check and thereby return false by default */
		return num::Abs(a) <= Type(num::ZeroPrecisionFactor * p);
	}

	/* compare the values for equality, given the corresponding precision */
	template <class Type>
	constexpr bool Cmp(Type a, Type b, Type p = num::Precision<Type>::Def) {
		if (std::isnan(a) || std::isnan(b))
			return false;
		if (a == 0)
			return num::Zero(b);
		if (b == 0)
			return num::Zero(a);
		const Type _a = num::Abs(a);
		const Type _b = num::Abs(b);
		return num::Abs(a - b) <= std::min(_a, _b) * p;
	}

	template <class Type>
	constexpr Type ToRadian(Type deg) {
		return deg * Type(num::Pi / 180);
	}

	template <class Type>
	constexpr Type ToDegree(Type deg) {
		return deg * Type(180 / num::Pi);
	}

	template <class Type>
	constexpr Type ToAngle(Type x, Type y) {
		Type deg = num::ToDegree(std::atan2(x, y));
		if (deg < 0)
			deg += 360;
		return deg;
	}

	/* compute the angle to add to [base] to reach [test] in degrees */
	template <class Type>
	constexpr Type AngleDiff(Type base, Type test) {
		Type diff = test - base;
		if (diff <= -180)
			diff += 360;
		else if (diff > 180)
			diff -= 360;
		return diff;
	}

	/* compute the absolute difference between [base] and [test] in degrees */
	template <class Type>
	constexpr Type AngleAbs(Type base, Type test) {
		Type diff = num::Abs(test - base);
		if (diff > 180)
			diff = 360 - diff;
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
}
