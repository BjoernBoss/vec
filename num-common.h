#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <utility>
#include <concepts>

namespace num {
	/*
	*	- right-handed system
	*	- counterclockwise rotations when the corresponding axis points towards the observer
	*	- angles are calculated in degrees
	*/

	template <std::floating_point Type> struct Const;
	template <> struct Const<float> {
		static constexpr float Precision = 0.00001f;
		static constexpr float Pi = 3.14159265f;
		static constexpr float ZeroPrecisionFactor = 0.01f;
	};
	template <> struct Const<double> {
		static constexpr double Precision = 0.00000001;
		static constexpr double Pi = 3.141592653589793;
		static constexpr double ZeroPrecisionFactor = 0.01;
	};

	/* float abs-function (not using std implementation to allow for constexpr) */
	template <std::floating_point Type>
	constexpr Type Abs(Type v) {
		return (v < 0 ? -v : v);
	}

	/* check if number can be considered zero */
	template <std::floating_point Type>
	constexpr bool Zero(Type a, Type p = num::Const<Type>::Precision) {
		/* dont check for nan as nan will fail this check and thereby return false by default */
		return num::Abs(a) <= (num::Const<Type>::ZeroPrecisionFactor * p);
	}

	/* compare the values for equality, given the corresponding precision */
	template <std::floating_point Type>
	bool Cmp(Type a, Type b, Type p = num::Const<Type>::Precision) {
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

	template <std::floating_point Type>
	constexpr Type ToRadian(Type deg) {
		return deg * (num::Const<Type>::Pi / 180);
	}

	template <std::floating_point Type>
	constexpr Type ToDegree(Type deg) {
		return deg * (180 / num::Const<Type>::Pi);
	}

	template <std::floating_point Type>
	constexpr Type ToAngle(Type x, Type y) {
		Type deg = num::ToDegree(std::atan2(x, y));
		if (deg < 0)
			deg += 360;
		return deg;
	}

	/* compute the angle to add to [base] to reach [test] in degrees */
	template <std::floating_point Type>
	constexpr Type AngleDiff(Type base, Type test) {
		Type diff = test - base;
		if (diff <= -180)
			diff += 360;
		else if (diff > 180)
			diff -= 360;
		return diff;
	}

	/* compute the absolute difference between [base] and [test] in degrees */
	template <std::floating_point Type>
	constexpr Type AngleAbs(Type base, Type test) {
		Type diff = num::Abs(test - base);
		if (diff > 180)
			diff = 360 - diff;
		return diff;
	}

	template <std::floating_point Type>
	struct Linear {
	public:
		Type s;
		Type t;

	public:
		constexpr Linear() : s{ 0 }, t{ 0 } {}
		constexpr Linear(float s, float t) : s{ s }, t{ t } {}
	};
}
