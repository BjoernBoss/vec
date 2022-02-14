#pragma once

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <utility>

/* define the float number helper */
namespace num {
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
}
