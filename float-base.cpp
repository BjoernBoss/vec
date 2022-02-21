#include "float-base.h"

/* implement the float zero comparison function */
bool num::Zero(float a, float p) {
	/* dont check for nan as nan will fail this check and thereby return false by default */
	return std::abs(a) <= ZeroPrecisionFactor * p;
}

/* implement the float comparison function */
bool num::Cmp(float a, float b, float p) {
	if (std::isnan(a) || std::isnan(b))
		return false;
	if (a == 0.0f)
		return Zero(b);
	if (b == 0.0f)
		return Zero(a);
	const float _a = std::abs(a);
	const float _b = std::abs(b);
	return std::abs(a - b) <= std::min(_a, _b) * p;
}

/* implement the angle conversion functions */
float num::ToAngle(float x, float y) {
	float deg = ToDegree(std::atan2(x, y));
	if (deg < 0)
		deg += 360.0f;
	return deg;
}

/* implement the angle comparison functions */
float num::AngleDiff(float base, float test) {
	float diff = test - base;
	if (diff <= -180.0f)
		diff += 360.0f;
	else if (diff > 180.0f)
		diff -= 360.0f;
	return diff;
}
float num::AngleAbs(float base, float test) {
	float diff = std::abs(test - base);
	if (diff > 180.0f)
		diff = 360.0f - diff;
	return diff;
}
