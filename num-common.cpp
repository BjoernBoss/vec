#include "num-common.h"

/* implement the float zero comparison function */
bool num::Zero(float a, float p) {
	/* dont check for nan as nan will fail this check and thereby return false by default */
	return std::abs(a) <= num::ZeroPrecisionFactor * p;
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

num::Linear::Linear() : s(0.0f), t(0.0f) {}
num::Linear::Linear(float _s, float _t) : s(_s), t(_t) {}

std::ostream& num::operator<<(std::ostream& out, const Linear& l) {
	return (out << "s: " << l.s << "| t: " << l.t);
}
std::wostream& num::operator<<(std::wostream& out, const Linear& l) {
	return (out << L"s: " << l.s << L"| t: " << l.t);
}
std::istream& num::operator>>(std::istream& in, Linear& l) {
	char pad0 = 0, pad1 = 0, pad2 = 0, pad3 = 0, pad4 = 0;
	in >> pad0 >> pad1 >> l.s >> pad2 >> pad3 >> pad4 >> l.t;
	if (pad0 != 's' || pad1 != ':' || pad2 != '|' || pad3 != 't' || pad4 != ':')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& num::operator>>(std::wistream& in, Linear& l) {
	wchar_t pad0 = 0, pad1 = 0, pad2 = 0, pad3 = 0, pad4 = 0;
	in >> pad0 >> pad1 >> l.s >> pad2 >> pad3 >> pad4 >> l.t;
	if (pad0 != L's' || pad1 != L':' || pad2 != L'|' || pad3 != L't' || pad4 != L':')
		in.setstate(std::ios::failbit);
	return in;

}
