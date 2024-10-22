#pragma once

#include "num-common.h"
#include "num-vec.h"
#include "num-line.h"
#include "num-plane.h"

constexpr num::Line num::Vec::line(const num::Vec& p) const {
	return num::Line{ *this, p - *this };
}

constexpr num::Plane num::Vec::plane(const num::Vec& p0, const num::Vec& p1) const {
	return num::Plane{ *this, p0 - *this, p1 - *this };
}

std::ostream& operator<<(std::ostream& out, const num::Vec& v) {
	return (out << v.x << ", " << v.y << ", " << v.z);
}
std::wostream& operator<<(std::wostream& out, const num::Vec& v) {
	return (out << v.x << L", " << v.y << L", " << v.z);
}
std::istream& operator>>(std::istream& in, num::Vec& v) {
	char pad0 = 0, pad1 = 0;
	in >> v.x >> pad0 >> v.y >> pad1 >> v.z;
	if (pad0 != ',' || pad1 != ',')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& operator>>(std::wistream& in, num::Vec& v) {
	wchar_t pad0 = 0, pad1 = 0;
	in >> v.x >> pad0 >> v.y >> pad1 >> v.z;
	if (pad0 != L',' || pad1 != L',')
		in.setstate(std::ios::failbit);
	return in;
}

std::ostream& operator<<(std::ostream& out, const num::Line& l) {
	return (out << l.o << " -> " << l.d);
}
std::wostream& operator<<(std::wostream& out, const num::Line& l) {
	return (out << l.o << L" -> " << l.d);
}
std::istream& operator>>(std::istream& in, num::Line& l) {
	char pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != '-' || pad1 != '>')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& operator>>(std::wistream& in, num::Line& l) {
	wchar_t pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != L'-' || pad1 != L'>')
		in.setstate(std::ios::failbit);
	return in;
}


std::ostream& operator<<(std::ostream& out, const num::Plane& p) {
	return (out << p.o << " -> " << p.a << " | " << p.b);
}
std::wostream& operator<<(std::wostream& out, const num::Plane& p) {
	return (out << p.o << L" -> " << p.a << L" | " << p.b);
}
std::istream& operator>>(std::istream& in, num::Plane& p) {
	char pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != '-' || pad1 != '>' || pad2 != '|')
		in.setstate(std::ios::failbit);
	return in;
}
std::wistream& operator>>(std::wistream& in, num::Plane& p) {
	wchar_t pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != L'-' || pad1 != L'>' || pad2 != L'|')
		in.setstate(std::ios::failbit);
	return in;

}
