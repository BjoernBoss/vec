/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright (c) 2024 Bjoern Boss Henrichsen */
#pragma once

#include <istream>
#include <ostream>

#include "num-common.h"
#include "num-vec.h"
#include "num-line.h"
#include "num-plane.h"

namespace num {
	using Constf = num::Const<float>;
	using Constd = num::Const<double>;

	using Vecf = num::Vec<float>;
	using Vecd = num::Vec<double>;

	using Linearf = num::Linear<float>;
	using Lineard = num::Linear<double>;

	using Linef = num::Line<float>;
	using Lined = num::Line<double>;

	using Planef = num::Plane<float>;
	using Planed = num::Plane<double>;
}

template <std::floating_point Type>
constexpr num::Line<Type> num::Vec<Type>::line(const num::Vec<Type>& p) const {
	return num::Line<Type>{ *this, p - *this };
}

template <std::floating_point Type>
constexpr num::Plane<Type> num::Vec<Type>::plane(const num::Vec<Type>& p0, const num::Vec<Type>& p1) const {
	return num::Plane<Type>{ *this, p0 - *this, p1 - *this };
}

template <class Type>
std::ostream& operator<<(std::ostream& out, const num::Vec<Type>& v) {
	return (out << '(' << v.x << ", " << v.y << ", " << v.z << ')');
}
template <class Type>
std::wostream& operator<<(std::wostream& out, const num::Vec<Type>& v) {
	return (out << L'(' << v.x << L", " << v.y << L", " << v.z << L')');
}
template <class Type>
std::istream& operator>>(std::istream& in, num::Vec<Type>& v) {
	char pad0 = 0, pad1 = 0, pad2 = 0, pad3 = 0;
	in >> pad0 >> v.x >> pad1 >> v.y >> pad2 >> v.z >> pad3;
	if (pad0 != '(' || pad1 != ',' || pad2 != ',' || pad3 != ')')
		in.setstate(std::ios::failbit);
	return in;
}
template <class Type>
std::wistream& operator>>(std::wistream& in, num::Vec<Type>& v) {
	wchar_t pad0 = 0, pad1 = 0, pad2 = 0, pad3 = 0;
	in >> pad0 >> v.x >> pad1 >> v.y >> pad2 >> v.z >> pad3;
	if (pad0 != L'(' || pad1 != L',' || pad2 != L',' || pad3 != L')')
		in.setstate(std::ios::failbit);
	return in;
}

template <class Type>
std::ostream& operator<<(std::ostream& out, const num::Line<Type>& l) {
	return (out << l.o << " -> " << l.d);
}
template <class Type>
std::wostream& operator<<(std::wostream& out, const num::Line<Type>& l) {
	return (out << l.o << L" -> " << l.d);
}
template <class Type>
std::istream& operator>>(std::istream& in, num::Line<Type>& l) {
	char pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != '-' || pad1 != '>')
		in.setstate(std::ios::failbit);
	return in;
}
template <class Type>
std::wistream& operator>>(std::wistream& in, num::Line<Type>& l) {
	wchar_t pad0 = 0, pad1 = 0;
	in >> l.o >> pad0 >> pad1 >> l.d;
	if (pad0 != L'-' || pad1 != L'>')
		in.setstate(std::ios::failbit);
	return in;
}

template <class Type>
std::ostream& operator<<(std::ostream& out, const num::Plane<Type>& p) {
	return (out << p.o << " -> " << p.a << " | " << p.b);
}
template <class Type>
std::wostream& operator<<(std::wostream& out, const num::Plane<Type>& p) {
	return (out << p.o << L" -> " << p.a << L" | " << p.b);
}
template <class Type>
std::istream& operator>>(std::istream& in, num::Plane<Type>& p) {
	char pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != '-' || pad1 != '>' || pad2 != '|')
		in.setstate(std::ios::failbit);
	return in;
}
template <class Type>
std::wistream& operator>>(std::wistream& in, num::Plane<Type>& p) {
	wchar_t pad0 = 0, pad1 = 0, pad2 = 0;
	in >> p.o >> pad0 >> pad1 >> p.a >> pad2 >> p.b;
	if (pad0 != L'-' || pad1 != L'>' || pad2 != L'|')
		in.setstate(std::ios::failbit);
	return in;

}
