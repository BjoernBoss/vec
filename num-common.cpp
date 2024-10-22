#include "num-common.h"

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
