# 3D Vector operations for C++
![C++](https://img.shields.io/badge/language-c%2B%2B20-blue?style=flat-square)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-brightgreen?style=flat-square)](LICENSE.txt)

Header-only library written in `C++20` to add support for various vector/line/plane functions from linear algebra. This includes testing operations, distance operations, as well as intersection functionalities.

The core types are `num::Vec<T>`, `num::Line<T>`, and `num::Plane<T>`, which each represent the corresponding concept form linear algebra. The types themselves are templated, and can be used with any `floating_point` type. For convenience, `num::Vecf` / `num::Vecd` ... are defined.

## Using the library
As this is a header-only library, no complex build setup is required. Simply clone the repository, and include `<vec/vec.h>`. The only further requirement is, that the library is compiled using `C++20`.

	$ git clone https://github.com/BjoernBoss/vec.git

## General Structure
The mathematics are based on a `right-handed system` where angles are calculated in `degrees`, and rotations are described `counterclockwise` when the corresponding axis points towards the observer.

All operations leave the corresponding object unchanged, but rather return the modified new version. For some functions there exist two variants, such as `num::Line::intersectf` and `num::Line::intersect`. These operations typically provide a `factor` version, which does not return the resulting `num::Vec` (in the case of the intersection), but rather the factors to scale the initial lines by to reach the intersection point.

For some operations, such as `num::Line::intersect`, there exist multiple variants, for example `num::Line::intersectX` / `num::Line::intersectY` / `num::Line::intersectZ`. These are typically optimized by reducing the input-space by one dimension, thereby reducing the mathematical complexity. 

All operations, which perform some form of testing, all take a precision as argument. It is used as the floating-point precision, which consideres two values identical. To prevent rounding-errors. This also extends to comparing objects, such as `num::Vec`. The `num::Vec` can, for example, be compared for being identical (i.e. all components are identical), or if two vectors match (i.e. they point into the same direction with the same magnitude, despite small imperfections).

## Example Usages

Example of computing the intersection between a line and a plane.

```C++
/* line passing through (4,4,4) and pointing parallel to the x-axis */
num::Linef _line{ num::Vecf{ 4.0f }, num::Vecf::AxisX() };

/* plane passing through the origin, extended into X/Y and Y/Z */
num::Planef _plane{ num::Vecf{ 1.0f, 1.0f, 0.0f }, num::Vecf{ 0.0f, 1.0f, 1.0f } };

/* intersection-point of the line and the plane */
num::Vecf _point = _plane.intersect(_line);
std::cout << _point << std::endl; // prints: (0, 4, 4)
```

Example of computing the intersection between two planes.

```C++
/* plane parallel to Y-Z plane and passing through (2, 0, 0) */
num::Planed _planeyz = num::Planed::AxisX(2.0);

/* plane parallel to X-Z plane and passing through (0, 2, 0) */
num::Planed _planexz = num::Planed::AxisY(2.0);

/* intersection-line of the two planes */
num::Lined _line = _planeyz.intersect(_planexz);
std::cout << _line << std::endl; // prints: (2, 2, 0) -> (0, 0, 1)
```
