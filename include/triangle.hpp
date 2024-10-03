#pragma once

#include "vec3d.hpp"

class Triangle
{
public:
    Vec3D triangle[3];
    Triangle() {}
    Triangle(const Vec3D &a, const Vec3D &b, const Vec3D &c) : triangle{a, b, c} {}
};
