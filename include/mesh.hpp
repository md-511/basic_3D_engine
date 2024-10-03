#pragma once

#include <vector>
#include "triangle.hpp"

class Mesh
{
public:
    std::vector<Triangle> tris;
    Mesh() {}
    Mesh(const std::vector<Triangle> &faces);
};
