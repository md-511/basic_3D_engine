#include "mesh.hpp"
#include "triangle.hpp"
#include <vector>

Mesh::Mesh(const std::vector<Triangle> &faces)
{
    // * copy the faces to tris
    for (Triangle t : faces)
    {
        tris.push_back(t);
    }
}
