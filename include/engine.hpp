#pragma once

#include <SDL2/SDL.h>
#include <vector>
#include <string>
#include "mesh.hpp"
#include "vec3d.hpp"

class Engine
{
    float projectionMat[4][4] = {0.0f};
    float  rotationMatX[4][4] = {0.0f};
    float  rotationMatY[4][4] = {0.0f};
    float  rotationMatZ[4][4] = {0.0f};

    std::vector<Mesh> objs;
    SDL_Window *window;
    SDL_Renderer *rndr;
    const int screenW, screenH, fps;
    float fov, aspectRatio, zFar, zNear;

    // ! Temporary
    float r_spdX = 0.0f, r_spdY = 0.0f, r_spdZ = 0.0f;
    Vec3D camera;

    int running;

    Mesh loadMesh(const std::vector<Triangle> &faces);
    Mesh loadMeshFromFile(const std::string path);

    void multiplyVecMat(Vec3D &i, Vec3D &o, float (*mat)[4]);
    void draw(void);

    // ! Maybe remove drawTraingle from the class and keep as a helper function in engine.cpp
    void drawTriangle(const Triangle &tri, SDL_Color clr, bool wireframe=false);

public:
    Engine(const int W=1000, const int H=800, const int FPS=60, float _fov=75, float _zFar=1000.0f, float _zNear=0.1f);
    void run(void);
    void handleEvents(void);
};
