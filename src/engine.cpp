#include "engine.hpp"
#include "mesh.hpp"
#include "triangle.hpp"
#include "vec3d.hpp"
#include <vector>
#include <SDL2/SDL.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#define PI 3.14159265359

#define DEBUG(x) std::cout << x << std::endl
#define ERROR std::cout << "Error: " << SDL_GetError() << std::endl

// * HELPER FUNCTIONS

void fixAngle(float &angle)
{
    if (angle < 0) { angle += 2*PI; return; }
    if (angle > 2*PI) { angle -= 2*PI; return; }
}

float getMag(Vec3D &vec)
{
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

void normalize(Vec3D &vec)
{
    float mag = getMag(vec);
    if (mag <= (float)1e-6) return;
    vec.x /= mag; vec.y /= mag; vec.z /= mag;
}

void crossProd(Vec3D &l1, Vec3D &l2, Vec3D &res)
{
    res.x = l1.y*l2.z - l1.z*l2.y;
    res.y = l1.z*l2.x - l1.x*l2.z;
    res.z = l1.x*l2.y - l1.y*l2.x;
}

float dotProd(Vec3D &A, Vec3D &B)
{
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

void swap(int &a, int &b)
{
    int tmp = a;
    a = b;
    b = tmp;
}

void fillTriangle(SDL_Renderer *renderer, float x1, float y1, float x2, float y2, float x3, float y3, SDL_Color clr)
{
    SDL_Vertex vertices[3] = {
        {{x1, y1}, clr, {0.0f, 0.0f}},
        {{x2, y2}, clr, {0.0f, 0.0f}},
        {{x3, y3}, clr, {0.0f, 0.0f}}
    };

    SDL_RenderGeometry(renderer, nullptr, vertices, 3, nullptr, 0);
}

Uint8 getOpacity(float corr)
{
    // corr *= -1;
    corr += 1.0f;
    return (Uint8) (corr * 0.5f * 255);
}

// * MAIN FUNCTIONS

Engine::Engine(const int W, const int H, const int FPS, float _fov, float _zFar, float _zNear) : zFar(_zFar), zNear(_zNear), fov(_fov), screenW(W), screenH(H), fps(FPS), window(nullptr), rndr(nullptr), running(1)
{
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        ERROR;
        SDL_Quit();
        exit(1);
    }

    window = SDL_CreateWindow("Rabbit", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, screenW, screenH, 0);
    if (!window) {
        ERROR;
        SDL_Quit();
        exit(1);
    }

    rndr = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!rndr) {
        ERROR;
        SDL_Quit();
        exit(1);
    }

    SDL_BlendMode customBlendMode = SDL_ComposeCustomBlendMode(
        SDL_BLENDFACTOR_ONE, SDL_BLENDFACTOR_ZERO,    // Color blend factors
        SDL_BLENDOPERATION_ADD,                       // Color blend operation
        SDL_BLENDFACTOR_ONE, SDL_BLENDFACTOR_ZERO,    // Alpha blend factors
        SDL_BLENDOPERATION_ADD                        // Alpha blend operation
    );
    // ! Temporary added for light check
    SDL_SetRenderDrawBlendMode(rndr, SDL_BLENDMODE_BLEND);

    fov = fov * PI / 180.0f;
    r_spdX = r_spdX * PI / 180.0f;
    r_spdY = r_spdY * PI / 180.0f;
    r_spdZ = r_spdZ * PI / 180.0f;

    aspectRatio = screenW / (float) screenH;

    projectionMat[0][0] = 1.0f / aspectRatio / tan(fov*0.5f);
    projectionMat[1][1] = 1.0f / tan(fov*0.5f);
    projectionMat[2][2] = zFar / (zFar - zNear);
    projectionMat[3][2] = -zFar * zNear / (zFar - zNear);
    projectionMat[2][3] = 1.0f;

    rotationMatX[0][0] =         1.0f;
    rotationMatX[1][1] =  cos(r_spdX);
    rotationMatX[1][2] = -sin(r_spdX);
    rotationMatX[2][1] =  sin(r_spdX);
    rotationMatX[2][2] =  cos(r_spdX);

    rotationMatY[0][0] =  cos(r_spdY);
    rotationMatY[0][2] =  sin(r_spdY);
    rotationMatY[1][1] =         1.0f;
    rotationMatY[2][0] = -sin(r_spdY);
    rotationMatY[2][2] =  cos(r_spdY);

    rotationMatZ[0][0] =  cos(r_spdZ);
    rotationMatZ[0][1] = -sin(r_spdZ);
    rotationMatZ[1][0] =  sin(r_spdZ);
    rotationMatZ[1][1] =  cos(r_spdZ);
    rotationMatZ[2][2] =         1.0f;

    // std::vector<Triangle> cubeMesh = {
    //     // SOUTH
    //     {{0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
    //     {{0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {1.0f, 0.0f, 0.0f}},
    //     // EAST
    //     {{1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 1.0f}},
    //     {{1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f}, {1.0f, 0.0f, 1.0f}},
    //     // TOP
    //     {{1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f}},
    //     {{1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f}, {1.0f, 1.0f, 1.0f}},
    //     // WEST
    //     {{0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
    //     {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 1.0f}},
    //     // NORTH
    //     {{1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 1.0f}},
    //     {{1.0f, 0.0f, 1.0f}, {1.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 1.0f}},
    //     // BOTTOM
    //     {{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}},
    //     {{1.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 1.0f}},
    // };

    camera.x = 0.0f; camera.y = 0.0f; camera.z = 0.0f;

    // loadMeshFromFile("../test.txt");

    // Mesh cube = loadMesh(cubeMesh);
    Mesh cube = loadMeshFromFile("../models/teapot.obj");

    objs.push_back(cube);
}

void Engine::run(void)
{

    Uint32 prevTime = SDL_GetTicks();
    float accumulatedTime = 0.0f;

    while (running) {
        Uint32 currTime = SDL_GetTicks();

        SDL_SetRenderDrawColor(rndr, 0, 0, 0, 255);
        SDL_RenderClear(rndr);

        accumulatedTime += currTime - prevTime;
        prevTime = currTime;

        if (accumulatedTime >= (float)(1000.0/fps)) {
            accumulatedTime -= (float)(1000.0/fps);

            // r_spdX = PI;
            // r_spdY = PI * 0.3f;
            // r_spdZ = PI * 0.1f;
            r_spdX += PI / 180.0f;
            // r_spdY += PI / 540.0f;
            r_spdZ += PI / 360.0f;

            // DEBUG("RX: " << r_spdX << " RY: " << r_spdY << " RZ: " << r_spdZ);

            fixAngle(r_spdX);
            fixAngle(r_spdY);
            fixAngle(r_spdZ);
            // DEBUG(r_spdX);

            rotationMatX[0][0] =         1.0f;
            rotationMatX[1][1] =  cos(r_spdX);
            rotationMatX[1][2] = -sin(r_spdX);
            rotationMatX[2][1] =  sin(r_spdX);
            rotationMatX[2][2] =  cos(r_spdX);

            rotationMatY[0][0] =  cos(r_spdY);
            rotationMatY[0][2] =  sin(r_spdY);
            rotationMatY[1][1] =         1.0f;
            rotationMatY[2][0] = -sin(r_spdY);
            rotationMatY[2][2] =  cos(r_spdY);

            rotationMatZ[0][0] =  cos(r_spdZ);
            rotationMatZ[0][1] = -sin(r_spdZ);
            rotationMatZ[1][0] =  sin(r_spdZ);
            rotationMatZ[1][1] =  cos(r_spdZ);
            rotationMatZ[2][2] =         1.0f;

            handleEvents();
        }

        draw();

        SDL_RenderPresent(rndr);
    }

    SDL_DestroyRenderer(rndr);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void Engine::handleEvents(void)
{
    SDL_Event event;
    while (SDL_PollEvent(&event))
    {
        switch(event.type)
        {
            case SDL_QUIT:
                running = 0;
                break;
        }
    }
}

Mesh Engine::loadMesh(const std::vector<Triangle> &faces)
{
    Mesh mesh(faces);
    return mesh;
}

Mesh Engine::loadMeshFromFile(const std::string path)
{
    std::ifstream objFile(path);
    std::string line;

    std::vector<Vec3D> vertexBuffer;
    std::vector<Triangle> faceBuffer;


    while (getline(objFile, line))
    {
        std::istringstream sstr(line);
        char dump;
        float a, b, c;
        sstr >> dump >> a >> b >> c;

        if (dump == 'v') {
            Vec3D vertex(a, b, c);
            vertexBuffer.push_back(vertex);
        } else if (dump == 'f') {
            faceBuffer.push_back( Triangle( vertexBuffer[(int)a-1], vertexBuffer[(int)b-1], vertexBuffer[(int)c-1] ) );
        }
    }

    objFile.close();

    return Mesh(faceBuffer);
}

void Engine::multiplyVecMat(Vec3D &i, Vec3D &o, float (*mat)[4])
{
    o.x = i.x * mat[0][0] + i.y * mat[1][0] + i.z * mat[2][0] + 1 * mat[3][0];
    o.y = i.x * mat[0][1] + i.y * mat[1][1] + i.z * mat[2][1] + 1 * mat[3][1];
    o.z = i.x * mat[0][2] + i.y * mat[1][2] + i.z * mat[2][2] + 1 * mat[3][2];
    float z = i.x * mat[0][3] + i.y * mat[1][3] + i.z * mat[2][3] + 1 * mat[3][3];
    
    if (z != 0.0f) {
        o.x /= z; o.y /= z; o.z /= z;
    }
}

void Engine::drawTriangle(const Triangle &tri, SDL_Color clr, bool wireframe)
{
    // DEBUG("called");
    fillTriangle(rndr, tri.triangle[0].x, tri.triangle[0].y, tri.triangle[1].x, tri.triangle[1].y, tri.triangle[2].x, tri.triangle[2].y, clr);

    if (wireframe) {
        SDL_SetRenderDrawColor(rndr, 0, 0, 0, 255);
        SDL_RenderDrawLine(rndr, tri.triangle[0].x, tri.triangle[0].y, tri.triangle[1].x, tri.triangle[1].y);
        SDL_RenderDrawLine(rndr, tri.triangle[1].x, tri.triangle[1].y, tri.triangle[2].x, tri.triangle[2].y);
        SDL_RenderDrawLine(rndr, tri.triangle[2].x, tri.triangle[2].y, tri.triangle[0].x, tri.triangle[0].y);
    }
}

void Engine::draw(void)
{
    Mesh cube = objs[0];

    std::vector<std::pair<Triangle, SDL_Color>> trianglesToDraw;

    for (auto &tri : cube.tris)
    {
        Triangle triProjected, triTranslated, triRotatedX, triRotatedXY, triRotatedXYZ;

        multiplyVecMat(tri.triangle[0], triRotatedX.triangle[0], rotationMatX);
        multiplyVecMat(tri.triangle[1], triRotatedX.triangle[1], rotationMatX);
        multiplyVecMat(tri.triangle[2], triRotatedX.triangle[2], rotationMatX);

        multiplyVecMat(triRotatedX.triangle[0], triRotatedXY.triangle[0], rotationMatY);
        multiplyVecMat(triRotatedX.triangle[1], triRotatedXY.triangle[1], rotationMatY);
        multiplyVecMat(triRotatedX.triangle[2], triRotatedXY.triangle[2], rotationMatY);

        multiplyVecMat(triRotatedXY.triangle[0], triRotatedXYZ.triangle[0], rotationMatZ);
        multiplyVecMat(triRotatedXY.triangle[1], triRotatedXYZ.triangle[1], rotationMatZ);
        multiplyVecMat(triRotatedXY.triangle[2], triRotatedXYZ.triangle[2], rotationMatZ);        

        triTranslated = triRotatedXYZ;

        triTranslated.triangle[0].z += 5.0f;
        triTranslated.triangle[1].z += 5.0f;
        triTranslated.triangle[2].z += 5.0f;

        Vec3D normal, l1, l2;
        
        l1.x = triTranslated.triangle[1].x - triTranslated.triangle[0].x;
        l1.y = triTranslated.triangle[1].y - triTranslated.triangle[0].y;
        l1.z = triTranslated.triangle[1].z - triTranslated.triangle[0].z;

        l2.x = triTranslated.triangle[2].x - triTranslated.triangle[0].x;
        l2.y = triTranslated.triangle[2].y - triTranslated.triangle[0].y;
        l2.z = triTranslated.triangle[2].z - triTranslated.triangle[0].z;

        crossProd(l1, l2, normal);

        // DEBUG("Mag: " << getMag(l1));

        normalize(normal);

        Vec3D dir;       // * all pts lie on same plane so can use anyone
        dir.x = triTranslated.triangle[0].x - camera.x;
        dir.y = triTranslated.triangle[0].y - camera.y;
        dir.z = triTranslated.triangle[0].z - camera.z;

        normalize(dir);

        if (dotProd(normal, dir) < 0.0f) {

            // * Simple Directional Light

            Vec3D light(-1.0f, 0.0f, 0.0f);
            normalize(light);

            float dp = dotProd(normal, light);

            // Map [-1,1] to [255,0].... In the future store the color information in triangle or maybe in Vec3D
            Uint8 opacity = getOpacity(dp);

            multiplyVecMat(triTranslated.triangle[0], triProjected.triangle[0], projectionMat);
            multiplyVecMat(triTranslated.triangle[1], triProjected.triangle[1], projectionMat);
            multiplyVecMat(triTranslated.triangle[2], triProjected.triangle[2], projectionMat);

            // * Scale UP

            triProjected.triangle[0].x += 1.0f;
            triProjected.triangle[0].y += 1.0f;
            // triProjected.triangle[0].z += 1.0f;

            triProjected.triangle[1].x += 1.0f;
            triProjected.triangle[1].y += 1.0f;
            // triProjected.triangle[0].z += 1.0f;

            triProjected.triangle[2].x += 1.0f;
            triProjected.triangle[2].y += 1.0f;
            // triProjected.triangle[2].z += 1.0f;

            triProjected.triangle[0].x *= 0.5f * (float) screenW;
            triProjected.triangle[0].y *= 0.5f * (float) screenH;

            triProjected.triangle[1].x *= 0.5f * (float) screenW;
            triProjected.triangle[1].y *= 0.5f * (float) screenH;

            triProjected.triangle[2].x *= 0.5f * (float) screenW;
            triProjected.triangle[2].y *= 0.5f * (float) screenH;

            // drawTriangle(triProjected, {opacity, 218, 245, opacity});
            trianglesToDraw.push_back({triProjected, {opacity, 240, 124, opacity}});
            drawTriangle(triProjected, {opacity, 240, 124, opacity});
            // drawTriangle(triProjected, {149, 218, 245, opacity});
            
            // DEBUG("Opacity: " << opacity);
        }

        // ! Sort based on z

        // std::sort(trianglesToDraw.begin(), trianglesToDraw.end(), [](std::pair<Triangle, SDL_Color> &t1, std::pair<Triangle, SDL_Color> &t2) {
        //     float z1 = (t1.first.triangle[0].z + t1.first.triangle[1].z + t1.first.triangle[2].z) / 3.0f;
        //     float z2 = (t2.first.triangle[0].z + t2.first.triangle[1].z + t2.first.triangle[2].z) / 3.0f;
        //     return z1 < z2;
        // });

        // for (auto &tri : trianglesToDraw)
        // {
        //     // DEBUG(tri.second);
        //     drawTriangle(tri.first, tri.second);
        // }
    }
}
