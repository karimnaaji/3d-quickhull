#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#define QUICKHULL_IMPLEMENTATION
#include "../quickhull.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <time.h>

#define ITERATION 100

void display_result(float time, float avg, unsigned int n, std::string file) {
    std::cout << std::setw(20) << std::setfill(' ') << file;
    std::cout.precision(4);
    std::cout << std::fixed;
    std::cout << std::setw(20) << std::setfill(' ') << time << " ms (time)";
    std::cout << std::setw(20) << std::setfill(' ') << avg << " ms (avg)";
    std::cout << std::setw(20) << std::setfill(' ') << "n: " << n;
    std::cout << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << std::endl;
        std::cout << "\tqh [obj_file|text_file]" << std::endl;
        std::cout << "\tobj_file: regular OBJ wavefront file format" << std::endl;
        std::cout << "\ttext_file:" << std::endl;
        std::cout << "\t\tvertex_0.x vertex_0.y vetex_0.y" << std::endl;
        std::cout << "\t\tvertex_1.x vertex_1.y vetex_1.y" << std::endl;
        std::cout << "\t\t[...]" << std::endl;
        std::cout << "\t\tvertex_n.x vertex_n.y vetex_n.y" << std::endl;
        return EXIT_FAILURE;
    }

    std::string path = std::string(argv[1]);
    std::string file = path.substr(path.find_last_of("/\\") + 1);
    std::string hull = "hull_" + file;
    if (file.substr(file.find_last_of(".") + 1) == "obj") {
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        std::string err;
        bool ret = tinyobj::LoadObj(shapes, materials, err, path.c_str(), NULL);
        if (!err.empty()) { std::cerr << err << std::endl; }
        if (!ret) { return EXIT_FAILURE; }
        std::vector<qh_vertex_t> vertices;
        vertices.reserve(shapes[0].mesh.positions.size() / 3);
        for (size_t f = 0; f < shapes[0].mesh.positions.size(); f += 3) {
            qh_vertex_t vertex;
            vertex.x = shapes[0].mesh.positions[f + 0];
            vertex.y = shapes[0].mesh.positions[f + 1];
            vertex.z = shapes[0].mesh.positions[f + 2];
            bool unique = true;
            for (const auto& v : vertices) {
                if (v.x == vertex.x && v.y == vertex.y && v.z == vertex.z) {
                    unique = false;
                    break;
                }
            }
            if (unique) {
                vertices.push_back(vertex);
            }
        }
        float avg = 0.f;
        float time = 0.f;
        int iteration = 0;
        qh_mesh_t m;
        do {
            clock_t start = clock();
            m = qh_quickhull3d(vertices.data(), vertices.size());
            clock_t end = clock();
            qh_free_mesh(m);
            time = (float)(end - start) / CLOCKS_PER_SEC * 1000;
            avg += time;
        } while (++iteration < ITERATION);

        display_result(time, avg / iteration, vertices.size(), file);
        m = qh_quickhull3d(vertices.data(), vertices.size());
        qh_mesh_export(&m, hull.c_str());
        qh_free_mesh(m);
    } else {
        int n = 0;
        std::vector<qh_vertex_t> vertices;
        std::ifstream ifs(path, std::ios::in);
        if (ifs) {
            std::ifstream count(path, std::ios::in);
            n = std::count(std::istreambuf_iterator<char>(count), std::istreambuf_iterator<char>(), '\n');
            vertices.resize(n);
            qh_vertex_t* v = vertices.data();
            int j = 0;
            while (!ifs.eof() && j != n) {
                ifs >> v->x >> v->y >> v->z;
                std::cout << v->x << " " << v->y << " " << v->z << std::endl;
                j++;
                v++;
            }
            ifs.close();

            float avg = 0.f;
            float time = 0.f;
            int iteration = 0;
            qh_mesh_t m;
            do {
                clock_t start = clock();
                qh_mesh_t m = qh_quickhull3d(vertices.data(), vertices.size());
                clock_t end = clock();
                qh_free_mesh(m);
                time = (float)(end - start) / CLOCKS_PER_SEC * 1000;
                avg += time;
            } while (++iteration < ITERATION);

            display_result(time, avg / iteration, vertices.size(), file);
            m = qh_quickhull3d(vertices.data(), vertices.size());
            qh_mesh_export(&m, hull.c_str());
            qh_free_mesh(m);
        }
    }

    return EXIT_SUCCESS;
}
