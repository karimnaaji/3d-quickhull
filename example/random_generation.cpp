#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#define QUICKHULL_DEBUG
#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"

static float rand_0_1() {
    return (float) rand() / RAND_MAX;
}

int main(int argc, char** argv) {
    //srand(time(NULL));
    const int n = 200;
    const int nmeshes = 100;
    qh_vertex_t vertices[n];
    float radius = 1.0;

    qh_mesh_t meshes[nmeshes];
    for (int i = 0; i < nmeshes; ++i) {
        for (int j = 0; j < n; ++j) {
            float a0 = (rand_0_1() * M_PI * 2);
            float a1 = (rand_0_1() * M_PI * 2);
            vertices[j].z = sin(a0) * radius;
            vertices[j].x = cos(a1) * cos(a0) * radius * rand_0_1();
            vertices[j].y = sin(a1) * cos(a0) * radius * rand_0_1();
        }

        clock_t start = clock();
        meshes[i] = qh_quickhull3d(vertices, n);
        clock_t end = clock();
        printf("Time (ms): %f\n", float(end - start) / CLOCKS_PER_SEC * 1000);
    }

    std::ofstream file("hull_mesh.obj");
    if (file.is_open()) {
        int vertexOffset = 0;
        int normalOffset = 0;
        for (int i = 0; i < nmeshes; ++i) {
            qh_mesh_t m = meshes[i];
            file << " " << "\n";
            file << "o " << std::to_string(i) << "\n";
            for (int i = 0; i < m.nvertices; ++i) {
                qh_vertex_t v = m.vertices[i];
                file << "v " << v.x << " " << v.y << " " << v.z << "\n";
            }
            for (int i = 0; i < m.nnormals; ++i) {
                qh_vec3_t n = m.normals[i];
                file << "vn " << n.x << " " << n.y << " " << n.z << "\n";
            }
            for (int i = 0, j = 0; i < m.nindices; i += 3, j++) {
                file << "f ";
                file << m.indices[i+0] + 1 + vertexOffset << "//";
                file << m.normalindices[j] + 1 + normalOffset << " ";
                file << m.indices[i+1] + 1 + vertexOffset << "//";
                file << m.normalindices[j] + 1 + normalOffset << " ";
                file << m.indices[i+2] + 1 + vertexOffset << "//";
                file << m.normalindices[j] + 1 + normalOffset << "\n";
            }
            vertexOffset += m.nvertices;
            normalOffset += m.nnormals;
        }
    }

    for (int i = 0; i < nmeshes; ++i) {
        qh_free_mesh(meshes[i]);
    }

    return 0;
}
