#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <time.h>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#define FAST

#ifndef FAST
#define QUICKHULL_DEBUG
#endif
#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"

static float rand_0_1() {
    return (float) rand() / RAND_MAX;
}

void dump(qh_vertex_t* vertices, unsigned int n, unsigned int failurestep) {
    std::ofstream file("dump" + std::to_string(failurestep) + ".txt");
    if (file.is_open()) {
        for (int i = 0; i < n; ++i) {
            qh_vertex_t v = vertices[i];
            file << v.x << " " << v.y << " " << v.z << "\n";
        }
    }
}

std::vector<qh_vec3_t> circle(int samples) {
    std::vector<qh_vec3_t> points;
    double step = (2.0 * M_PI) / samples;
    for (int i = 0; i < samples; ++i) {
        points.push_back({(float)cos(i * step), (float)sin(i * step), 0.0});
    }
    return points;
}

bool load_obj(std::string path, std::vector<qh_vertex_t>& vertices) {
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    bool ret = tinyobj::LoadObj(shapes, materials, err, path.c_str(), NULL);
    if (!ret) {
        return false;
    }
    if (!err.empty()) {
        printf("tinyobj::LoadObj error: %s\n", err.c_str());
    }
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
    return true;
}

void test_obj_file(std::string path, std::string name) {
    std::vector<qh_vertex_t> vertices;
    float epsilon;
    qh_context_t context;
    if (!load_obj(path, vertices)) {
        printf("Loading %s failed\n", path.c_str());
        REQUIRE(false);
    }
    printf("Loaded %ld vertices [%s]\n", vertices.size(), path.c_str());
    epsilon = qh__compute_epsilon(vertices.data(), vertices.size());
    qh__init_context(&context, vertices.data(), vertices.size());
    qh__remove_vertex_duplicates(&context, epsilon);
    qh__build_tetrahedron(&context, epsilon);
    unsigned int failurestep = 0;
#ifdef FAST
    qh__build_hull(&context, epsilon);
#else
    qh__build_hull(&context, epsilon, -1, &failurestep);
#endif
    int valid = qh__test_hull(&context, epsilon, 0);
    qh_mesh_t m = qh_quickhull3d(vertices.data(), vertices.size());
    if (!valid) {
        qh_mesh_export(&m, std::string(name + "_failure.obj").c_str());
        dump(vertices.data(), vertices.size(), failurestep);
    } else {
        qh_mesh_export(&m, std::string(name + "_hull.obj").c_str());
        // dump(vertices.data(), vertices.size(), 0);
    }
    REQUIRE(valid);
    for (int i = 0; i < m.nindices; ++i) {
        REQUIRE(m.indices[i] >= 0);
        REQUIRE(m.indices[i] < m.nvertices);
        REQUIRE(!std::isnan(m.vertices[m.indices[i]].x));
        REQUIRE(!std::isnan(m.vertices[m.indices[i]].y));
        REQUIRE(!std::isnan(m.vertices[m.indices[i]].z));
    }
    for (int i = 0; i < m.nnormals; ++i) {
        REQUIRE(!std::isnan(m.normals[m.normalindices[i]].x));
        REQUIRE(!std::isnan(m.normals[m.normalindices[i]].y));
        REQUIRE(!std::isnan(m.normals[m.normalindices[i]].z));
    }
    qh_free_mesh(m);
    qh__free_context(&context);
}

// -----------------------------------------------------------------------------

TEST_CASE("200 meshes on a sphere", "quickhull.h") {
    return;
    const int n = 200;
    const int nmeshes = 100;
    qh_vertex_t vertices[n];
    float radius = 1.0;
    for (int i = 0; i < nmeshes; ++i) {
        for (int j = 0; j < n; ++j) {
            float a0 = (rand_0_1() * M_PI * 2);
            float a1 = (rand_0_1() * M_PI * 2);
            vertices[j].z = sin(a0) * radius;
            vertices[j].x = cos(a1) * cos(a0) * radius * rand_0_1();
            vertices[j].y = sin(a1) * cos(a0) * radius * rand_0_1();
        }
        qh_context_t context;
        float epsilon;

        epsilon = qh__compute_epsilon(vertices, n);
        qh__init_context(&context, vertices, n);
        qh__remove_vertex_duplicates(&context, epsilon);
        qh__build_tetrahedron(&context, epsilon);
        unsigned int failurestep = 0;
#ifdef FAST
        qh__build_hull(&context, epsilon);
#else
        qh__build_hull(&context, epsilon, -1, &failurestep);
#endif
        int valid = qh__test_hull(&context, epsilon, 0);
        if (!valid) {
            qh_mesh_t m = qh_quickhull3d(vertices, n);
            std::cout << "Saving failure mesh" << std::endl;
            qh_mesh_export(&m, std::string("failure_mesh" + std::to_string(i) + ".obj").c_str());
            dump(vertices, n, failurestep);
            qh_free_mesh(m);
        }
        REQUIRE(valid);
        qh__free_context(&context);
    }
}

TEST_CASE("Two circle shapes", "quickhull.h") {
    for (int i = 10; i < 1000; ++i) {
        auto c0 = circle(6 + i / 10);
        auto c1 = circle(5 + i);
        for (int j = 0; j < c1.size(); ++j) {
            c1[j].z += 1.0;
        }
        auto v(c0);
        v.insert(v.end(), c1.begin(), c1.end());
        qh_vertex_t* vertices = v.data();
        unsigned int n = v.size();
        qh_context_t context;
        float epsilon = qh__compute_epsilon(vertices, n);
        qh__init_context(&context, vertices, n);
        qh__remove_vertex_duplicates(&context, epsilon);
        qh__build_tetrahedron(&context, epsilon);
        unsigned int failurestep = 0;
#ifdef FAST
        qh__build_hull(&context, epsilon);
#else
        qh__build_hull(&context, epsilon, -1, &failurestep);
#endif
        int valid = qh__test_hull(&context, epsilon, 0);
        if (!valid) {
            qh_mesh_t m = qh_quickhull3d(vertices, n);
            std::cout << "Saving failure mesh" << std::endl;
            qh_mesh_export(&m, std::string("failure_mesh_two_circles" + std::to_string(i) + ".obj").c_str());
            dump(vertices, n, i);
            qh_free_mesh(m);
        }
        REQUIRE(valid);
        qh__free_context(&context);
    }
}

TEST_CASE("19295.24642.16.obj", "quickhull.h") {
    test_obj_file("models/19295.24642.16.obj", "tile0");
}

TEST_CASE("19296.24630.16.obj", "quickhull.h") {
    test_obj_file("models/19296.24630.16.obj", "tile0");
}

TEST_CASE("19296.24641.16.obj", "quickhull.h") {
    test_obj_file("models/19296.24641.16.obj", "tile1");
}

#if 0 // fails with nan
TEST_CASE("19294.24642.16.obj", "quickhull.h") {
    test_obj_file("models/19294.24642.16.obj", "tile3");
}
#endif

TEST_CASE("19292.24642.16.obj", "quickhull.h") {
    test_obj_file("models/19292.24642.16.obj", "tile4");
}

TEST_CASE("19294.24640.16.obj", "quickhull.h") {
    test_obj_file("models/19294.24640.16.obj", "tile5");
}

TEST_CASE("19294.24644.16.obj", "quickhull.h") {
    test_obj_file("models/19294.24644.16.obj", "tile6");
}

TEST_CASE("column.obj", "quickhull.h") {
    test_obj_file("models/column.obj", "column");
}

TEST_CASE("bunny.obj", "quickhull.h") {
    test_obj_file("models/bunny.obj", "bunny");
}

TEST_CASE("suzanne.obj", "quickhull.h") {
    test_obj_file("models/suzanne.obj", "suzanne");
}

TEST_CASE("cube.obj", "quickhull.h") {
    test_obj_file("models/cube.obj", "cube");
}

TEST_CASE("tree.obj", "quickhull.h") {
    test_obj_file("models/tree.obj", "tree");
}

TEST_CASE("platform.obj", "quickhull.h") {
    test_obj_file("models/platform.obj", "platform");
}

TEST_CASE("tree1b_lod1_2.obj", "quickhull.h") {
    test_obj_file("models/tree1b_lod1_2.obj", "tree_lod");
}

TEST_CASE("sponza_cooked.obj", "quickhull.h") {
    test_obj_file("models/sponza_cooked.obj", "sponza");
}

TEST_CASE("banner.obj", "quickhull.h") {
    test_obj_file("models/banner.obj", "banner");
}

TEST_CASE("sphere.obj", "quickhull.h") {
    test_obj_file("models/sphere.obj", "sphere");
}

