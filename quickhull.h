//
// LICENCE:
//  The MIT License (MIT)
//
//  Copyright (c) 2016 Karim Naaji, karim.naaji@gmail.com
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all
//  copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE
//
// REFERENCES:
//  [1] http://box2d.org/files/GDC2014/DirkGregorius_ImplementingQuickHull.pdf
//  [2] http://www.cs.smith.edu/~orourke/books/compgeom.html
//  [3] http://www.flipcode.com/archives/The_Half-Edge_Data_Structure.shtml
//  [4] http://doc.cgal.org/latest/HalfedgeDS/index.html
//  [5] http://thomasdiewald.com/blog/?p=1888
//  [6] https://fgiesen.wordpress.com/2012/02/21/half-edge-based-mesh-representations-theory/
//
// HOWTO:
//  #define QUICKHULL_IMPLEMENTATION
//  #define QUICKHULL_DEBUG // Only if assertions need to be checked
//  #include "quickhull.h"
//
// HISTORY:
//  - version 1.0: Initial
//
// TODO:
//  - use float* from public interface
//  - reduce memory usage

#ifndef QUICKHULL_H
#define QUICKHULL_H

// ------------------------------------------------------------------------------------------------
// QUICKHULL PUBLIC API
//

typedef struct qh_vertex {
    union {
        float v[3];
        struct {
            float x;
            float y;
            float z;
        };
    };
} qh_vertex_t;

typedef qh_vertex_t qh_vec3_t;

typedef struct qh_mesh {
    qh_vertex_t* vertices;
    qh_vec3_t* normals;
    unsigned int* indices;
    unsigned int* normalindices;
    size_t nindices;
    size_t nvertices;
    size_t nnormals;
} qh_mesh_t;

qh_mesh_t qh_quickhull3d(qh_vertex_t* vertices, unsigned int nvertices);

void qh_free_mesh(qh_mesh_t mesh);

//
// END QUICKHULL PUBLIC API
// ------------------------------------------------------------------------------------------------

#endif // QUICKHULL_H

#ifdef QUICKHULL_IMPLEMENTATION

#include <math.h> // sqrt & fabs

// Quickhull helpers, define your own if needed
#ifndef QUICKHULL_HELPERS
#include <stdlib.h> // malloc, free, realloc
#define QUICKHULL_HELPERS 1
#define QH_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define QH_REALLOC(T, P, N) ((T*)realloc(P, sizeof(T) * N))
#define QH_FREE(T) free(T)
#define QH_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#ifdef QUICKHULL_DEBUG
#define QH_ASSERT(STMT) if (!(STMT)) { *(int *)0 = 0; }
#else
#define QH_ASSERT(STMT)
#endif // QUICKHULL_DEBUG
#endif // QUICKHULL_HELPERS

#ifndef QH_FLT_MAX
#define QH_FLT_MAX 1e+37F
#endif

#ifndef QH_VERTEX_SET_SIZE
#define QH_VERTEX_SET_SIZE 128
#endif

typedef long qh_index_t;

typedef struct qh_half_edge {
    qh_index_t opposite_he;     // index of the opposite half edge
    qh_index_t next_he;         // index of the next half edge
    qh_index_t previous_he;     // index of the previous half edge
    qh_index_t he;              // index of the current half edge
    qh_index_t to_vertex;       // index of the next vertex
    qh_index_t adjacent_face;   // index of the ajacent face
} qh_half_edge_t;

typedef struct qh_index_set {
    qh_index_t* indices;
    size_t size;
    size_t capacity;
} qh_index_set_t;

typedef struct qh_face {
    qh_vec3_t normal;
    qh_index_t edges[3];
    qh_index_t face;
    float sdist;
    qh_index_set_t iset;
    qh_vertex_t centroid;
    char valid;
    int visitededges;
} qh_face_t;

typedef struct qh_index_stack {
    qh_index_t* begin;
    size_t size;
} qh_index_stack_t;

typedef struct qh_context {
    qh_vertex_t* vertices;
    qh_half_edge_t* edges;
    qh_face_t* faces;
    size_t nedges;
    size_t nvertices;
    size_t nfaces;
    qh_index_stack_t facestack;
    qh_index_stack_t scratch;
    qh_index_stack_t horizonedges;
    qh_index_stack_t newhorizonedges;
} qh_context_t;

void qh__find_6eps(qh_vertex_t* vertices, unsigned int nvertices, qh_index_t* eps)
{
    qh_vertex_t* ptr = vertices;

    float minxy = +QH_FLT_MAX;
    float minxz = +QH_FLT_MAX;
    float minyz = +QH_FLT_MAX;

    float maxxy = -QH_FLT_MAX;
    float maxxz = -QH_FLT_MAX;
    float maxyz = -QH_FLT_MAX;

    int i = 0;
    for (i = 0; i < 6; ++i) {
        eps[i] = 0;
    }

    for (i = 0; i < nvertices; ++i) {
        if (ptr->z < minxy) {
            eps[0] = i;
            minxy = ptr->z;
        }
        if (ptr->y < minxz) {
            eps[1] = i;
            minxz = ptr->y;
        }
        if (ptr->x < minyz) {
            eps[2] = i;
            minyz = ptr->x;
        }
        if (ptr->z > maxxy) {
            eps[3] = i;
            maxxy = ptr->z;
        }
        if (ptr->y > maxxz) {
            eps[4] = i;
            maxxz = ptr->y;
        }
        if (ptr->x > maxyz) {
            eps[5] = i;
            maxyz = ptr->x;
        }
        ptr++;
    }
}

float qh__vertex_segment_length2(qh_vertex_t* p, qh_vertex_t* a, qh_vertex_t* b)
{
    float dx = b->x - a->x;
    float dy = b->y - a->y;
    float dz = b->z - a->z;

    float d = dx * dx + dy * dy + dz * dz;

    float x = a->x;
    float y = a->y;
    float z = a->z;

    if (d != 0) {
        float t = ((p->x - a->x) * dx +
            (p->y - a->y) * dy +
            (p->z - a->z) * dz) / d;

        if (t > 1) {
            x = b->x;
            y = b->y;
            z = b->z;
        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
            z += dz * t;
        }
    }

    dx = p->x - x;
    dy = p->y - y;
    dz = p->z - z;

    return dx * dx + dy * dy + dz * dz;
}

void qh__vec3_sub(qh_vec3_t* a, qh_vec3_t* b)
{
    a->x -= b->x;
    a->y -= b->y;
    a->z -= b->z;
}

void qh__vec3_add(qh_vec3_t* a, qh_vec3_t* b)
{
    a->x += b->x;
    a->y += b->y;
    a->z += b->z;
}

void qh__vec3_multiply(qh_vec3_t* a, float v)
{
    a->x *= v;
    a->y *= v;
    a->z *= v;
}

float qh__vec3_length2(qh_vec3_t* v)
{
    return v->x * v->x + v->y * v->y + v->z * v->z;
}

float qh__vec3_dot(qh_vec3_t* v1, qh_vec3_t* v2)
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

void qh__vec3_normalize(qh_vec3_t* v)
{
    qh__vec3_multiply(v, 1.f / sqrt(qh__vec3_length2(v)));
}

void qh__find_2dps_6eps(qh_vertex_t* vertices, qh_index_t* eps, int* ii, int* jj)
{
    int i, j;
    float max = -QH_FLT_MAX;

    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) {
            qh_vertex_t d;
            float d2;

            if (i == j) {
                continue;
            }

            d = vertices[eps[i]];
            qh__vec3_sub(&d, &vertices[eps[j]]);
            d2 = qh__vec3_length2(&d);

            if (d2 > max) {
                *ii = i;
                *jj = j;
                max = d2;
            }
        }
    }
}

qh_vec3_t qh__vec3_cross(qh_vec3_t* v1, qh_vec3_t* v2)
{
    qh_vec3_t cross;

    cross.x = v1->y * v2->z - v1->z * v2->y;
    cross.y = v1->z * v2->x - v1->x * v2->z;
    cross.z = v1->x * v2->y - v1->y * v2->x;

    return cross;
}

qh_vertex_t qh__face_centroid(qh_index_t vertices[3], qh_context_t* context)
{
    qh_vertex_t centroid;
    int i;

    centroid.x = centroid.y = centroid.z = 0.0;
    for (i = 0; i < 3; ++i) {
        qh__vec3_add(&centroid, context->vertices + vertices[i]);
    }

    qh__vec3_multiply(&centroid, 1.0 / 3.0);

    return centroid;
}

float qh__dist_point_plane(qh_vertex_t* v, qh_vec3_t* normal, float sdist)
{
    return fabs(qh__vec3_dot(v, normal) - sdist);
}

void qh__init_half_edge(qh_half_edge_t* half_edge) {
    half_edge->adjacent_face = -1;
    half_edge->he = -1;
    half_edge->next_he = -1;
    half_edge->opposite_he = -1;
    half_edge->to_vertex = -1;
    half_edge->previous_he = -1;
}

qh_half_edge_t* qh__next_edge(qh_context_t* context)
{
    qh_half_edge_t* edge = context->edges + context->nedges;

    QH_ASSERT(context->nedges + 1 < context->maxedges);

    qh__init_half_edge(edge);

    edge->he = context->nedges;
    context->nedges++;

    QH_ASSERT(context->nedges < context->maxedges);

    return edge;
}

qh_face_t* qh__next_face(qh_context_t* context)
{
    qh_face_t* face = context->faces + context->nfaces;

    QH_ASSERT(context->nfaces + 1 < context->maxfaces);

    face->face = context->nfaces;
    face->iset.indices = NULL;
    context->nfaces++;

    QH_ASSERT(context->nedges < context->maxfaces);

    return face;
}

qh_vec3_t qh__edge_vec3(qh_half_edge_t* edge, qh_context_t* context)
{
    qh_half_edge_t prevhe = context->edges[edge->previous_he];
    qh_vec3_t v0, v1;

    v0 = context->vertices[prevhe.to_vertex];
    v1 = context->vertices[edge->to_vertex];

    qh__vec3_sub(&v1, &v0);
    qh__vec3_normalize(&v1);

    return v1;
}

void qh__face_init(qh_face_t* face,
    qh_index_t vertices[3],
    qh_context_t* context)
{
    qh_half_edge_t* e0 = qh__next_edge(context);
    qh_half_edge_t* e1 = qh__next_edge(context);
    qh_half_edge_t* e2 = qh__next_edge(context);
    qh_vec3_t v0, v1;
    qh_vertex_t centroid, normal;

    e2->to_vertex = vertices[0];
    e0->to_vertex = vertices[1];
    e1->to_vertex = vertices[2];

    e0->next_he = e1->he;
    e2->previous_he = e1->he;
    face->edges[1] = e1->he;

    e1->next_he = e2->he;
    e0->previous_he = e2->he;
    face->edges[2] = e2->he;
    v1 = qh__edge_vec3(e2, context);

    e2->next_he = e0->he;
    e1->previous_he = e0->he;
    face->edges[0] = e0->he;
    v0 = qh__edge_vec3(e0, context);

    e2->adjacent_face = face->face;
    e1->adjacent_face = face->face;
    e0->adjacent_face = face->face;

    qh__vec3_multiply(&v1, -1.f);
    normal = qh__vec3_cross(&v0, &v1);

    qh__vec3_normalize(&normal);
    centroid = qh__face_centroid(vertices, context);

    face->centroid = centroid;
    face->sdist = qh__vec3_dot(&normal, &centroid);
    face->normal = normal;
    face->iset.indices = QH_MALLOC(qh_index_t, QH_VERTEX_SET_SIZE);
    face->iset.capacity = QH_VERTEX_SET_SIZE;
    face->iset.size = 0;
    face->visitededges = 0;
    face->valid = 1;
}

void qh__tetrahedron_basis(qh_context_t* context, qh_index_t vertices[3])
{
    qh_index_t eps[6];
    int i, j, k, l;
    float max = -QH_FLT_MAX;

    qh__find_6eps(context->vertices, context->nvertices, eps);
    qh__find_2dps_6eps(context->vertices, eps, &j, &k);

    for (i = 0; i < 6; ++i) {
        float d2;

        if (i == j || i == k) {
            continue;
        }

        d2 = qh__vertex_segment_length2(context->vertices + eps[i],
            context->vertices + eps[j],
            context->vertices + eps[k]);

        if (d2 > max) {
            max = d2;
            l = i;
        }
    }

    vertices[0] = eps[j];
    vertices[1] = eps[k];
    vertices[2] = eps[l];
}

void qh__push_stack(qh_index_stack_t* stack, qh_index_t index)
{
    stack->begin[stack->size] = index;
    stack->size++;
}

qh_index_t qh__pop_stack(qh_index_stack_t* stack)
{
    qh_index_t top = -1;

    if (stack->size > 0) {
        top = stack->begin[stack->size - 1];
        stack->size--;
    }

    return top;
}

qh_index_t qh__furthest_point_from_plane(qh_context_t* context,
    qh_index_t* indices,
    int nindices,
    qh_vec3_t* normal,
    float sdist)
{
    int i, j;
    float max = -QH_FLT_MAX;

    for (i = 0; i < nindices; ++i) {
        qh_index_t index = indices ? *(indices + i) : i;
        float dist = qh__dist_point_plane(context->vertices + index, normal, sdist);

        if (dist > max) {
            j = i;
            max = dist;
        }
    }

    return j;
}

int qh__face_can_see_vertex(qh_face_t* face, qh_vertex_t* v) {
    qh_vec3_t tov = *v;

    qh__vec3_sub(&tov, &face->centroid);
    return qh__vec3_dot(&tov, &face->normal) > 0;
}

static inline void qh__assert_half_edge(qh_half_edge_t* edge, qh_context_t* context) {
    QH_ASSERT(edge->opposite_he != -1);
    QH_ASSERT(edge->he != -1);
    QH_ASSERT(edge->adjacent_face != -1);
    QH_ASSERT(edge->next_he != -1);
    QH_ASSERT(edge->previous_he != -1);
    QH_ASSERT(edge->to_vertex != -1);
    QH_ASSERT(context->edges[edge->opposite_he].to_vertex != edge->to_vertex);
}

static inline void qh__assert_face(qh_face_t* face, qh_context_t* context)
{
    int i;

    for (i = 0; i < 3; ++i) {
        qh__assert_half_edge(context->edges + face->edges[i], context);
    }

    QH_ASSERT(face->valid);
}

void qh__build_hull(qh_context_t* context)
{
    qh_index_t topface = qh__pop_stack(&context->facestack);
    int i, j, k;

    while (topface != -1) {
        qh_face_t* face = context->faces + topface;
        qh_index_t fvi, apex;
        qh_vertex_t* fv;
        int reversed = 0;

        if (!face->valid) {
            topface = qh__pop_stack(&context->facestack);
            continue;
        }

        fvi = qh__furthest_point_from_plane(context, face->iset.indices,
            face->iset.size, &face->normal, face->sdist);
        fv = context->vertices + *(face->iset.indices + fvi);

        qh__assert_face(face, context);

        // Reset visited flag for faces
        {
            for (i = 0; i < context->nfaces; ++i) {
                context->faces[i].visitededges = 0;
            }
        }

        // Find horizon edge
        {
            qh_index_t tovisit = topface;
            qh_face_t* facetovisit = context->faces + tovisit;

            // Release scratch
            context->scratch.size = 0;

            while (tovisit != -1) {
                if (facetovisit->visitededges >= 3) {
                    tovisit = qh__pop_stack(&context->scratch);
                    facetovisit = context->faces + tovisit;
                } else {
                    qh_index_t edgeindex = facetovisit->edges[facetovisit->visitededges];
                    qh_half_edge_t* edge;
                    qh_half_edge_t* oppedge;
                    qh_face_t* adjface;

                    facetovisit->visitededges++;

                    edge = context->edges + edgeindex;
                    oppedge = context->edges + edge->opposite_he;
                    adjface = context->faces + oppedge->adjacent_face;

                    qh__assert_half_edge(oppedge, context);
                    qh__assert_half_edge(edge, context);
                    qh__assert_face(adjface, context);

                    if (!qh__face_can_see_vertex(adjface, fv)) {
                        qh__push_stack(&context->horizonedges, edge->he);
                    } else {
                        facetovisit->valid = 0;
                        qh__push_stack(&context->scratch, adjface->face);
                    }
                }
            }
        }

        apex = face->iset.indices[fvi];

        // Sort horizon edges in CCW order
        {
            qh_vertex_t triangle[3];
            int vindex = 0;
            qh_vec3_t v0, v1, toapex;
            qh_vertex_t n;

            for (i = 0; i < context->horizonedges.size; ++i) {
                qh_index_t he0 = context->horizonedges.begin[i];
                qh_index_t he0vert = context->edges[he0].to_vertex;
                qh_index_t phe0 = context->edges[he0].previous_he;
                qh_index_t phe0vert = context->edges[phe0].to_vertex;

                for (j = i + 2; j < context->horizonedges.size; ++j) {
                    qh_index_t he1 = context->horizonedges.begin[j];
                    qh_index_t he1vert = context->edges[he1].to_vertex;
                    qh_index_t phe1 = context->edges[he1].previous_he;
                    qh_index_t phe1vert = context->edges[phe1].to_vertex;

                    if (phe1vert == he0vert || phe0vert == he1vert) {
                        QH_SWAP(qh_index_t, context->horizonedges.begin[j],
                                context->horizonedges.begin[i + 1]);
                        break;
                    }
                }

                if (vindex < 3) {
                    triangle[vindex++] = context->vertices[context->edges[he0].to_vertex];
                }
            }

            // Detect first triangle face ordering
            v0 = triangle[0];
            v1 = triangle[2];

            qh__vec3_sub(&v0, &triangle[1]);
            qh__vec3_sub(&v1, &triangle[1]);

            n = qh__vec3_cross(&v0, &v1);

            // Get the vector to the apex
            toapex = triangle[0];
            qh__vec3_sub(&toapex, context->vertices + apex);

            reversed = qh__vec3_dot(&n, &toapex) < 0.f;
        }

        // Create new faces
        {
            qh_index_t top = qh__pop_stack(&context->horizonedges);
            qh_index_t last = qh__pop_stack(&context->horizonedges);
            qh_index_t first = top;
            int looped = 0;

            QH_ASSERT(context->newhorizonedges.size == 0);

            // Release scratch
            context->scratch.size = 0;

            while (!looped) {
                qh_half_edge_t* prevhe;
                qh_half_edge_t* nexthe;
                qh_half_edge_t* oppedge;
                qh_index_t verts[3];
                qh_face_t* newface;

                if (last == -1) {
                    looped = 1;
                    last = first;
                }

                prevhe = context->edges + last;
                nexthe = context->edges + top;

                if (reversed) {
                    QH_SWAP(qh_half_edge_t*, prevhe, nexthe);
                }

                verts[0] = prevhe->to_vertex;
                verts[1] = nexthe->to_vertex;
                verts[2] = apex;

                context->faces[nexthe->adjacent_face].valid = 0;

                oppedge = context->edges + nexthe->opposite_he;
                newface = qh__next_face(context);

                qh__face_init(newface, verts, context);

                oppedge->opposite_he = context->edges[newface->edges[0]].he;
                context->edges[newface->edges[0]].opposite_he = oppedge->he;

                qh__push_stack(&context->scratch, newface->face);
                qh__push_stack(&context->newhorizonedges, newface->edges[0]);

                top = last;
                last = qh__pop_stack(&context->horizonedges);
            }
        }

        // Attach point sets to newly created faces
        {
            for (k = 0; k < context->nfaces; ++k) {
                qh_face_t* f = context->faces + k;

                if (f->valid || f->iset.size == 0) {
                    continue;
                }

                QH_ASSERT(f->visitededges <= 3);

                for (i = 0; i < f->iset.size; ++i) {
                    qh_index_t vertex = f->iset.indices[i];
                    qh_vertex_t* v = context->vertices + vertex;
                    qh_face_t* dface = NULL;

                    for (j = 0; j < context->scratch.size; ++j) {
                        qh_face_t* newface = context->faces + context->scratch.begin[j];
                        qh_half_edge_t* e0 = context->edges + newface->edges[0];
                        qh_half_edge_t* e1 = context->edges + newface->edges[1];
                        qh_half_edge_t* e2 = context->edges + newface->edges[2];
                        qh_vertex_t cv;

                        if (e0->to_vertex == vertex ||
                            e1->to_vertex == vertex ||
                            e2->to_vertex == vertex) {
                            continue;
                        }

                        cv = *v;
                        qh__vec3_sub(&cv, &newface->centroid);

                        if (qh__vec3_dot(&cv, &newface->normal) >= 0) {
                            dface = newface;
                            break;
                        }
                    }

                    if (dface) {
                        if (dface->iset.size + 1 >= dface->iset.capacity) {
                            dface->iset.capacity *= 2;
                            dface->iset.indices = QH_REALLOC(qh_index_t,
                                dface->iset.indices, dface->iset.capacity);
                        }

                        dface->iset.indices[dface->iset.size++] = vertex;
                    }
                }

                f->iset.size = 0;
            }
        }

        // Link new faces together
        {
            for (i = 0; i < context->newhorizonedges.size; ++i) {
                qh_index_t phe0, nhe1;
                qh_half_edge_t* he0;
                qh_half_edge_t* he1;
                int ii;

                if (reversed) {
                    ii = (i == 0) ? context->newhorizonedges.size - 1 : (i-1);
                } else {
                    ii = (i+1) % context->newhorizonedges.size;
                }

                phe0 = context->edges[context->newhorizonedges.begin[i]].previous_he;
                nhe1 = context->edges[context->newhorizonedges.begin[ii]].next_he;

                he0 = context->edges + phe0;
                he1 = context->edges + nhe1;

                QH_ASSERT(he1->to_vertex == apex);
                QH_ASSERT(he0->opposite_he == -1);
                QH_ASSERT(he1->opposite_he == -1);

                he0->opposite_he = he1->he;
                he1->opposite_he = he0->he;
            }

            context->newhorizonedges.size = 0;
        }

        // Push new face to stack
        {
            for (i = 0; i < context->scratch.size; ++i) {
                qh_face_t* face = context->faces + context->scratch.begin[i];

                if (face->iset.size > 0) {
                    qh__push_stack(&context->facestack, face->face);
                }
            }

            // Release scratch
            context->scratch.size = 0;
        }

        topface = qh__pop_stack(&context->facestack);
    }
}

qh_face_t* qh__build_tetrahedron(qh_context_t* context)
{
    int i, j;
    qh_index_t vertices[3];
    qh_index_t apex;
    qh_face_t* faces;
    qh_vertex_t normal, centroid, vapex;

    // Get the initial tetrahedron basis (first face)
    qh__tetrahedron_basis(context, &vertices[0]);

    // Find apex from the tetrahedron basis
    {
        float sdist;
        qh_vec3_t v0, v1;

        v0 = context->vertices[vertices[1]];
        v1 = context->vertices[vertices[2]];

        qh__vec3_sub(&v0, context->vertices + vertices[0]);
        qh__vec3_sub(&v1, context->vertices + vertices[0]);

        normal = qh__vec3_cross(&v0, &v1);
        qh__vec3_normalize(&normal);

        centroid = qh__face_centroid(vertices, context);
        sdist = qh__vec3_dot(&normal, &centroid);

        apex = qh__furthest_point_from_plane(context, NULL,
            context->nvertices, &normal, sdist);
        vapex = context->vertices[apex];

        qh__vec3_sub(&vapex, &centroid);

        // Whether the face is looking towards the apex
        if (qh__vec3_dot(&vapex, &normal) > 0) {
            QH_SWAP(qh_index_t, vertices[1], vertices[2]);
        }
    }

    faces = qh__next_face(context);
    qh__face_init(&faces[0], vertices, context);

    // Build faces from the tetrahedron basis to the apex
    {
        qh_index_t facevertices[3];
        for (i = 0; i < 3; ++i) {
            qh_half_edge_t* edge = context->edges + faces[0].edges[i];
            qh_half_edge_t prevedge = context->edges[edge->previous_he];
            qh_face_t* face = faces+i+1;
            qh_half_edge_t* e0;

            facevertices[0] = edge->to_vertex;
            facevertices[1] = prevedge.to_vertex;
            facevertices[2] = apex;

            qh__next_face(context);
            qh__face_init(face, facevertices, context);

            e0 = context->edges + faces[i+1].edges[0];
            edge->opposite_he = e0->he;
            e0->opposite_he = edge->he;
        }
    }

    // Attach half edges to faces tied to the apex
    {
        for (i = 0; i < 3; ++i) {
            qh_face_t* face;
            qh_face_t* nextface;
            qh_half_edge_t* e1;
            qh_half_edge_t* e2;

            j = (i+2) % 3;

            face = faces+i+1;
            nextface = faces+j+1;

            e1 = context->edges + face->edges[1];
            e2 = context->edges + nextface->edges[2];

            QH_ASSERT(e1->opposite_he == -1);
            QH_ASSERT(e2->opposite_he == -1);

            e1->opposite_he = e2->he;
            e2->opposite_he = e1->he;

            qh__assert_half_edge(e1, context);
            qh__assert_half_edge(e2, context);
        }
    }

    // Create initial point set; every point is
    // attached to the first face it can see
    {
        for (i = 0; i < context->nvertices; ++i) {
            qh_vertex_t* v;
            qh_face_t* dface = NULL;

            if (vertices[0] == i || vertices[1] == i || vertices[2] == i) {
                continue;
            }

            v = context->vertices+i;

            for (j = 0; j < 4; ++j) {
                qh_face_t* face = context->faces+j;
                qh_vertex_t cv = *v;
                qh__vec3_sub(&cv, &face->centroid);

                if (qh__vec3_dot(&cv, &face->normal) >= 0) {
                    dface = face;
                }
            }

            if (dface) {
                if (dface->iset.size + 1 >= dface->iset.capacity) {
                    dface->iset.capacity *= 2;
                    dface->iset.indices = QH_REALLOC(qh_index_t,
                        dface->iset.indices, dface->iset.capacity);
                }

                dface->iset.indices[dface->iset.size++] = i;
            }
        }
    }

    // Add initial tetrahedron faces to the face stack
    for (i = 0; i < 4; ++i) {
        context->faces[i].valid = 1;
        qh__assert_face(context->faces + i, context);
        if (context->faces[i].iset.size != 0) {
            qh__push_stack(&context->facestack, i);
        }
    }

    QH_ASSERT(context->nedges == context->nfaces * 3);
    QH_ASSERT(context->nfaces == 4);

    return faces;
}

void qh__init_context(qh_context_t* context, qh_vertex_t* vertices, unsigned int nvertices)
{
    // TODO:
    // size_t nedges = 3 * nvertices - 6;
    // size_t nfaces = 2 * nvertices - 4;
    size_t nfaces = nvertices * (nvertices - 1);
    size_t nedges = nfaces * 3;

    context->edges = QH_MALLOC(qh_half_edge_t, nedges);
    context->faces = QH_MALLOC(qh_face_t, nfaces);
    context->facestack.begin = QH_MALLOC(qh_index_t, nfaces);
    context->scratch.begin = QH_MALLOC(qh_index_t, nfaces);
    context->horizonedges.begin = QH_MALLOC(qh_index_t, nedges);
    context->newhorizonedges.begin = QH_MALLOC(qh_index_t, nedges);

    context->vertices = vertices;
    context->nvertices = nvertices;
    context->nedges = 0;
    context->nfaces = 0;
    context->facestack.size = 0;
    context->scratch.size = 0;
    context->horizonedges.size = 0;
    context->newhorizonedges.size = 0;
}

void qh__free_context(qh_context_t* context)
{
    int i;

    for (i = 0; i < context->nfaces; ++i) {
        QH_FREE(context->faces[i].iset.indices);
        context->faces[i].iset.size = 0;
    }

    context->nvertices = 0;
    context->nfaces = 0;

    QH_FREE(context->edges);

    QH_FREE(context->faces);
    QH_FREE(context->facestack.begin);
    QH_FREE(context->scratch.begin);
    QH_FREE(context->horizonedges.begin);
    QH_FREE(context->newhorizonedges.begin);
}

void qh_free_mesh(qh_mesh_t mesh)
{
    QH_FREE(mesh.vertices);
    QH_FREE(mesh.indices);
    QH_FREE(mesh.normalindices);
}

qh_mesh_t qh_quickhull3d(qh_vertex_t* vertices, unsigned int nvertices)
{
    qh_mesh_t m;
    qh_context_t context;
    unsigned int* indices;
    unsigned int nfaces = 0, i, index, nindices;

    qh__init_context(&context, vertices, nvertices);

    // Build the initial tetrahedron
    qh__build_tetrahedron(&context);

    // Build the convex hull
    qh__build_hull(&context);

    for (i = 0; i < context.nfaces; ++i) {
        if (context.faces[i].valid) { nfaces++; }
    }

    nindices = nfaces * 3;

    indices = QH_MALLOC(unsigned int, nindices);

    m.normals = QH_MALLOC(qh_vertex_t, nfaces);
    m.nnormals = nfaces;
    m.normalindices = QH_MALLOC(unsigned int, nfaces);

    // Assign normals
    {
        index = 0;

        for (i = 0; i < context.nfaces; ++i) {
            if (!context.faces[i].valid) { continue; }

            m.normals[index] = context.faces[i].normal;
            index++;
        }
    }

    // Assign normal indices and default set of indices
    {
        index = 0;

        for (i = 0; i < context.nfaces; ++i) {
            if (!context.faces[i].valid) { continue; }

            m.normalindices[index] = index;

            qh_half_edge_t e0 = context.edges[context.faces[i].edges[0]];
            qh_half_edge_t e1 = context.edges[context.faces[i].edges[1]];
            qh_half_edge_t e2 = context.edges[context.faces[i].edges[2]];

            indices[index*3+0] = e0.to_vertex;
            indices[index*3+1] = e1.to_vertex;
            indices[index*3+2] = e2.to_vertex;

            index++;
        }
    }

    // Rearrange vertices to remove unused ones
    {
        int* ni = QH_MALLOC(int, context.nfaces * 3);

        m.vertices = QH_MALLOC(qh_vertex_t, nindices);
        m.indices = QH_MALLOC(unsigned int, nindices);
        m.nindices = nindices;
        m.nvertices = 0;

        for (i = 0; i < context.nfaces * 3; ++i) {
            ni[i] = -1;
        }

        for (i = 0; i < nfaces * 3; ++i) {
            unsigned int oi = indices[i];

            if (ni[oi] == -1) {
                m.vertices[m.nvertices++] = context.vertices[oi];
                ni[oi] = m.nvertices - 1;
            }

            m.indices[i] = ni[oi];
        }

        QH_FREE(ni);
    }

    QH_FREE(indices);

    qh__free_context(&context);

    return m;
}

#endif // QUICKHULL_IMPLEMENTATION
