#ifndef QUICKHULL_DEBUG_H
#define QUICKHULL_DEBUG_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"

static void printf_vertex(qh_vertex_t* v) {
    printf("%f %f %f\n", v->x, v->y, v->z);
}

static float rand_0_1() {
    return (float) rand() / RAND_MAX;
}

static qh_vertex_t rand_vertex(float scale) {
    qh_vertex_t v;
    v.x = (2.f * rand_0_1() - 1.f) * scale;
    v.y = (2.f * rand_0_1() - 1.f) * scale;
    v.z = (2.f * rand_0_1() - 1.f) * scale;
    return v;
}

#endif
