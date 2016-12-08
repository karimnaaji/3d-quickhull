# 3d-quickhull
_Header only 3d quickhull in ANSI C_

![](images/quickhull.png)

Usage
=====

To use this library, simply include `quickhull.h` once with the `QUICKHULL_IMPLEMENTATION` define in a `.cpp` file.

```c
#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"
```

The usage of the library is quite simple, generate or gather a set of points, and call `qh_quickhull3d`. The result is a mesh with a set of indexed normals and vertices ready to upload in a GPU.

```c
const int n = 100;
qh_vertex_t vertices[n];

for (int i = 0; i < n; ++i) {
    float a0 = (rand_0_1() * M_PI * 2);
    float a1 = (rand_0_1() * M_PI * 2);
    vertices[i].z = sin(a0) * radius;
    vertices[i].x = cos(a1) * cos(a0) * rand_0_1();
    vertices[i].y = sin(a1) * cos(a0) * rand_0_1();
}

qh_mesh_t mesh = qh_quickhull3d(vertices, n);

// ...

qh_free_mesh(mesh);

```

Example
=======


- Open Frameworks wrapper: https://github.com/satcy/ofx3DQuickHull
- Obj to hull: https://github.com/karimnaaji/3d-quickhull/blob/master/example/command_line.cpp

If you're interested in low-polygon rendering, using quickhull as a base for mesh triangulation can give such results:

![](images/mesh_quickhull.png)
![](http://karim.naaji.fr/images/lowpoly0.png)
