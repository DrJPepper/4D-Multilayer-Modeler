#pragma once

#include "util.h"

class Scene {
    public:
        Scene();
        vector<vtkSphereSource *> spheres;
        vector<vtkLineSource *> lines;
        vector<vtkActor *> sphereActors;
        vector<vtkActor *> lineActors;
};
